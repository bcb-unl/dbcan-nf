/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dbcanmicrobiome_pipeline'

// new modules in the pipeline
include { KRAKEN2_BUILDSTANDARD           } from '../modules/nf-core/kraken2/buildstandard/main'
include { KRAKEN2_KRAKEN2                 } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS  } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { TRIMGALORE                      } from '../modules/nf-core/trimgalore/main'
include { MEGAHIT                         } from '../modules/nf-core/megahit/main'
include { PYRODIGAL                       } from '../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_FAA            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF            } from '../modules/nf-core/gunzip/main'
include { RUNDBCAN_DATABASE               } from '../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_EASYSUBSTRATE          } from '../modules/nf-core/rundbcan/easysubstrate/main'
include { RUNDBCAN_UTILS_CAL_COVERAGE  as RUNDBCAN_UTILS_CAL_COVERAGE_DNA  } from '../modules/local/dbcan_utils_cal_coverage'
include { RUNDBCAN_UTILS_CAL_COVERAGE  as RUNDBCAN_UTILS_CAL_COVERAGE_RNA  } from '../modules/local/dbcan_utils_cal_coverage'
include { BWA_INDEX                     } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX               as SAMTOOLS_INDEX_DNA         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX               as SAMTOOLS_INDEX_RNA         } from '../modules/nf-core/samtools/index/main'
include { RUNDBCAN_UTILS_CAL_ABUND     as RUNDBCAN_UTILS_CAL_ABUND_DNA    } from '../modules/local/dbcan_utils_cal_abund'
include { RUNDBCAN_UTILS_CAL_ABUND     as RUNDBCAN_UTILS_CAL_ABUND_RNA    } from '../modules/local/dbcan_utils_cal_abund'
include { RUNDBCAN_PLOT_BAR            as RUNDBCAN_PLOT_BAR_DNA        } from '../modules/local/dbcan_plot'
include { RUNDBCAN_PLOT_BAR            as RUNDBCAN_PLOT_BAR_RNA        } from '../modules/local/dbcan_plot'
include { SEQTK_SAMPLE                 } from '../modules/nf-core/seqtk/sample/main'
include { COMBINE_PAIRED_READS        } from '../modules/local/combine_raw_reads_before_coassembly'
// new subworkflows
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS as FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS as FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                as FASTQC_TRIMGALORE_DNA         } from '../subworkflows/local/fastqc_trimgalore'
include { FASTQC_TRIMGALORE                as FASTQC_TRIMGALORE_RNA         } from '../subworkflows/local/fastqc_trimgalore'
include { BWA_INDEX_MEM                as BWA_INDEX_MEM_DNA          } from '../subworkflows/local/bwa_index_mem'
include { BWA_INDEX_MEM                as BWA_INDEX_MEM_RNA          } from '../subworkflows/local/bwa_index_mem'
include { CGC_DEPTH_PLOT                   as CGC_DEPTH_PLOT_DNA            } from '../subworkflows/local/cgc_depth_plot'
include { CGC_DEPTH_PLOT                   as CGC_DEPTH_PLOT_RNA            } from '../subworkflows/local/cgc_depth_plot'

//prepare the project parameters and databases

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow DBCANMICROBIOME {
    take:
        ch_samplesheet

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        // DNA channel
        ch_samplesheet_dna = ch_samplesheet
            .map { meta, fastqs, transcriptomes ->
                def fq_list = fastqs ?: []
                def dna_meta = [ id: meta.id+'_dna', single_end: meta.single_end ]
                tuple(dna_meta, fq_list)
            }
            .filter { meta, fq_list -> fq_list && fq_list.size() > 0 }

        // RNA channel, only use RNA samples when not subsample or coassembly
        // Store flag to avoid repeated condition checks
        def use_rna = !params.subsample && !params.coassembly
        
        if (use_rna) {
            ch_samplesheet_rna = ch_samplesheet
                .map { meta, fastqs, transcriptomes ->
                    def t_list = transcriptomes ?: []
                    def rna_meta = [ id: meta.id+'_rna', single_end: meta.single_end ]
                    tuple(rna_meta, t_list)
                }
                .filter { meta, t_list -> t_list && t_list.size() > 0 }
        } else {
            ch_samplesheet_rna = Channel.empty()
        }


        //tmp view channels
        //ch_samplesheet_dna.view()
        //ch_samplesheet_rna.view()



        // process database parameters
        if (params.dbcan_db) {
            ch_dbcan_db_final = Channel.fromPath(params.dbcan_db, checkIfExists: true)
        } else {
            RUNDBCAN_DATABASE()
            ch_dbcan_db_final = RUNDBCAN_DATABASE.out.dbcan_db
            ch_versions = ch_versions.mix(RUNDBCAN_DATABASE.out.versions)
        }

        // Process kraken database
        if (!params.skip_kraken_extraction) {
            if (params.kraken_db) {
                ch_db_for_kraken2 = Channel.fromPath(params.kraken_db, type: 'dir', checkIfExists: true)
            } else {
                KRAKEN2_BUILDSTANDARD(false)
                ch_db_for_kraken2 = KRAKEN2_BUILDSTANDARD.out.db
            }
        }


        // 2. pair-end process
        FASTQC_TRIMGALORE_DNA (
            ch_samplesheet_dna,
            params.skip_fastqc,
            params.skip_trimming
        )

        // Process RNA channels if enabled
        if (use_rna) {
            FASTQC_TRIMGALORE_RNA (
                ch_samplesheet_rna,
                params.skip_fastqc,
                params.skip_trimming
            )
            ch_trimmed_reads_rna = FASTQC_TRIMGALORE_RNA.out.reads
        } else {
            ch_trimmed_reads_rna = Channel.empty()
        }

        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQC_TRIMGALORE_DNA.out.fastqc_zip.map{ it[1][0] }
        )

        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE_DNA.out.versions)
        ch_trimmed_reads_dna = FASTQC_TRIMGALORE_DNA.out.reads

        // 3. combine all reads for the kraken2
        //subworkflow to extract kraken2 reads

        if (!params.skip_kraken_extraction) {
            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA (
                    ch_trimmed_reads_dna,
                    ch_db_for_kraken2,
                    params.kraken_tax)

            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA (
                    ch_trimmed_reads_rna,
                    ch_db_for_kraken2,
                    params.kraken_tax)

            ch_extracted_from_kraken2_reads_dna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.extracted_kraken2_reads
            ch_extracted_from_kraken2_reads_rna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.extracted_kraken2_reads

            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.multiqc_files)
            ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.versions )
            ch_megahit_input_dna = ch_extracted_from_kraken2_reads_dna
            //ch_megahit_input_dna.view()
        // if skip kraken2 extraction, use trimmed reads
        } else {
            ch_megahit_input_dna = ch_trimmed_reads_dna
            .map { meta, reads_list ->
            if (reads_list.size() == 2) {
                tuple(meta, reads_list[0], reads_list[1])
            } else {
                tuple(meta, reads_list[0], null)
            }
                }

            ch_extracted_from_kraken2_reads_rna = ch_trimmed_reads_rna
            .map { meta, reads_list ->
                if (reads_list.size() == 2) {
                    tuple(meta, reads_list[0], reads_list[1])
                } else {
                    tuple(meta, reads_list[0], null)
                }
            }
        }
            //ch_extracted_from_kraken2_reads_rna.view()
        // MODULE: Megahit to assemble metagenomics
        // make sure all reads are the list
    if (params.subsample && !params.coassembly) {
        // Subsample each read (single/paired)
        ch_subsample_input_dna = ch_megahit_input_dna
        .flatMap { meta, read1, read2 ->
            def reads = read2 ? [read1, read2] : [read1]
            reads.collect { tuple(meta, it, params.subsample_size as Integer) }
        }

        SEQTK_SAMPLE(ch_subsample_input_dna)

    ch_subsampled_dna = SEQTK_SAMPLE.out.reads
        .map { meta, reads -> tuple(meta, reads) }
        .groupTuple(by: 0)
        .map { meta, reads_list ->
            def subsample_meta = meta.clone()
            subsample_meta.id = meta.id + '_subsample'
            def sorted_reads = reads_list.sort { it.toString() }
            if (sorted_reads.size() == 2) {
                tuple(subsample_meta, sorted_reads[0], sorted_reads[1])
            } else {
                tuple(subsample_meta, sorted_reads[0], null)
            }
        }
        ch_megahit_input_final_dna = ch_subsampled_dna

    } else if (params.coassembly) {
        // check at least 2 samples
        ch_megahit_input_dna.count().subscribe { n ->
            if (n < 2) {
                error "Co-assembly mode requires at least 2 samples, but only ${n} sample(s) are present. Please do not use the --coassembly parameter with a single sample."
            }
        }
        //ch_megahit_input_dna.view()
        // combine all reads for co-assembly
        // This snippet prepares input for co-assembly in a DSL2-compliant way.
        // Input: ch_megahit_input_dna (tuple(meta, read1, read2))
        // Output: ch_coassembly_input (tuple(meta, all_read1, all_read2 OR null))

        // Save original samples for later read mapping
        // Ensure meta.id has _dna suffix for consistency
        ch_original_samples_dna = ch_megahit_input_dna
            .map { meta, r1, r2 ->
                def new_meta = meta.clone()
                if (!new_meta.id.endsWith('_dna')) {
                    new_meta.id = new_meta.id + '_dna'
                }
                tuple(new_meta, r1, r2)
            }

        // [meta, r1, r2, meta, r1, r2, ...]
        ch_coassembly_input_dna = ch_megahit_input_dna
            .collect()
            .map { flat_list ->
                def tuples = flat_list.collate(3) //
                def all_read1 = tuples.collect { it[1] }
                def all_read2 = tuples.collect { it[2] }.findAll { it }
                def meta = [ id: 'coassembly', single_end: all_read2.isEmpty() ]
                if (meta.single_end) {
                    tuple(meta, all_read1, null)
                } else {
                    tuple(meta, all_read1, all_read2)
                }
            }

        COMBINE_PAIRED_READS(ch_coassembly_input_dna)
        ch_megahit_input_final_dna = COMBINE_PAIRED_READS.out.reads
        .map { meta, reads ->
            if (meta.single_end) {
                tuple(meta, reads[0], null)
            } else {
                tuple(meta, reads[0], reads[1])
            }
        }
    } else {
        // Use original reads
        ch_megahit_input_final_dna = ch_megahit_input_dna
        // For non-coassembly mode, set ch_original_samples_dna to same as input
        ch_original_samples_dna = ch_megahit_input_dna
    }

        // pair-end for MEGAHIT
        MEGAHIT ( ch_megahit_input_final_dna )
        ch_versions = ch_versions.mix ( MEGAHIT.out.versions )
        ch_megahit_contigs_dna = MEGAHIT.out.contigs

        //
        // MODULE: Pyrodigal to find genes in metagenomics data
        // Will make it as subworkflow
        PYRODIGAL ( ch_megahit_contigs_dna, 'gff' )
        ch_versions = ch_versions.mix(PYRODIGAL.out.versions)
        ch_faa = PYRODIGAL.out.faa
        ch_gff = PYRODIGAL.out.annotations

        //
        // MODULE: GUNZIP to unzip faa and gff files
        // Will make it as subworkflow from pyrodigal to gunzip

        GUNZIP_FAA (ch_faa)
        GUNZIP_GFF (ch_gff)
        ch_gunzip_faa = GUNZIP_FAA.out.gunzip
        ch_gunzip_gff = GUNZIP_GFF.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP_FAA.out.versions)
        
        // Save original coassembly gff and faa for RUNDBCAN_EASYSUBSTRATE in coassembly mode
        if (params.coassembly) {
            ch_gunzip_gff_coassembly = ch_gunzip_gff
            ch_gunzip_faa_coassembly = ch_gunzip_faa
        } else {
            ch_gunzip_gff_coassembly = ch_gunzip_gff
            ch_gunzip_faa_coassembly = ch_gunzip_faa
        }

        //
        // MODULE: run_dbCAN to find CAZymes and CGCs in metagenomics data
        //Will write it into subworkflow later
        // In coassembly mode, run on coassembly results first, then replicate to samples

        ch_gunzip_gff_with_type = ch_gunzip_gff_coassembly.map { meta, gff -> tuple(meta, gff, 'prodigal') }

        RUNDBCAN_EASYSUBSTRATE (
            ch_gunzip_faa_coassembly,
            ch_gunzip_gff_with_type,
            ch_dbcan_db_final
        )

        ch_versions = ch_versions.mix(RUNDBCAN_EASYSUBSTRATE.out.versions)
        ch_dbcan_results = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results
        
        // In coassembly mode, replicate coassembly results to each original sample
        if (params.coassembly) {
            // In coassembly mode, there should be only one dbcan result (from coassembly)
            // Use combine (not cross) to replicate single result to each original sample
            // ch_original_samples_dna already has _dna suffix, so we can use it directly
            ch_dbcan_results = ch_dbcan_results
                .combine(ch_original_samples_dna.map { meta, r1, r2 -> meta })
                .map { coassembly_meta, dbcan_results, original_meta ->
                    // Use original sample meta (already has _dna suffix)
                    tuple(original_meta, dbcan_results)
                }
            
            // Also replicate gff and faa to each original sample for downstream processing
            ch_gunzip_gff = ch_gunzip_gff
                .combine(ch_original_samples_dna.map { meta, r1, r2 -> meta })
                .map { coassembly_meta, gff, original_meta ->
                    // Use original sample meta (already has _dna suffix)
                    tuple(original_meta, gff)
                }
            
            ch_gunzip_faa = ch_gunzip_faa
                .combine(ch_original_samples_dna.map { meta, r1, r2 -> meta })
                .map { coassembly_meta, faa, original_meta ->
                    // Use original sample meta (already has _dna suffix)
                    tuple(original_meta, faa)
                }
        }
        
        ch_gunzip_faa_rna = ch_gunzip_faa
            .map { meta, faa ->
                def new_meta = meta.clone()
                if (new_meta.id.endsWith('_dna')) {
                    new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
                }
                tuple(new_meta, faa)
            }

        ch_gunzip_gff_rna = ch_gunzip_gff
            .map { meta, gff ->
                def new_meta = meta.clone()
                if (new_meta.id.endsWith('_dna')) {
                    new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
                }
                tuple(new_meta, gff)
            }
        
        ch_dbcan_results_rna = ch_dbcan_results
            .map { meta, dbcan_results ->
                def new_meta = meta.clone()
                if (new_meta.id.endsWith('_dna')) {
                    new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
                }
                tuple(new_meta, dbcan_results)
            }
        //ch_megahit_input_final_dna.view()

        // In coassembly mode, use original samples for read mapping
        // ch_original_samples_dna already has _dna suffix, so we can use it directly
        if (params.coassembly) {
            ch_bwameme_input_dna = ch_original_samples_dna
        } else {
            ch_bwameme_input_dna = ch_megahit_input_final_dna
        }
        //ch_bwameme_input_dna.view()
        ch_bwameme_input_rna = ch_extracted_from_kraken2_reads_rna
        //ch_bwameme_input_rna.view()

        //
        // In coassembly mode, create index once using coassembly contigs, then replicate to each sample
        // For non-coassembly, use contigs directly
        if (params.coassembly) {
            // In coassembly mode, there should be only one contig set (from coassembly)
            // Create index once for coassembly contigs
            BWA_INDEX(ch_megahit_contigs_dna)
            
            // In coassembly mode, replicate index to each sample with sample meta.id
            // BWA_MEM uses meta.id to find index directory (bwa_${meta.id}), so we need
            // to use sample meta.id for all channels to match, but the actual index
            // directory is from coassembly. We'll create symbolic links in the process.
            ch_index_dna_coassembly = BWA_INDEX.out.index
            
            // Use combine (not cross) to replicate single coassembly result to each original sample
            // combine creates all combinations, cross requires matching keys
            // ch_original_samples_dna already has _dna suffix, so we can use it directly
            ch_index_dna = ch_index_dna_coassembly
                .combine(ch_original_samples_dna.map { meta, r1, r2 -> meta })
                .map { coassembly_meta, index_dir, original_meta ->
                    // Use sample meta for join and index lookup (already has _dna suffix)
                    def sample_meta = original_meta.clone()
                    // Store coassembly meta.id in a custom field for later use
                    sample_meta.coassembly_index_id = coassembly_meta.id
                    tuple(sample_meta, index_dir)
                }
            
            // Replicate contigs to each original sample for join operation
            ch_megahit_contigs_dna_for_join = ch_megahit_contigs_dna
                .combine(ch_original_samples_dna.map { meta, r1, r2 -> meta })
                .map { contigs_meta, contigs, original_meta ->
                    // Use original sample meta (already has _dna suffix)
                    tuple(original_meta, contigs)
                }
        } else {
            BWA_INDEX(ch_megahit_contigs_dna)
            ch_index_dna = BWA_INDEX.out.index
            ch_megahit_contigs_dna_for_join = ch_megahit_contigs_dna
        }

        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

        //ch_bwameme_input_dna.view()
        //
        // Read mapping with BWA-MEME
        // Use meta.id as join key to ensure proper matching
        ch_bwameme_input_dna_keyed = ch_bwameme_input_dna.map { meta, r1, r2 ->
            tuple(meta.id, meta, r1, r2)
        }
        
        // All channels must use the same meta.id for matching in BWA_INDEX_MEM
        ch_index_dna_keyed_for_join = ch_index_dna.map { meta, idx ->
            tuple(meta.id, meta, idx)
        }
        
        ch_contigs_dna_keyed_for_join = ch_megahit_contigs_dna_for_join.map { meta, fasta ->
            tuple(meta.id, meta, fasta)
        }
        
        ch_bwa_mem_input = ch_bwameme_input_dna_keyed
            .join(ch_index_dna_keyed_for_join, by: 0)
            .join(ch_contigs_dna_keyed_for_join, by: 0)
            .map { key, meta, r1, r2, meta_idx, idx, meta_ctg, fasta ->
                tuple(meta, r1, r2, meta_idx, idx, fasta)
            }


        ch_fastq = ch_bwa_mem_input.map { t -> tuple(t[0], t[1], t[2]) } // (meta, read1, read2)
        ch_index = ch_bwa_mem_input.map { t -> tuple(t[3], t[4]) }      // (meta, index) - meta.id matches sample, but index dir is coassembly
        ch_genome_fna = ch_bwa_mem_input.map { t -> tuple(t[0], t[5]) } // (meta, fasta)

        BWA_INDEX_MEM_DNA(
            ch_index,
            ch_genome_fna,
            ch_fastq
        )

        //
        // Prepare channels for RNA mapping (reuse DNA index prefix)
        //

        // Create join key based on sample name (remove _rna/_dna suffix)
        ch_bwameme_input_rna_keyed = ch_bwameme_input_rna.map { meta, r1, r2 ->
            def key = meta.id.replaceFirst(/_rna$/, '')
            tuple(key, meta, r1, r2)
        }

        ch_index_dna_keyed = ch_index_dna.map { meta, idx ->
            def key = meta.id.replaceFirst(/_dna$/, '')
            tuple(key, meta, idx)
        }

        // In coassembly mode, use ch_megahit_contigs_dna_for_join which has been replicated to each sample
        // In non-coassembly mode, use ch_megahit_contigs_dna directly
        ch_contigs_dna_keyed = ch_megahit_contigs_dna_for_join.map { meta, fasta ->
            def key = meta.id.replaceFirst(/_dna$/, '')
            tuple(key, meta, fasta)
        }

        // Join RNA reads with DNA index/contigs using base key
        ch_bwa_mem_input_rna = ch_bwameme_input_rna_keyed
            .join(ch_index_dna_keyed, by: 0)
            .join(ch_contigs_dna_keyed, by: 0)
            .map { key, meta_rna, r1, r2, meta_idx, idx, meta_ctg, fasta ->
                // Use DNA meta (with *_dna prefix) for BWA to match the generated index prefix
                tuple(meta_idx, r1, r2, idx, fasta)
            }

        // Split into subworkflow inputs
        ch_fastq_rna       = ch_bwa_mem_input_rna.map { t -> tuple(t[0], t[1], t[2]) }
        ch_index_rna_final = ch_bwa_mem_input_rna.map { t -> tuple(t[0], t[3]) }
        ch_genome_fna_rna  = ch_bwa_mem_input_rna.map { t -> tuple(t[0], t[4]) }

        // Call RNA subworkflow (will use *_dna prefix index internally)
        BWA_INDEX_MEM_RNA (
            ch_index_rna_final,
            ch_genome_fna_rna,
            ch_fastq_rna
        )

        // Change output meta from *_dna back to *_rna for subsequent join
        ch_bam_bai_rna = BWA_INDEX_MEM_RNA.out.ch_bam_bai.map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
            tuple(new_meta, bam, bai)
        }

        //
        ch_bam_bai_dna = BWA_INDEX_MEM_DNA.out.ch_bam_bai

        // Use meta.id as join key to avoid meta object mismatch issues
        ch_gff_bam_bai_dna = ch_gunzip_gff
            .map { meta, gff -> tuple(meta.id, meta, gff) }
            .join(ch_bam_bai_dna.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }, by: 0)
            .map { key, meta_gff, gff, meta_bam, bam, bai ->
                tuple(meta_gff, gff, bam, bai)
            }

        ch_gff_bam_bai_rna = ch_gunzip_gff_rna
            .map { meta, gff -> tuple(meta.id, meta, gff) }
            .join(ch_bam_bai_rna.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }, by: 0)
            .map { key, meta_gff, gff, meta_bam, bam, bai ->
                tuple(meta_gff, gff, bam, bai)
            }

        RUNDBCAN_UTILS_CAL_COVERAGE_DNA(
            ch_gff_bam_bai_dna.map { meta, gff, bam, bai -> tuple(meta, gff) },
            ch_gff_bam_bai_dna.map { meta, gff, bam, bai -> tuple(meta, bam, bai) }
        )

        RUNDBCAN_UTILS_CAL_COVERAGE_RNA(
            ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, gff) },
            ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, bam, bai) }
        )

        ch_versions = ch_versions.mix(RUNDBCAN_UTILS_CAL_COVERAGE_DNA.out.versions)
        ch_coverage_dna = RUNDBCAN_UTILS_CAL_COVERAGE_DNA.out.depth_txt
        ch_coverage_rna = RUNDBCAN_UTILS_CAL_COVERAGE_RNA.out.depth_txt

        CGC_DEPTH_PLOT_DNA(
            ch_dbcan_results,
            ch_bam_bai_dna,
            ch_dbcan_db_final
        )

        CGC_DEPTH_PLOT_RNA(
            ch_dbcan_results_rna,
            ch_bam_bai_rna,
            ch_dbcan_db_final
        )
        ch_versions = ch_versions.mix(CGC_DEPTH_PLOT_DNA.out.versions)
        ch_cgc_depth_dna = CGC_DEPTH_PLOT_DNA.out.tsv
        ch_cgc_depth_rna = CGC_DEPTH_PLOT_RNA.out.tsv

        // Use meta.id as join key to avoid meta object mismatch issues
        ch_abund_input_dna = ch_coverage_dna
            .map { meta, depth -> tuple(meta.id, meta, depth) }
            .join(ch_dbcan_results.map { meta, results -> tuple(meta.id, meta, results) }, by: 0)
            .map { key, meta1, depth, meta2, results -> tuple(meta1, depth, results) }

        ch_abund_input_rna = ch_coverage_rna
            .map { meta, depth -> tuple(meta.id, meta, depth) }
            .join(ch_dbcan_results_rna.map { meta, results -> tuple(meta.id, meta, results) }, by: 0)
            .map { key, meta1, depth, meta2, results -> tuple(meta1, depth, results) }


        RUNDBCAN_UTILS_CAL_ABUND_DNA(
            ch_abund_input_dna
        )
        RUNDBCAN_UTILS_CAL_ABUND_RNA(
            ch_abund_input_rna
        )
        //RUNDBCAN_UTILS_CAL_ABUND_DNA.out.abund_dir.view()
        ch_plot_input_dna = RUNDBCAN_UTILS_CAL_ABUND_DNA.out.abund_dir
        .toList()
        .map { tuples ->
            tuple(
                tuples*.getAt(0).id,
                tuples*.getAt(1)
            )
        }


        ch_plot_input_rna = RUNDBCAN_UTILS_CAL_ABUND_RNA.out.abund_dir
        .toList()
        .map { tuples ->
            tuple(
                tuples*.getAt(0).id,
                tuples*.getAt(1)
            )
        }
        //ch_plot_input.view()
        //ch_plot_input_rna.view()
        //
        RUNDBCAN_PLOT_BAR_DNA(ch_plot_input_dna)
        RUNDBCAN_PLOT_BAR_RNA(ch_plot_input_rna)
        ch_versions = ch_versions.mix(RUNDBCAN_PLOT_BAR_DNA.out.versions)

    //
    // Collate and save software versions
    //
    // softwareVersionsToYAML only accepts paths; convert ch_versions to paths and filter empty
    ch_versions_paths = ch_versions
        .map { v ->
            v instanceof List && v.size() >= 2 ? v[1] : v
        }
        .filter { it != null }
        .map { file(it, checkIfExists: true) }

    softwareVersionsToYAML(ch_versions_paths)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'dbcanmicrobiome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }



    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
        multiqc_report = MULTIQC.out.report.toList()
        versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
