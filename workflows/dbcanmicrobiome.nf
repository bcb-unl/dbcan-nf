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
include { BWAMEME_INDEX_MEM                as BWAMEME_INDEX_MEM_DNA          } from '../subworkflows/local/bwameme_index_mem'
include { BWAMEME_INDEX_MEM                as BWAMEME_INDEX_MEM_RNA          } from '../subworkflows/local/bwameme_index_mem'
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
                def dna_meta = [ id: meta.id, type:'dna', single_end: meta.single_end ]
                tuple(dna_meta, fq_list)
            }
            .filter { meta, fq_list -> fq_list && fq_list.size() > 0 }

        // RNA channel
        ch_samplesheet_rna = ch_samplesheet
            .map { meta, fastqs, transcriptomes ->
                def t_list = transcriptomes ?: []
                def rna_meta = [ id: meta.id, type:'rna', single_end: meta.single_end ]
                tuple(rna_meta, t_list)
            }
            .filter { meta, t_list -> t_list && t_list.size() > 0 }
            .ifEmpty { Channel.empty() }


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

        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQC_TRIMGALORE_DNA.out.fastqc_zip.map{ it[1][0] }
        )

        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE_DNA.out.versions)
        ch_trimmed_reads = FASTQC_TRIMGALORE_DNA.out.reads

        // 3. combine all reads for the kraken2
        //subworkflow to extract kraken2 reads

        if (!params.skip_kraken_extraction) {
            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA (
                    ch_trimmed_reads,
                    ch_db_for_kraken2,
                    params.kraken_tax)
            ch_extracted_from_kraken2_reads = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.extracted_kraken2_reads

            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.multiqc_files)
            ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.versions )
            ch_megahit_input_dna = ch_extracted_from_kraken2_reads
            //ch_megahit_input_dna.view()
        // if skip kraken2 extraction, use trimmed reads
        } else {
            ch_megahit_input_dna = ch_trimmed_reads
            .map { meta, reads_list ->
            if (reads_list.size() == 2) {
                tuple(meta, reads_list[0], reads_list[1])
            } else {
                tuple(meta, reads_list[0], null)
            }
                }
                }
        // MODULE: Megahit to assemble metagenomics
        // make sure all reads are the list
    if (params.subsample && !params.coassembly) {
        // Subsample each read (single/paired)
        ch_subsample_input = ch_megahit_input_dna
        .flatMap { meta, read1, read2 ->
            def reads = read2 ? [read1, read2] : [read1]
            reads.collect { tuple(meta, it, params.subsample_size as Integer) }
        }

        SEQTK_SAMPLE(ch_subsample_input)

        ch_subsampled = SEQTK_SAMPLE.out.reads
            .map { meta, reads -> tuple(meta, reads) }
            .groupTuple(by: 0)
            .map { meta, reads_list ->
                if (reads_list.size() == 2) {
                    tuple(meta, reads_list[0], reads_list[1])
                } else {
                    tuple(meta, reads_list[0], null)
                }
            }
        ch_megahit_input_final = ch_subsampled

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

        // [meta, r1, r2, meta, r1, r2, ...]
        ch_coassembly_input = ch_megahit_input_dna
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

        COMBINE_PAIRED_READS(ch_coassembly_input)
        ch_megahit_input_final = COMBINE_PAIRED_READS.out.reads
        //ch_megahit_input_final.view()

    } else {
        // Use original reads
        ch_megahit_input_final = ch_megahit_input_dna
    }

        // pair-end for MEGAHIT
        MEGAHIT ( ch_megahit_input_dna )
        ch_versions = ch_versions.mix ( MEGAHIT.out.versions )
        ch_megahit_contigs = MEGAHIT.out.contigs

        //
        // MODULE: Pyrodigal to find genes in metagenomics data
        // Will make it as subworkflow
        PYRODIGAL ( ch_megahit_contigs, 'gff' )
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

        //
        // MODULE: run_dbCAN to find CAZymes and CGCs in metagenomics data
        //Will write it into subworkflow later

        ch_gunzip_gff_with_type = ch_gunzip_gff.map { meta, gff -> tuple(meta, gff, 'prodigal') }

        RUNDBCAN_EASYSUBSTRATE (
            ch_gunzip_faa,
            ch_gunzip_gff_with_type,
            ch_dbcan_db_final
        )

        ch_versions = ch_versions.mix(RUNDBCAN_DATABASE.out.versions)
        ch_dbcan_results = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results
        ch_bwameme_input_dna = ch_megahit_input_dna
            .map { meta, read1, read2 ->
                if (read2) {
                    tuple(meta, [read1, read2])
                } else {
                    tuple(meta, [read1])
                }
            }


        //ch_bwameme_input_dna.view()
        //
        // Read mapping with BWA-MEME
        BWAMEME_INDEX_MEM_DNA (
            ch_megahit_contigs,
            ch_bwameme_input_dna
        )
        //
        ch_bam_bai_dna = BWAMEME_INDEX_MEM_DNA.out.ch_bam_bai
        ch_versions = ch_versions.mix(BWAMEME_INDEX_MEM_DNA.out.versions)

        ch_gff_bam_bai_dna = ch_gunzip_gff
            .join(ch_bam_bai_dna, by: 0)
            .map { meta, gff, bam, bai ->
                tuple(meta, gff, bam, bai)
            }

        RUNDBCAN_UTILS_CAL_COVERAGE_DNA(
            ch_gff_bam_bai_dna.map { meta, gff, bam, bai -> tuple(meta, gff) },
            ch_gff_bam_bai_dna.map { meta, gff, bam, bai -> tuple(meta, bam, bai) }
        )


        ch_versions = ch_versions.mix(RUNDBCAN_UTILS_CAL_COVERAGE_DNA.out.versions)
        ch_coverage_dna = RUNDBCAN_UTILS_CAL_COVERAGE_DNA.out.depth_txt

        CGC_DEPTH_PLOT_DNA(
            RUNDBCAN_EASYSUBSTRATE.out.dbcan_results,
            ch_bam_bai_dna,
            ch_dbcan_db_final
        )

        ch_versions = ch_versions.mix(CGC_DEPTH_PLOT_DNA.out.versions)
        ch_cgc_depth_dna = CGC_DEPTH_PLOT_DNA.out.tsv

        ch_abund_input_dna = ch_coverage_dna
            .join(ch_dbcan_results, by: 0)
            .map { meta, cgc_depth, dbcan_results -> tuple(meta, cgc_depth, dbcan_results) }

        RUNDBCAN_UTILS_CAL_ABUND_DNA(
            ch_abund_input_dna
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
        //ch_plot_input.view()

        //
        RUNDBCAN_PLOT_BAR_DNA(ch_plot_input_dna)
        //ch_samplesheet_rna.view()
// only if ch_samplesheet_rna has content
    // rna processing
    FASTQC_TRIMGALORE_RNA (
        ch_samplesheet_rna,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_TRIMGALORE_RNA.out.fastqc_zip.map{ it[1][0] }
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE_RNA.out.versions)
    ch_trimmed_reads_rna = FASTQC_TRIMGALORE_RNA.out.reads

    // 3. combine all reads for the kraken2
    //subworkflow to extract kraken2 reads
        if (!params.skip_kraken_extraction) {
            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA (
                    ch_trimmed_reads_rna,
                    ch_db_for_kraken2,
                    params.kraken_tax)
            ch_extracted_from_kraken2_reads_rna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.extracted_kraken2_reads

            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.multiqc_files)
            ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.versions )
            extract_kraken2_reads_rna = ch_extracted_from_kraken2_reads_rna
            //ch_megahit_input_dna.view()
        // if skip kraken2 extraction, use trimmed reads
        } else {
            extract_kraken2_reads_rna = ch_trimmed_reads_rna
        }

    BWAMEME_INDEX_MEM_RNA (
        ch_megahit_contigs,
        extract_kraken2_reads_rna,
    )
    ch_bam_bai_rna = BWAMEME_INDEX_MEM_RNA.out.ch_bam_bai
    ch_versions = ch_versions.mix(BWAMEME_INDEX_MEM_RNA.out.versions)

    ch_gff_bam_bai_rna = ch_gunzip_gff
            .map { meta, gff -> tuple(meta.id, meta, gff) }
            .join(
                ch_bam_bai_rna.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) },
                by: 0
            )
            .map { id, meta_gff, gff, meta_bam, bam, bai ->
                tuple(meta_gff, gff, bam, bai)
            }

    RUNDBCAN_UTILS_CAL_COVERAGE_RNA(
            ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, gff) },
            ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, bam, bai) }
        )


    ch_versions = ch_versions.mix(RUNDBCAN_UTILS_CAL_COVERAGE_RNA.out.versions)
    ch_coverage_rna = RUNDBCAN_UTILS_CAL_COVERAGE_RNA.out.depth_txt

    ch_dbcan_results_rna = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results
        .map { meta, dbcan_results -> tuple(meta, dbcan_results) }

    CGC_DEPTH_PLOT_RNA (
        ch_dbcan_results_rna,
        ch_bam_bai_rna,
        ch_dbcan_db_final
    )
    ch_versions = ch_versions.mix(CGC_DEPTH_PLOT_RNA.out.versions)
    ch_cgc_depth_rna = CGC_DEPTH_PLOT_RNA.out.tsv

    ch_abund_input_rna = ch_coverage_rna
        .map { meta, cgc_depth -> tuple(meta.id, meta, cgc_depth) }
        .join(
            ch_dbcan_results.map { meta, dbcan_results -> tuple(meta.id, meta, dbcan_results) },
            by: 0
        )
        .map { id, meta_cov, cgc_depth, meta_dbcan, dbcan_results ->
            tuple(meta_cov, cgc_depth, dbcan_results)
        }

    RUNDBCAN_UTILS_CAL_ABUND_RNA(
        ch_abund_input_rna
    )

    ch_plot_input_rna = RUNDBCAN_UTILS_CAL_ABUND_RNA.out.abund_dir
        .toList()
        .map { tuples ->
            tuple(
                tuples*.getAt(0).id,
                tuples*.getAt(1)
            )
        }

    //RUNDBCAN_PLOT_BAR_RNA(ch_plot_input_rna)







    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
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
