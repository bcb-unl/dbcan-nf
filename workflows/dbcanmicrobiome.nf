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

// new subworkflows
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS as FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS as FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                as FASTQC_TRIMGALORE_DNA         } from '../subworkflows/local/fastqc_trimgalore'
include { FASTQC_TRIMGALORE                as FASTQC_TRIMGALORE_RNA         } from '../subworkflows/local/fastqc_trimgalore'
include { BWAMEME_INDEX_MEM                as BWAMEME_INDEX_MEM_DNA          } from '../subworkflows/local/bwameme_index_mem'
include { BWAMEME_INDEX_MEM                as BWAMEME_INDEX_MEM_RNA          } from '../subworkflows/local/bwameme_index_mem'
include { CGC_DEPTH                        as CGC_DEPTH_DNA            } from '../subworkflows/local/cgc_depth'
include { CGC_DEPTH                        as CGC_DEPTH_RNA            } from '../subworkflows/local/cgc_depth'

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
                def dna_meta = [ id: meta.id + '_dna', single_end: meta.single_end ]
                tuple(dna_meta, fq_list)
            }
            .filter { meta, fq_list -> fq_list && fq_list.size() > 0 }

        // RNA channel
        ch_samplesheet_rna = ch_samplesheet
            .map { meta, fastqs, transcriptomes ->
                def t_list = transcriptomes ?: []
                def rna_meta = [ id: meta.id + '_rna' ]
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

        CGC_DEPTH_DNA(
            RUNDBCAN_EASYSUBSTRATE.out.dbcan_results,
            ch_bam_bai_dna
        )

        ch_versions = ch_versions.mix(CGC_DEPTH_DNA.out.versions)
        ch_cgc_depth_dna = CGC_DEPTH_DNA.out.tsv

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
if (ch_samplesheet_rna) {
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
    extract_kraken2_reads_rna = FASTQC_TRIMGALORE_RNA.out.reads
    // 3. combine all reads for the kraken2
    //subworkflow to extract kraken2 reads
    if (!params.skip_kraken_extraction) {
        // TODO: build db process
        // Use nf-core subworkflow at https://nf-co.re/subworkflows/fasta_build_add_kraken2/
        FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA (
            extract_kraken2_reads_rna,
            ch_db_for_kraken2,
            params.kraken_tax
        )
        extract_kraken2_reads_rna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.extracted_kraken2_reads

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.multiqc_files)
        ch_versions = ch_versions.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.versions)
    }

    BWAMEME_INDEX_MEM_RNA (
        ch_megahit_contigs,
        extract_kraken2_reads_rna,
    )
    ch_bam_bai_rna = BWAMEME_INDEX_MEM_RNA.out.ch_bam_bai
    ch_versions = ch_versions.mix(BWAMEME_INDEX_MEM_RNA.out.versions)


        ch_gff_bam_bai_rna = ch_gunzip_gff
            .join(ch_bam_bai_rna, by: 0)
            .map { tuple ->
            def (meta, gff, bam_tuple) = tuple
            def (meta2, bam, bai) = bam_tuple
            tuple(meta, gff, bam, bai)
        }

    RUNDBCAN_UTILS_CAL_COVERAGE_RNA(
        ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, gff) },
        ch_gff_bam_bai_rna.map { meta, gff, bam, bai -> tuple(meta, bam, bai) }
    )

    ch_versions = ch_versions.mix(RUNDBCAN_UTILS_CAL_COVERAGE_RNA.out.versions)
    ch_coverage_rna = RUNDBCAN_UTILS_CAL_COVERAGE_RNA.out.depth_txt

    CGC_DEPTH_RNA (
        RUNDBCAN_EASYSUBSTRATE.out.dbcan_results,
        ch_bam_bai_rna
    )
    ch_versions = ch_versions.mix(CGC_DEPTH_RNA.out.versions)
    ch_cgc_depth_rna = CGC_DEPTH_RNA.out.tsv

    ch_abund_input_rna = ch_coverage_rna
        .join(ch_dbcan_results, by: 0)
        .map { meta, cgc_depth, dbcan_results -> tuple(meta, cgc_depth, dbcan_results) }

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
}






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
