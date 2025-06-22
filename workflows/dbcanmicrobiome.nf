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

include { FLYE                            } from '../modules/nf-core/flye/main'
include { RUNDBCAN_UTILS_CAL_COVERAGE     } from '../modules/local/dbcan_utils_cal_coverage'
include { SAMTOOLS_INDEX                  } from '../modules/nf-core/samtools/index/main'
include { RUNDBCAN_UTILS_CAL_ABUND        } from '../modules/local/dbcan_utils_cal_abund'
include { RUNDBCAN_PLOT_BAR               } from '../modules/local/dbcan_plot'

// new subworkflows
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                } from '../subworkflows/local/fastqc_trimgalore'
include { BWAMEME_INDEX_MEM                } from '../subworkflows/local/bwameme_index_mem'
include { CGC_DEPTH                        } from '../subworkflows/local/cgc_depth'

//prepare the project parameters and databases

ch_kraken2_db_file = params.kraken_db ?
    Channel.fromPath(params.kraken_db, checkIfExists: true) :
    Channel.empty()
ch_dbcan_db = params.dbcan_db ? Channel.fromPath(params.dbcan_db, checkIfExists: true) : Channel.empty()
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DBCANMICROBIOME {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC_TrimGalore
    //

    // 1. split channels based on the type of reads
    // DNA channel
    ch_samplesheet
    .map { meta, fastqs, transcriptomes ->
        def fq_list = fastqs ?: []
        def dna_meta = [ id: meta.id, single_end: meta.single_end ]
        tuple(dna_meta, fq_list)
    }
    .filter { meta, fq_list -> fq_list && fq_list.size() > 0 }
    .set { ch_samplesheet_dna }

    // RNA channel
    ch_samplesheet
    .map { meta, fastqs, transcriptomes ->
        def t_list = transcriptomes ?: []
        def rna_meta = [ id: meta.id ]
        tuple(rna_meta, t_list)
    }
    .filter { meta, t_list -> t_list && t_list.size() > 0 }
    .set { ch_samplesheet_rna }


    //ch_samplesheet_dna.view()
    //ch_samplesheet_rna.view()

    // 2. pair-end process
    FASTQC_TRIMGALORE (
        ch_samplesheet_dna,
        params.skip_fastqc,
        params.skip_trimming
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_TRIMGALORE.out.fastqc_zip.map{ it[1][0] }
    )

    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    extract_kraken2_reads_pe = FASTQC_TRIMGALORE.out.reads

    // 3. combine all reads for the kraken2
    extract_kraken2_reads = extract_kraken2_reads_pe

    //
    // MODULE: Kraken2 Build Database
    //
    if (!params.skip_kraken_extraction) {
        if (ch_kraken2_db_file.isEmpty()) {
            KRAKEN2_BUILDSTANDARD (false)
            ch_db_for_kraken2 = KRAKEN2_BUILDSTANDARD.out.db
        } else {
            ch_db_for_kraken2 = ch_kraken2_db_file
        }
    } else {
        ch_db_for_kraken2 = ch_kraken2_db_file
    }

    //subworkflow to extract kraken2 reads

    if (!params.skip_kraken_extraction) {
        // TODO: build db process
        // Use nf-core subworkflow at https://nf-co.re/subworkflows/fasta_build_add_kraken2/

        FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS (
                extract_kraken2_reads,
                ch_db_for_kraken2,
                params.tax_id
        ).extract_kraken2_reads.set { extract_kraken2_reads }

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS.out.multiqc_files)
        ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS.out.versions )
    }


    // MODULE: Megahit to assemble metagenomics

    // make sure all reads are the list
    extract_kraken2_reads_fixed = extract_kraken2_reads

    // pair-end for MEGAHIT
    ch_megahit_input = extract_kraken2_reads_fixed
        .filter { meta, reads -> !meta.single_end && reads.size() == 2 && reads[0] && reads[1] }
        .map { meta, reads -> tuple(meta, reads[0], reads[1]) }


    MEGAHIT ( ch_megahit_input )

    ch_versions = ch_versions.mix ( MEGAHIT.out.versions )
    ch_megahit_contigs = MEGAHIT.out.contigs

    //
    // MODULE: Pyrodigal to find genes in metagenomics data
    // Will make it as subworkflow
    PYRODIGAL ( ch_megahit_contigs, 'gff' )
    ch_versions = ch_versions.mix(PYRODIGAL.out.versions)
    ch_faa = PYRODIGAL.out.faa
    ch_gff = PYRODIGAL.out.annotations

    GUNZIP_FAA (ch_faa)
    GUNZIP_GFF (ch_gff)
    ch_gunzip_faa = GUNZIP_FAA.out.gunzip
    ch_gunzip_gff = GUNZIP_GFF.out.gunzip


    //
    // MODULE: run_dbCAN to find CAZymes and CGCs in metagenomics data
    //Will write it into subworkflow later
    RUNDBCAN_DATABASE ()
    if (!ch_dbcan_db.isEmpty()) {
        ch_dbcan_db = Channel.fromPath(ch_dbcan_db, checkIfExists: true)
    } else {
        ch_dbcan_db = RUNDBCAN_DATABASE.out.dbcan_db
    }

    ch_gunzip_gff_with_type = ch_gunzip_gff.map { meta, gff -> tuple(meta, gff, 'prodigal') }

    RUNDBCAN_EASYSUBSTRATE (
        ch_gunzip_faa,
        ch_gunzip_gff_with_type,
        ch_dbcan_db
    )

    ch_versions = ch_versions.mix(RUNDBCAN_DATABASE.out.versions)
    ch_dbcan_results = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results

    //
    // Read mapping with BWA-MEME
    BWAMEME_INDEX_MEM (
        ch_megahit_contigs,
        extract_kraken2_reads_fixed,
    )
    //
    ch_bam_bai = BWAMEME_INDEX_MEM.out.ch_bam_bai
    ch_versions = ch_versions.mix(BWAMEME_INDEX_MEM.out.versions)



    RUNDBCAN_UTILS_CAL_COVERAGE (
        ch_gunzip_gff,
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(RUNDBCAN_UTILS_CAL_COVERAGE.out.versions)
    ch_coverage = RUNDBCAN_UTILS_CAL_COVERAGE.out.depth_txt


    CGC_DEPTH(
        RUNDBCAN_EASYSUBSTRATE.out.dbcan_results,
        ch_bam_bai
    )

    ch_versions = ch_versions.mix(CGC_DEPTH.out.versions)
    ch_cgc_depth = CGC_DEPTH.out.tsv

    ch_abund_input = ch_coverage
        .join(ch_dbcan_results, by: 0)
        .map { meta, cgc_depth, dbcan_results -> tuple(meta, cgc_depth, dbcan_results) }

    RUNDBCAN_UTILS_CAL_ABUND(
    ch_abund_input
    )

    RUNDBCAN_UTILS_CAL_ABUND.out.abund_dir.view()
    ch_plot_input = RUNDBCAN_UTILS_CAL_ABUND.out.abund_dir
    .toList()
    .map { tuples ->
        tuple(
            tuples*.getAt(0).id,
            tuples*.getAt(1)
        )
    }
    //ch_plot_input.view()

    //
    RUNDBCAN_PLOT_BAR(ch_plot_input)
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

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
