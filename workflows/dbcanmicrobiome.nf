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

// new subworkflows
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                } from '../subworkflows/local/fastqc_trimgalore'


//prepare the project parameters and databases

ch_kraken2_db_file = params.kraken_db ? path(params.kraken_db, type: 'dir', checkIfExists: true) : []
ch_kraken2_tax     = params.kraken_tax ?: []
ch_dbcan_db        = params.dbcan_db ? path(params.dbcan_db, checkIfExists: true) : []

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
    ch_samplesheet_pe = ch_samplesheet.filter { meta, reads -> !meta.single_end }
    ch_samplesheet_se = ch_samplesheet.filter { meta, reads -> meta.single_end }

    // 2. pair-end process
    FASTQC_TRIMGALORE (
        ch_samplesheet_pe,
        params.skip_fastqc,
        params.skip_trimming
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_TRIMGALORE.out.fastqc_zip.map{ it[1][0] }
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    extract_kraken2_reads_pe = FASTQC_TRIMGALORE.out.reads
    //3. single end process
    extract_kraken2_reads_se = ch_samplesheet_se

    // 4. combine all reads for the kraken2
    extract_kraken2_reads = extract_kraken2_reads_pe.mix(extract_kraken2_reads_se)

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
    extract_kraken2_reads_fixed = extract_kraken2_reads.map { meta, reads ->
        tuple(meta, reads instanceof List ? reads : [reads])
    }

    // pair-end for MEGAHIT
    ch_megahit_input = extract_kraken2_reads_fixed
        .filter { meta, reads -> !meta.single_end && reads.size() == 2 && reads[0] && reads[1] }
        .map { meta, reads -> tuple(meta, reads[0], reads[1]) }

    // single-end for FLYE
    ch_flye_input = extract_kraken2_reads_fixed
        .filter { meta, reads -> meta.single_end && reads.size() == 1 && reads[0] }
        .map { meta, reads -> tuple(meta, reads[0]) }
    ch_flye_mode  = Channel.value(params.flye_mode ?: "--pacbio-hifi")

    MEGAHIT ( ch_megahit_input )
    ch_flye_input.view { "FLYE input: $it" }
    FLYE ( ch_flye_input, ch_flye_mode )
    FLYE.out.fasta.view { "FLYE contigs: $it" }
    ch_versions = ch_versions.mix ( MEGAHIT.out.versions )
    ch_megahit_contigs = MEGAHIT.out.contigs
    ch_flye_contigs = FLYE.out.fasta

    // combine all contigs
    ch_all_contigs = ch_megahit_contigs.mix(ch_flye_contigs)

    //
    // MODULE: Pyrodigal to find genes in metagenomics data
    // Will make it as subworkflow
    PYRODIGAL ( ch_all_contigs, 'gff' )
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
