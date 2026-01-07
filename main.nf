#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/dbcanmicrobiome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/dbcanmicrobiome
    Website: https://nf-co.re/dbcanmicrobiome
    Slack  : https://nfcore.slack.com/channels/dbcanmicrobiome
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DBCANMICROBIOME  } from './workflows/dbcanmicrobiome'
include { DBCANMICROBIOMELONG } from './workflows/dbcanmicrobiome_long'
include { DBCANMICROBIOMEASSEMBLYFREE } from './workflows/dbcanmicrobiome_assembly_free'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_dbcanmicrobiome_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_dbcanmicrobiome_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_DBCANMICROBIOME {
    take:
        samplesheet

    main:
        if (params.type == 'shortreads') {
            DBCANMICROBIOME(samplesheet)
            multiqc_report = DBCANMICROBIOME.out.multiqc_report
        } else if (params.type == 'longreads') {
            DBCANMICROBIOMELONG(samplesheet)
            multiqc_report = DBCANMICROBIOMELONG.out.multiqc_report
        } else if (params.type == 'assemfree') {
            DBCANMICROBIOMEASSEMBLYFREE(samplesheet)
            multiqc_report = DBCANMICROBIOMEASSEMBLYFREE.out.multiqc_report
        } else {
            log.error "Invalid value for --type parameter. Please choose from 'shortreads', 'longreads', or 'assemfree'."
            exit 1
        }

    emit:
        multiqc_report = multiqc_report
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_DBCANMICROBIOME (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_DBCANMICROBIOME.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
