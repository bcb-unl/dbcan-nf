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
include { SEQTK_SEQ                      } from '../modules/nf-core/seqtk/seq/main'
include { DIAMOND_BLASTX                  } from '../modules/nf-core/diamond/blastx/main'
include { RUNDBCAN_ASMFREE_CAL_ABUND      } from '../modules/local/dbcan_asmfree_cal_abund'
include { DOWNLOAD_CAZY_DMND              } from '../modules/local/download_cazy_dmnd/main'
include { RUNDBCAN_PLOT_BAR            as RUNDBCAN_PLOT_BAR_DNA        } from '../modules/local/dbcan_plot'
// new subworkflows
include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS as FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                as FASTQC_TRIMGALORE_DNA         } from '../subworkflows/local/fastqc_trimgalore'

//prepare the project parameters and databases

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow DBCANMICROBIOMEASSEMBLYFREE {
    take:
        ch_samplesheet

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        // DNA channel
        // Note: For assemfree type, samplesheet format is tuple(meta, fastqs) - only 2 elements
        ch_samplesheet_dna = ch_samplesheet
            .map { meta, fastqs ->
                def fq_list = fastqs ?: []
                def dna_meta = [ id: meta.id+'_dna', single_end: meta.single_end ]
                tuple(dna_meta, fq_list)
            }
            .filter { meta, fq_list -> fq_list && fq_list.size() > 0 }

        // Download CAZy DIAMOND database for assembly-free analysis
        // For assembly-free, we only need CAZy.dmnd, so download it directly instead of full dbCAN database
        DOWNLOAD_CAZY_DMND()
        ch_cazy_dmnd = DOWNLOAD_CAZY_DMND.out.cazy_dmnd
        ch_cazyid_subfam_mapping = DOWNLOAD_CAZY_DMND.out.cazyid_subfam_mapping
        ch_versions = ch_versions.mix(DOWNLOAD_CAZY_DMND.out.versions)

        // Prepare diamond database channel for DIAMOND_BLASTX
        // DIAMOND_BLASTX expects tuple(meta, db) where db can be a file or directory
        ch_diamond_db = ch_cazy_dmnd
            .map { cazy_dmnd_file ->
                // Pass the directory containing the .dmnd file
                // The module script will find the .dmnd file in this directory
                tuple([ id: 'diamond_db' ], cazy_dmnd_file.parent)
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

        // 1. FastQC + TrimGalore
        FASTQC_TRIMGALORE_DNA (
            ch_samplesheet_dna,
            params.skip_fastqc,
            params.skip_trimming
        )

        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQC_TRIMGALORE_DNA.out.fastqc_zip.map{ it[1][0] }
        )

        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE_DNA.out.versions)
        ch_trimmed_reads_dna = FASTQC_TRIMGALORE_DNA.out.reads

        // 2. Kraken2 extraction (optional)
        if (!params.skip_kraken_extraction) {
            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA (
                    ch_trimmed_reads_dna,
                    ch_db_for_kraken2,
                    params.kraken_tax)

            ch_extracted_from_kraken2_reads_dna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.extracted_kraken2_reads

            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.multiqc_files)
            ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_DNA.out.versions )
            ch_assembly_free_input_dna = ch_extracted_from_kraken2_reads_dna
        } else {
            // if skip kraken2 extraction, use trimmed reads
            ch_assembly_free_input_dna = ch_trimmed_reads_dna
        }

        // 3. Convert FASTQ to FASTA using seqtk seq
        // Note: ext.args '-A' should be set in modules.config to convert FASTQ to FASTA
        // TRIMGALORE output format: tuple(meta, path) where path can match multiple files
        // When path matches multiple files, Nextflow returns them as a list
        // We need to expand this to separate tuples for R1 and R2 for paired-end reads
        ch_fasta_for_blastx_dna = ch_assembly_free_input_dna
            .flatMap { meta, reads_list ->
                // TRIMGALORE output: reads_list is a list when multiple files match
                // Convert to list if it's not already a list
                def reads = reads_list instanceof List ? reads_list : [reads_list]
                reads = reads.findAll { it != null && it != false }
                
                if (reads.size() > 1) {
                    // Paired-end: process R1 and R2 separately
                    def r1_meta = meta.clone()
                    r1_meta.id = "${meta.id}_R1"
                    def r2_meta = meta.clone()
                    r2_meta.id = "${meta.id}_R2"
                    [
                        tuple(r1_meta, reads[0]),
                        tuple(r2_meta, reads[1])
                    ]
                } else if (reads.size() == 1) {
                    // Single-end: process once
                    [tuple(meta, reads[0])]
                } else {
                    // No valid reads, return empty list
                    []
                }
            }

        // Convert to FASTA
        SEQTK_SEQ(ch_fasta_for_blastx_dna)
        ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)
        ch_fasta_dna = SEQTK_SEQ.out.fastx

        // 4. Diamond blastx: align reads to CAZyme database
        // For paired-end, we process R1 and R2 separately to get two different blastx results
        ch_fasta_for_blastx_dna_final = ch_fasta_dna
        ch_db_for_blastx_dna = ch_diamond_db

        DIAMOND_BLASTX(
            ch_fasta_for_blastx_dna_final,
            ch_db_for_blastx_dna,
            'txt',
            'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen '
        )

        ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)
        ch_diamond_results_dna = DIAMOND_BLASTX.out.txt

        // 5. dbcan_asmfree: calculate abundance (all types in one module)
        // Prepare reads for dbcan_asmfree (need original reads, not fasta)
        ch_reads_for_asmfree_dna = ch_assembly_free_input_dna
            .map { meta, reads_list ->
                // For paired-end, use first read file (R1) as --raw_reads
                // reads_list is a list of files for paired-end, or a single file for single-end
                def reads_file = reads_list instanceof List ? reads_list[0] : reads_list
                tuple(meta, reads_file)
            }

        // Group diamond results by sample (for paired-end, we have R1 and R2 results)
        // Extract base sample ID (remove _R1 or _R2 suffix)
        ch_diamond_results_grouped = ch_diamond_results_dna
            .map { meta, diamond_result ->
                def base_id = meta.id.replaceFirst(/_R[12]$/, '')
                tuple(base_id, meta, diamond_result)
            }
            .groupTuple()
            .map { base_id, metas, results ->
                // Determine if paired-end based on number of results
                def is_paired = results.size() > 1
                def original_meta = metas[0].clone()
                original_meta.id = base_id
                
                if (is_paired) {
                    // Paired-end: return tuple with both R1 and R2 results
                    // Keep base_id as first element for join operation
                    tuple(base_id, original_meta, results[0], results[1])
                } else {
                    // Single-end: return tuple with single result (duplicate for compatibility)
                    tuple(base_id, original_meta, results[0], results[0])
                }
            }

        // Join diamond results with reads
        ch_reads_for_asmfree_keyed = ch_reads_for_asmfree_dna
            .map { meta, reads_file ->
                tuple(meta.id, meta, reads_file)
            }

        ch_asmfree_input_dna = ch_diamond_results_grouped
            .join(ch_reads_for_asmfree_keyed, by: 0)
            .map { key, meta_diamond, diamond_r1, diamond_r2, meta_reads, reads_file ->
                tuple(meta_diamond, diamond_r1, diamond_r2, meta_reads, reads_file)
            }

        // Prepare inputs for dbcan_asmfree
        // ch_asmfree_input_dna already contains: tuple(meta_diamond, diamond_r1, diamond_r2, meta_reads, reads_file)
        // Extract diamond and reads channels separately for module inputs
        ch_diamond_for_asmfree_dna = ch_asmfree_input_dna.map { meta, diamond_r1, diamond_r2, reads_meta, reads_file -> tuple(meta, diamond_r1, diamond_r2) }
        ch_reads_for_asmfree_dna_final = ch_asmfree_input_dna.map { meta, diamond_r1, diamond_r2, reads_meta, reads_file -> tuple(meta, reads_file) }
        
        // Prepare db directory for dbcan_asmfree (contains CAZyID_subfam_mapping.tsv)
        // Get the parent directory from the mapping file
        // DOWNLOAD_CAZY_DMND only runs once, so ch_cazyid_subfam_mapping has only one value
        // Use .first() to convert to value channel
        ch_db_dir_for_asmfree = ch_cazyid_subfam_mapping.map { it.parent }.first()
        
        // Calculate all abundance types (fam, subfam, EC, substrate) in one module
        // This is similar to RUNDBCAN_UTILS_CAL_ABUND in the assembly-based workflow
        RUNDBCAN_ASMFREE_CAL_ABUND(
            ch_diamond_for_asmfree_dna,
            ch_reads_for_asmfree_dna_final,
            ch_db_dir_for_asmfree,
            'FPKM'
        )

        ch_versions = ch_versions.mix(RUNDBCAN_ASMFREE_CAL_ABUND.out.versions)
        ch_asmfree_abund_dna = RUNDBCAN_ASMFREE_CAL_ABUND.out.abund_dir

        // 6. Plotting (similar to assembly-based workflow)
        ch_plot_input_dna = ch_asmfree_abund_dna
            .toList()
            .map { tuples ->
                tuple(
                    tuples*.getAt(0).id,
                    tuples*.getAt(1)
                )
            }

        RUNDBCAN_PLOT_BAR_DNA(ch_plot_input_dna)
        ch_versions = ch_versions.mix(RUNDBCAN_PLOT_BAR_DNA.out.versions)

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
