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
include { FLYE                            } from '../modules/nf-core/flye/main'
include { PYRODIGAL                       } from '../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_FAA            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF            } from '../modules/nf-core/gunzip/main'
include { RUNDBCAN_DATABASE               } from '../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_EASYSUBSTRATE          } from '../modules/nf-core/rundbcan/easysubstrate/main'
include { RUNDBCAN_UTILS_CAL_COVERAGE  as RUNDBCAN_UTILS_CAL_COVERAGE_DNA  } from '../modules/local/dbcan_utils_cal_coverage'
include { RUNDBCAN_UTILS_CAL_COVERAGE  as RUNDBCAN_UTILS_CAL_COVERAGE_RNA  } from '../modules/local/dbcan_utils_cal_coverage'
include { BWA_INDEX                     } from '../modules/nf-core/bwa/index/main'
include { MINIMAP2_ALIGN                } from '../modules/nf-core/minimap2/align/main'
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


workflow DBCANMICROBIOMELONG {
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

        // RNA channel - 只在提供了 transcriptomes 时创建通道,否则为空
        ch_samplesheet_rna = ch_samplesheet
            .map { meta, fastqs, transcriptomes ->
                def t_list = transcriptomes ?: []
                if (t_list && t_list.size() > 0) {
                    def rna_meta = [ id: meta.id+'_rna', single_end: meta.single_end ]
                    tuple(rna_meta, t_list)
                } else {
                    null
                }
            }
            .filter { it != null }


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
        // DNA: long reads 不做 QC/trim，直接透传
        ch_trimmed_reads_dna = ch_samplesheet_dna

        // RNA：仅当存在 RNA 样本时再做 QC/trim
        FASTQC_TRIMGALORE_RNA (
            ch_samplesheet_rna,
            params.skip_fastqc,
            params.skip_trimming
        )
        ch_trimmed_reads_rna = FASTQC_TRIMGALORE_RNA.out.reads
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE_RNA.out.versions)

        if (!params.skip_kraken_extraction) {
            FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA (
                    ch_trimmed_reads_rna,
                    ch_db_for_kraken2,
                    params.kraken_tax)

            ch_extracted_from_kraken2_reads_rna = FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.extracted_kraken2_reads
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.multiqc_files)
            ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS_RNA.out.versions )

        } else {
            ch_extracted_from_kraken2_reads_rna = ch_trimmed_reads_rna
                .map { meta, reads_list ->
                    if (reads_list.size() == 2) {
                        tuple(meta, reads_list[0], reads_list[1])
                    } else {
                        tuple(meta, reads_list[0], null)
                    }
                }
        }

        // FLYE 仅接受 (meta, reads)；长读取第一个文件
        ch_flye_input_dna = ch_trimmed_reads_dna
            .map { meta, reads_list ->
                tuple(meta, reads_list[0])
            }

        // Use original reads
        ch_flye_input_final_dna = ch_flye_input_dna


        // pair-end for FLYE
        FLYE ( ch_flye_input_final_dna, '--pacbio-hifi' )
        ch_versions = ch_versions.mix ( FLYE.out.versions )
        ch_flye_contigs_dna = FLYE.out.fasta

        //
        // MODULE: Pyrodigal to find genes in metagenomics data
        // Will make it as subworkflow
        PYRODIGAL ( ch_flye_contigs_dna, 'gff' )
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

        //
        // MODULE: run_dbCAN to find CAZymes and CGCs in metagenomics data
        //Will write it into subworkflow later

        ch_gunzip_gff_with_type = ch_gunzip_gff.map { meta, gff -> tuple(meta, gff, 'prodigal') }


        RUNDBCAN_EASYSUBSTRATE (
            ch_gunzip_faa,
            ch_gunzip_gff_with_type,
            ch_dbcan_db_final
        )

        ch_versions = ch_versions.mix(RUNDBCAN_EASYSUBSTRATE.out.versions)
        ch_dbcan_results = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results
        ch_dbcan_results_rna = RUNDBCAN_EASYSUBSTRATE.out.dbcan_results
            .map { meta, dbcan_results ->
                def new_meta = meta.clone()
                if (new_meta.id.endsWith('_dna')) {
                    new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
                }
                tuple(new_meta, dbcan_results)
            }
        //ch_flye_input_final_dna.view()

        ch_bwameme_input_dna = ch_flye_input_final_dna
        //ch_bwameme_input_dna.view()
        ch_bwameme_input_rna = ch_extracted_from_kraken2_reads_rna
        //ch_bwameme_input_rna.view()

        //
        BWA_INDEX(
            ch_flye_contigs_dna
        )

        ch_index_dna = BWA_INDEX.out.index

        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

        //ch_bwameme_input_dna.view()
        //
        // Read mapping with BWA-MEME
        ch_bwa_mem_input = ch_bwameme_input_dna
            .join(ch_index_dna)
            .join(ch_flye_contigs_dna)
            .map { meta, read1, index_dir, fasta ->
                tuple(meta, read1, index_dir, fasta)
            }
        // Long reads 只有单个文件,直接使用
            ch_minimap2_input = ch_bwameme_input_dna
                .join(ch_flye_contigs_dna, by: 0)
                .map { meta, reads, fasta ->
                    tuple(meta, reads, fasta)
                }

            ch_fastq = ch_minimap2_input.map { meta, reads, fasta -> tuple(meta, reads) }
            ch_genome_fna = ch_minimap2_input.map { meta, reads, fasta -> tuple(meta, fasta) }

            MINIMAP2_ALIGN (
                ch_fastq,
                ch_genome_fna,
                true,
                "bai",
                false,
                true,
            )

        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // DNA 的 BAM/BAI
        ch_bam_bai_dna = MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index, by: 0)
            .map { meta, bam, bai ->
                tuple(meta, bam, bai)
            }

        // 为 RNA 准备 BWA 输入 - 使用 DNA 的索引（因为它们使用相同的 contigs）
        // 基于样本名的 join key（去掉 _rna/_dna）
        ch_bwameme_input_rna_keyed = ch_bwameme_input_rna.map { meta, r1, r2 ->
            def key = meta.id.replaceFirst(/_rna$/, '')
            tuple(key, meta, r1, r2)
        }

        ch_index_dna_keyed = ch_index_dna.map { meta, idx ->
            def key = meta.id.replaceFirst(/_dna$/, '')
            tuple(key, meta, idx)
        }

        ch_contigs_dna_keyed = ch_flye_contigs_dna.map { meta, fasta ->
            def key = meta.id.replaceFirst(/_dna$/, '')
            tuple(key, meta, fasta)
        }

        // join RNA reads with DNA index/contigs using base key
        ch_bwa_mem_input_rna = ch_bwameme_input_rna_keyed
            .join(ch_index_dna_keyed, by: 0, remainder: false)
            .join(ch_contigs_dna_keyed, by: 0, remainder: false)
            .map { key, meta_rna, r1, r2, meta_idx, idx, meta_ctg, fasta ->
                // 让 BWA 使用 DNA meta（含 *_dna 前缀），以匹配已生成的索引前缀
                tuple(meta_idx, r1, r2, idx, fasta)
            }

        // BWA_MEM requires two paths (read1 and read2), so if read2 is null (single-end), use read1 as read2
        ch_fastq_rna = ch_bwa_mem_input_rna.map { meta, r1, r2, idx, fasta -> 
            def read2_final = r2 ?: r1  // Use read1 as read2 if read2 is null
            tuple(meta, r1, read2_final)
        }
        ch_index_rna_final = ch_bwa_mem_input_rna.map { meta, r1, r2, idx, fasta -> tuple(meta, idx) }
        ch_genome_fna_rna = ch_bwa_mem_input_rna.map { meta, r1, r2, idx, fasta -> tuple(meta, fasta) }

        BWA_INDEX_MEM_RNA (
            ch_index_rna_final,
            ch_genome_fna_rna,
            ch_fastq_rna
        )

        // 将输出 meta 从 *_dna 改回 *_rna，便于后续 join
        ch_bam_bai_rna = BWA_INDEX_MEM_RNA.out.ch_bam_bai.map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.id = new_meta.id.replaceFirst(/_dna$/, '_rna')
            tuple(new_meta, bam, bai)
        }
        ch_versions = ch_versions.mix(BWA_INDEX_MEM_RNA.out.versions)
        ch_gff_bam_bai_dna = ch_gunzip_gff
            .join(ch_bam_bai_dna, by: 0)
            .map { meta, gff, bam, bai ->
                tuple(meta, gff, bam, bai)
            }

        ch_gff_bam_bai_rna = ch_gunzip_gff_rna
            .join(ch_bam_bai_rna, by: 0)
            .map { meta, gff, bam, bai ->
                tuple(meta, gff, bam, bai)
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

        ch_abund_input_dna = ch_coverage_dna
            .join(ch_dbcan_results, by: 0)
            .map { meta, cgc_depth, dbcan_results -> tuple(meta, cgc_depth, dbcan_results) }

        ch_abund_input_rna = ch_coverage_rna
            .join(ch_dbcan_results_rna, by: 0)
            .map { meta, cgc_depth, dbcan_results -> tuple(meta, cgc_depth, dbcan_results) }


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
    // softwareVersionsToYAML 仅接受路径；将 ch_versions 统一成路径并过滤空
    def ch_versions_paths = ch_versions
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
