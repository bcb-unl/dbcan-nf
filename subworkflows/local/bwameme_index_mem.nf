#!/usr/bin/env nextflow

//
// Read mapping with BWA-MEME
//
include { BWAMEME_INDEX  } from '../../modules/nf-core/bwameme/index/main.nf'
include { BWAMEME_MEM    } from '../../modules/nf-core/bwameme/mem/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow BWAMEME_INDEX_MEM {
    take:
        ch_genome_fna  // tuple(meta, genome_fna)
        ch_fastq       // tuple(meta, fastq)

    main:

        //1. build index
        BWAMEME_INDEX(
            ch_genome_fna
        )
        // 2. map reads
        BWAMEME_MEM(
            ch_fastq,
            BWAMEME_INDEX.out.index,
            ch_genome_fna.map { meta, fna -> tuple(meta, fna) },
            true,   // sort_bam
            3072,    // mbuffer (MB)
            2        // samtools_threads
        )
        // 3. index bam
        SAMTOOLS_INDEX(
            BWAMEME_MEM.out.bam
        )
        // 4. prepare output channel with bam and its index (bai or csi)
        BWAMEME_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai ->
                    [ meta, bam, bai ]
            }
        .set { ch_bam_bai }

    emit:
        bam = BWAMEME_MEM.out.bam
        bam_index = SAMTOOLS_INDEX.out.bai
        ch_bam_bai = ch_bam_bai
        versions = BWAMEME_MEM.out.versions.mix(BWAMEME_INDEX.out.versions)
}
