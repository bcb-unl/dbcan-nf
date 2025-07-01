#!/usr/bin/env nextflow

//
// Read mapping with BWA-MEME
//
include { BWA_INDEX  } from '../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM    } from '../../modules/nf-core/bwa/mem/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow BWA_INDEX_MEM {
    take:
        ch_genome_fna  // tuple(meta, genome_fna)
        ch_fastq       // tuple(meta, fastq)

    main:

        //1. build index
        BWA_INDEX(
            ch_genome_fna
        )
        // 2. map reads
        BWA_MEM(
            ch_fastq,
            BWA_INDEX.out.index,
            ch_genome_fna.map { meta, fna -> tuple(meta, fna) },
            'sort',   // sort_bam
        )
        // 3. index bam
        SAMTOOLS_INDEX(
            BWA_MEM.out.bam
        )
        // 4. prepare output channel with bam and its index (bai or csi)
        BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai ->
                    [ meta, bam, bai ]
            }
        .set { ch_bam_bai }

    emit:
        bam = BWA_MEM.out.bam
        bam_index = SAMTOOLS_INDEX.out.bai
        ch_bam_bai = ch_bam_bai
        versions = BWA_MEM.out.versions.mix(BWA_INDEX.out.versions)
}
