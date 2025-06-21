include { EXTRACT_CGC_RANGES } from '../../modules/local/extract_cgc_range.nf'
include { SAMTOOLS_DEPTH     } from '../../modules/nf-core/samtools/depth/main.nf'

workflow CGC_DEPTH {
    take:
        dbcan_folder

        ch_bam_bai  // tuple(meta, bam, bai)

    main:
        // 1. extract cgc_ranges.tsv
        ch_cgc_standard_out = dbcan_folder.map { meta, folder -> tuple(meta, file("${folder}/cgc_standard_out.tsv")) }
        ch_substrate_prediction = dbcan_folder.map { meta, folder -> tuple(meta, file("${folder}/substrate_prediction.tsv")) }

        EXTRACT_CGC_RANGES(
            ch_cgc_standard_out,
            ch_substrate_prediction
        )

        // 2. read cgc_ranges.tsv and generate region channel
        ch_cgc_regions = EXTRACT_CGC_RANGES.out.cgc_ranges
            .flatMap { meta, cgc_ranges_file ->
                cgc_ranges_file.splitCsv(header:true, sep:'\t')
                    .collect { row ->
                        def region = "${row['Contig ID']}:${row['Gene Start']}-${row['Gene Stop']}"
                        tuple(meta, row['CGCID'], region)
                    }
            }
         //ch_cgc_regions.view() // Uncomment to debug the cgc_regions channel
         //ch_bam_bai.view()

        // 3. combine cgc_regions with bam_channel to prepare input for samtools depth
        //  bam_channel: tuple(meta, bam, bai)
        ch_cgc_depth_input = ch_cgc_regions
            .combine(ch_bam_bai)
            .map { meta, cgcid, region, meta2, bam, bai ->
                // meta: sample meta
                // cgcid: CGC id
                // region: region string
                // meta2: sample meta (again, from bam_bai)
                // bam, bai: bam/bai files
                tuple(meta, cgcid, region, bam, bai)
            }

        //ch_cgc_depth_input.view()

        SAMTOOLS_DEPTH(
            ch_cgc_depth_input.map{ meta, cgcid, region, bam, bai -> tuple(meta, bam, bai) },
            ch_cgc_depth_input.map{ meta, cgcid, region, bam, bai -> tuple([cgcid: cgcid], region) }
        )

    emit:
        tsv = SAMTOOLS_DEPTH.out.tsv
        versions = SAMTOOLS_DEPTH.out.versions
}
