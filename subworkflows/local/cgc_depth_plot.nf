include { EXTRACT_CGC_RANGES } from '../../modules/local/extract_cgc_range.nf'
include { SAMTOOLS_DEPTH     } from '../../modules/nf-core/samtools/depth/main.nf'
include { RUNDBCAN_PLOT_CGC  } from '../../modules/local/dbcan_plot_cgc.nf'

workflow CGC_DEPTH_PLOT {
    take:
        dbcan_folder // tuple(meta, path)
        ch_bam_bai   // tuple(meta, bam, bai)
        dbcan_db     // tuple(meta, path)

    main:
        // 1. extract CGC ranges from dbcan results
        EXTRACT_CGC_RANGES(
            dbcan_folder
        )

        // 2. read CGC ranges and create (meta, cgcid, region) tuples
        //
        ch_cgc_regions = EXTRACT_CGC_RANGES.out.cgc_ranges
            .flatMap { meta, cgc_ranges_file ->
                //
                if (cgc_ranges_file.size() > 0) {
                    cgc_ranges_file.splitCsv(header:true, sep:'\t')
                        .findAll { row -> row['Contig ID'] && row['Gene Start'] && row['Gene Stop'] } // 确保行数据完整
                        .collect { row ->
                            def region = "${row['Contig ID']}:${row['Gene Start']}-${row['Gene Stop']}"
                            tuple(meta, row['CGCID'], region)
                        }
                } else {
                    [] //
                }
            }
        //ch_cgc_regions.view()
        //ch_bam_bai.view()
        // combine to (meta, cgcid, region, bam, bai)
        ch_cgc_depth_input = ch_cgc_regions
            .combine(ch_bam_bai, by:0) //

        //ch_cgc_depth_input.view()

        // 4. calculate depth for each CGC region in each sample
        SAMTOOLS_DEPTH(
            ch_cgc_depth_input
        )

        // 5. link depth results with dbcan results for plotting
        ch_readscount = SAMTOOLS_DEPTH.out.tsv // : tuple(meta, tsv_path)

        ch_plot_input = ch_readscount
           .combine(dbcan_folder, by: 0) // connect with dbcan results, tuple(meta, cgcid, tsv_path) + tuple(meta, dbcan_path)
            .map {meta, readscount, dbcan_results  ->
                def cgc_id        = readscount.getBaseName().replaceAll('-', /\|/) // get cgcid from tsv filename
                tuple(meta, dbcan_results, cgc_id, readscount)
            }
        //ch_plot_input.view()

        // 6. plot
        RUNDBCAN_PLOT_CGC(
            ch_plot_input,
            dbcan_db
        )

    emit:
        cgc_abund_pdf = RUNDBCAN_PLOT_CGC.out.cgc_abund_pdf
        tsv = SAMTOOLS_DEPTH.out.tsv
        versions = SAMTOOLS_DEPTH.out.versions
}
