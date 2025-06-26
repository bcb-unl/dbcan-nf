include { EXTRACT_CGC_RANGES } from '../../modules/local/extract_cgc_range.nf'
include { SAMTOOLS_DEPTH     } from '../../modules/nf-core/samtools/depth/main.nf'
include { RUNDBCAN_PLOT_CGC  } from '../../modules/local/dbcan_plot_cgc.nf'

workflow CGC_DEPTH_PLOT {
    take:
        dbcan_folder // tuple(meta, path)
        ch_bam_bai   // tuple(meta, bam, bai)

    main:
        // 1. 从dbcan结果中提取每个CGC的坐标范围
        EXTRACT_CGC_RANGES(
            dbcan_folder
        )

        // 2. 读取范围文件，并为每个CGC创建独立的channel元素
        //    这是你的原始逻辑，是正确的
        ch_cgc_regions = EXTRACT_CGC_RANGES.out.cgc_ranges
            .flatMap { meta, cgc_ranges_file ->
                // 增加健壮性：如果文件为空则跳过，避免报错
                if (cgc_ranges_file.size() > 0) {
                    cgc_ranges_file.splitCsv(header:true, sep:'\t')
                        .findAll { row -> row['Contig ID'] && row['Gene Start'] && row['Gene Stop'] } // 确保行数据完整
                        .collect { row ->
                            def region = "${row['Contig ID']}:${row['Gene Start']}-${row['Gene Stop']}"
                            tuple(meta, row['CGCID'], region)
                        }
                } else {
                    [] // 如果文件为空，返回空列表
                }
            }

        // 3. 用 groupTuple+join+flatMap 实现一对多配对
        ch_cgc_depth_input = ch_cgc_regions
            .groupTuple() // 按meta分组，得到 (meta, [cgcid], [region])
            .join(ch_bam_bai, by: 0) // 按meta join bam/bai
            .flatMap { meta, cgcids, regions, bam, bai ->
                def results = []
                for (int i = 0; i < cgcids.size(); i++) {
                    results << tuple(meta, bam, bai, cgcids[i], regions[i])
                }
                return results
            }
        ch_cgc_depth_input.view()
        // 4. 为每个CGC计算深度
        SAMTOOLS_DEPTH(
            ch_cgc_depth_input
        )

        // 5. 将计算出的深度结果与dbcan文件夹连接，为绘图做准备
        ch_readscount = SAMTOOLS_DEPTH.out.tsv // 输出: tuple(meta, cgcid, tsv_path)

        ch_plot_input = ch_readscount
            .join(dbcan_folder, by: 0) // 按meta连接，这也是一个"一对多"连接
            .map { tsv_tuple, dbcan_tuple ->
                def meta          = tsv_tuple[0]
                def cgc_id        = tsv_tuple[1]
                def readscount    = tsv_tuple[2]
                def dbcan_results = dbcan_tuple[1]
                tuple(meta, dbcan_results, cgc_id, readscount)
            }

        // 6. 对每个CGC进行绘图
        RUNDBCAN_PLOT_CGC(
            ch_plot_input
        )

    emit:
        plot_dir = RUNDBCAN_PLOT_CGC.out.plot_dir
        tsv = SAMTOOLS_DEPTH.out.tsv
        versions = SAMTOOLS_DEPTH.out.versions.mix(RUNDBCAN_PLOT_CGC.out.versions)
}
