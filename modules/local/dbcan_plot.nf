process RUNDBCAN_PLOT_BAR {
    tag "${meta_list.join('_')}"
    label 'process_medium'

    conda "bioconda::dbcan=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta_list), path(abund_dirs)

    output:
    tuple val(meta_list), path("*pdf"), emit: plot_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    abund_dirs && abund_dirs.size() > 0

    script:
    def args = task.ext.args ?: ''
    def meta_ids = meta_list instanceof List ? meta_list : [meta_list]
    def prefix = task.ext.prefix ?: "${meta_ids.join('_')}_pdf"
    def sample_names = meta_ids.join(',')
    def input_files_fam = abund_dirs.collect { "${it}/fam_abund.out" }.join(',')
    def input_files_subfam = abund_dirs.collect { "${it}/subfam_abund.out" }.join(',')
    def input_files_ec = abund_dirs.collect { "${it}/EC_abund.out" }.join(',')
    def input_files_fam_substrate = abund_dirs.collect { "${it}/fam_substrate_abund.out" }.join(',')

    """

    mkdir -p ${prefix}

    plots.py heatmap_plot \\
        --samples ${sample_names} \\
        -i ${input_files_fam_substrate} \\
        --show_abund \\
        --top 20 \\
        ${args}

    mv heatmap.pdf ${prefix}

    plots.py bar_plot \\
        --samples ${sample_names} \\
        --vertical_bar \\
        --top 20 \\
        -i ${input_files_fam} \\
        --pdf ${prefix}/fam.pdf \\
        ${args}

    plots.py bar_plot \\
        --samples ${sample_names} \\
        --vertical_bar \\
        --top 20 \\
        -i ${input_files_subfam} \\
        --pdf ${prefix}/subfam.pdf \\
        ${args}

    plots.py bar_plot \\
        --samples ${sample_names} \\
        --vertical_bar \\
        --top 20 \\
        -i ${input_files_ec} \\
        --pdf ${prefix}/ec.pdf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
