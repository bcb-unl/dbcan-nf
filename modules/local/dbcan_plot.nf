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
    tuple val(meta_list), path("*_barpdf"), emit: bar_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def meta_ids = meta_list //
    def prefix = task.ext.prefix ?: "${meta_ids.join('_')}_barpdf"
    def sample_names = meta_ids.join(',')
    def input_files = abund_dirs.collect { "${it}/fam_abund.out" }.join(',')
    """
    mkdir -p ${prefix}
    plots.py bar_plot \\
        --samples ${sample_names} \\
        --vertical_bar \\
        --top 20 \\
        -i ${input_files} \\
        --pdf ${prefix}/fam.pdf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
