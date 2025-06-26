process RUNDBCAN_PLOT_CGC {
    tag "${meta.id} cgc_plot"
    label 'process_medium'

    conda "bioconda::dbcan=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dbcan_results), val(cgc_id), path(readscount)

    output:
    tuple val(meta), path("*_cgc_plot"), emit: plot_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_cgc_plot"

    """
    mkdir -p ${prefix}_cgc_plot

    plots.py CGC_synteny_coverage_plot \\
        -i ${dbcan_results} \\
        --cgc_id '${cgc_id}' \\
        --readscount ${readscount} \\
        ${args}
    mv *syntenic-cov.pdf ${prefix}_cgc_plot/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
