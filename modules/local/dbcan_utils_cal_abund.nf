process RUNDBCAN_UTILS_CAL_ABUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cgc_depth)
    tuple val(meta2), path(dbcan_results)

    output:
    tuple val(meta), path("*_abund")                  , emit: abund_dir
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_abund"
    def type = task.ext.type ?: 'TPM'
    """
    mkdir -p ${prefix}
    cd ${prefix}

    dbcan_utils fam_abund \\
        -bt ../${cgc_depth} \\
        -i ../${dbcan_results} \\
        -a ${type} \\
        ${args}

    dbcan_utils fam_substrate_abund \\
        -bt ../${cgc_depth} \\
        -i ../${dbcan_results} \\
        -a ${type} \\
        ${args}

    dbcan_utils CGC_abund  \\
        -bt ../${cgc_depth} \\
        -i ../${dbcan_results} \\
        -a ${type} \\
        ${args}

    dbcan_utils CGC_substrate_abund \\
        -bt ../${cgc_depth} \\
        -i ../${dbcan_results} \\
        -a ${type} \\
        ${args}
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
