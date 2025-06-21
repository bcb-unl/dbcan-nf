process RUNDBCAN_EASYSUBSTRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(input_raw_data)
    tuple val(meta2), path(input_gff), val(gff_type)
    path  dbcan_db

    output:
    tuple val(meta), path("${prefix}_dbcan/")            , emit: dbcan_results
    path  "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir ${prefix}_dbcan \\
        --input_gff ${input_gff} \\
        --gff_type ${gff_type} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir -p ${prefix}_dbcan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
