process RUNDBCAN_UTILS_CAL_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    tuple val(meta2), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.depth.txt")               , emit: depth_txt
    path  "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    dbcan_utils cal_coverage \\
    -g ${gff} \\
    -i ${bam} \\
    -o ${prefix}.depth.txt \\
    --overlap_base_ratio 0.2 \\
    --mapping_quality 30 \\
    --identity 0.98 \\
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
    touch ${prefix}_depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
