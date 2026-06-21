process SEQTK_SAMPLE_PAIRED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(read1), path(read2), val(sample_spec)

    output:
    tuple val(meta), path(read1), path(read2), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(args ==~ /.*\ -s\ ?[0-9]+.*/)) {
        args += " -s100"
    }
    if ( sample_spec == null ) {
        error "SEQTK_SAMPLE_PAIRED must have a sample_spec value included"
    }
    """
    seqtk \\
        sample \\
        $args \\
        ${read1} \\
        ${sample_spec} \\
        | seqtk seq -n \\
        | sort -u > read_names_r1.txt

    sed 's|/1\$|/2|' read_names_r1.txt > read_names_r2.txt

    seqtk \\
        subseq \\
        ${read1} \\
        read_names_r1.txt \\
        | gzip --no-name > ${prefix}_\$(basename ${read1})

    seqtk \\
        subseq \\
        ${read2} \\
        read_names_r2.txt \\
        | gzip --no-name > ${prefix}_\$(basename ${read2})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}_\$(basename ${read1})
    echo "" | gzip > ${prefix}_\$(basename ${read2})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
