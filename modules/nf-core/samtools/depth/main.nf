process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(cgcid), val(intervals)

    output:
    tuple val(meta), val(cgcid), path("*_samtools_depth/*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_id = meta.sample ?: meta.id
    def outdir = "${sample_id}_samtools_depth"
    def prefix = cgcid.replaceAll(/\|/, '-')
    def positions = intervals ? "-r ${intervals}" : ""
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    mkdir -p ${outdir}
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        $positions \\
        -o ${outdir}/${prefix}.tsv \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
