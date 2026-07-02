process RUNDBCAN_ASMFREE_FAM_ABUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.2--pyhdfd78af_0' }"
    
    input:
    tuple val(meta), path(diamond_r1), path(diamond_r2)  // For paired-end: two different results; for single-end: same result twice
    tuple val(meta2), path(reads)
    path db_dir     // Directory containing CAZyID_subfam_mapping.tsv
    val normalize   // Normalization method, e.g., 'FPKM'

    output:
    tuple val(meta), path("*_abund") , emit: abund_dir
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    def normalize_val = normalize ?: 'FPKM'
    def is_paired = !meta.single_end
    def paf_args = is_paired ? "-paf1 ../${diamond_r1} -paf2 ../${diamond_r2}" : "-paf1 ../${diamond_r1}"
    
    """
    mkdir -p ${prefix}
    ln -sfn \$(readlink -f ${db_dir}) ${prefix}/db
    cd ${prefix}

    dbcan_asmfree diamond_fam_abund \\
        ${paf_args} \\
        --raw_reads ../${reads} \\
        -n ${normalize_val} \\
        -o fam_abund.out \\
        ${args}

    rm -f db
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    """
    mkdir -p ${prefix}
    touch ${prefix}/fam_abund.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
