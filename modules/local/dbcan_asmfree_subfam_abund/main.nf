process RUNDBCAN_ASMFREE_SUBFAM_ABUND {
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
    
    // Determine if paired-end based on meta.single_end
    def is_paired = !meta.single_end
    
    // Copy db directory into per-sample folder so dbcan_asmfree can find CAZyID_subfam_mapping.tsv
    def cmd = """
    mkdir -p ${prefix}
    cp -r ${db_dir} ${prefix}/db || true
    cd ${prefix}
    dbcan_asmfree diamond_subfam_abund"""
    if (is_paired) {
        // Paired-end: use two different blastx results (paf1 and paf2)
        // We are now inside ${prefix}, so use paths relative to parent work dir
        cmd += " -paf1 ../${diamond_r1} -paf2 ../${diamond_r2}"
    } else {
        // Single-end: only paf1 (diamond_r1 and diamond_r2 are the same)
        cmd += " -paf1 ../${diamond_r1}"
    }
    cmd += " --raw_reads ../${reads}"
    cmd += " -n ${normalize_val}"
    cmd += " -o ${prefix}"
    cmd += " ${args}"
    
    """
    ${cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}

