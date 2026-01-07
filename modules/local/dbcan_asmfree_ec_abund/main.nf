process RUNDBCAN_ASMFREE_EC_ABUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.2--pyhdfd78af_0' }"
    
    input:
    tuple val(meta), path(input_dir)  // Input directory from subfam_abund
    path db_dir     // Directory containing subfam_EC_mapping.tsv

    output:
    tuple val(meta), path("*_abund") , emit: abund_dir
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    
    // Copy db directory into per-sample folder so dbcan_asmfree can find subfam_EC_mapping.tsv
    def cmd = """
    mkdir -p ${prefix}
    cp -r ${db_dir} ${prefix}/db || true
    cd ${prefix}

    dbcan_asmfree diamond_EC_abund"""
    cmd += " -i ../${input_dir}"
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

