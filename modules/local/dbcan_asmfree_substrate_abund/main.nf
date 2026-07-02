process RUNDBCAN_ASMFREE_SUBSTRATE_ABUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.2--pyhdfd78af_0' }"
    
    input:
    tuple val(meta), path(subfam_dir), path(fam_dir), path(ec_dir)  // All abundance directories combined
    path db_dir     // Directory containing fam-substrate-mapping.tsv

    output:
    tuple val(meta), path("*_abund") , emit: abund_dir
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    
    """
    mkdir -p ${prefix}
    ln -sfn \$(readlink -f ${db_dir}) ${prefix}/db
    cd ${prefix}

    dbcan_asmfree diamond_substrate_abund \\
        -i ../${subfam_dir}/subfam_abund.out \\
        -o substrate_abund.out \\
        ${args}

    rm -f db

    # Symlink other abundance outputs for RUNDBCAN_PLOT_BAR (expects all .out files in one directory)
    if [ -f ../${fam_dir}/fam_abund.out ]; then
        ln -sf \$(readlink -f ../${fam_dir}/fam_abund.out) fam_abund.out
    fi
    if [ -f ../${subfam_dir}/subfam_abund.out ]; then
        ln -sf \$(readlink -f ../${subfam_dir}/subfam_abund.out) subfam_abund.out
    fi
    if [ -f ../${ec_dir}/EC_abund.out ]; then
        ln -sf \$(readlink -f ../${ec_dir}/EC_abund.out) EC_abund.out
    fi

    if [ -f substrate_abund.out ]; then
        ln -sf substrate_abund.out fam_substrate_abund.out
    fi

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
    touch ${prefix}/substrate_abund.out
    touch ${prefix}/fam_substrate_abund.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
