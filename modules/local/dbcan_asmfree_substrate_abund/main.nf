process RUNDBCAN_ASMFREE_SUBSTRATE_ABUND {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.2--pyhdfd78af_0' }"
    
    input:
    tuple val(meta), path(subfam_dir), path(fam_dir), path(ec_dir)  // All abundance directories combined
    // Note: input_dir is the same as subfam_dir (substrate_abund uses subfam_abund output as input)
    path db_dir     // Directory containing fam-substrate-mapping.tsv

    output:
    tuple val(meta), path("*_abund") , emit: abund_dir
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_abund"
    
    // Copy db directory into per-sample folder so dbcan_asmfree can find fam-substrate-mapping.tsv
    def cmd = """
    mkdir -p ${prefix}
    cp -r ${db_dir} ${prefix}/db || true
    cd ${prefix}
    dbcan_asmfree diamond_substrate_abund"""
    cmd += " -i ../${subfam_dir}"
    cmd += " -o ${prefix}"
    cmd += " ${args}"
    
    // Create symlinks to other abundance files in the output directory
    // This makes the output directory compatible with RUNDBCAN_PLOT_BAR
    // which expects all .out files in one directory
    // Note: We are inside ${prefix} directory, so other directories are in parent directory
    cmd += """
    
    # Create symlinks to other abundance .out files in the output directory
    # Find and link fam_abund.out from fam_dir
    if [ -d "../${fam_dir}" ]; then
        FAM_OUT=\$(find ../${fam_dir} -name "*.out" -type f | head -n1)
        if [ -n "\$FAM_OUT" ]; then
            ln -sf \$(readlink -f "\$FAM_OUT") fam_abund.out
        fi
    fi
    
    # Find and link subfam_abund.out from subfam_dir
    if [ -d "../${subfam_dir}" ]; then
        SUBFAM_OUT=\$(find ../${subfam_dir} -name "*.out" -type f | head -n1)
        if [ -n "\$SUBFAM_OUT" ]; then
            ln -sf \$(readlink -f "\$SUBFAM_OUT") subfam_abund.out
        fi
    fi
    
    # Find and link EC_abund.out from ec_dir
    if [ -d "../${ec_dir}" ]; then
        EC_OUT=\$(find ../${ec_dir} -name "*.out" -type f | head -n1)
        if [ -n "\$EC_OUT" ]; then
            ln -sf \$(readlink -f \"\$EC_OUT\") EC_abund.out
        fi
    fi
    
    # Rename substrate_abund.out to fam_substrate_abund.out (dbcan_plot expects this name)
    if [ -f substrate_abund.out ]; then
        ln -sf substrate_abund.out fam_substrate_abund.out
    elif [ -f \$(find . -name "*.out" -type f | head -n1) ]; then
        SUBSTRATE_OUT=\$(find . -name "*.out" -type f | head -n1)
        ln -sf \$(basename \"\$SUBSTRATE_OUT\") fam_substrate_abund.out
    fi
    """
    
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

