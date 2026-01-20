process DOWNLOAD_CAZY_DMND {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.16--h13889ed_0' :
        'quay.io/biocontainers/diamond:2.1.16--h13889ed_0' }"

    output:
    path "CAZy.dmnd", emit: cazy_dmnd
    path "db/CAZyID_subfam_mapping.tsv", emit: cazyid_subfam_mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cazy_url = task.ext.cazy_url ?: 'https://dbcan.s3.us-west-2.amazonaws.com/dbcan_asmfree/CAZyDB.07242025.fa'
    def cazyid_subfam_url = task.ext.cazyid_subfam_url ?: 'https://dbcan.s3.us-west-2.amazonaws.com/dbcan_asmfree/CAZyID_subfam_mapping.tsv'
    def fam_sub_mapping_url = task.ext.fam_sub_mapping_url ?: 'https://dbcan.s3.us-west-2.amazonaws.com/dbcan_asmfree/fam-substrate-mapping.tsv'
    def subfam_ec_mapping_url = task.ext.ec_mapping_url ?: 'https://dbcan.s3.us-west-2.amazonaws.com/dbcan_asmfree/subfam_EC_mapping.tsv'
    """
    # Download CAZyDB FASTA file
    wget -c ${cazy_url} -O CAZyDB.fa
    
    # Download CAZyID_subfam_mapping.tsv and create db directory
    mkdir -p db
    wget -c ${cazyid_subfam_url} -O db/CAZyID_subfam_mapping.tsv
    wget -c ${fam_sub_mapping_url} -O db/fam-substrate-mapping.tsv
    wget -c ${subfam_ec_mapping_url} -O db/subfam_EC_mapping.tsv
    
    # Build DIAMOND database from CAZyDB
    diamond makedb \\
        --in CAZyDB.fa \\
        --db CAZy \\
        --threads ${task.cpus} \\
        ${args}

    # Clean up FASTA file to save space
    rm -f CAZyDB.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch CAZy.dmnd
    mkdir -p db
    touch db/CAZyID_subfam_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}

