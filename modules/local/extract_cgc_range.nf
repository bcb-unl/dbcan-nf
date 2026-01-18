process EXTRACT_CGC_RANGES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dbcan_results)


    output:
    tuple val(meta), path("*cgc_ranges.tsv"), emit: cgc_ranges


    when:



    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_cgc_ranges.py ${dbcan_results}/cgc_standard_out.tsv ${dbcan_results}/substrate_prediction.tsv ${prefix}_cgc_ranges.tsv
    """

}
