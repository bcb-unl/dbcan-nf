process EXTRACT_CGC_RANGES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::dbcan=5.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.1.2--pyhdfd78af_0' :
        'biocontainers/dbcan:5.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cgc_standard_out)
    tuple val(meta2), path(substrate_prediction)


    output:
    tuple val(meta), path("cgc_ranges.tsv"), emit: cgc_ranges

    when:

    script:
    """
    extract_cgc_ranges.py ${cgc_standard_out} ${substrate_prediction} cgc_ranges.tsv
    """

}
