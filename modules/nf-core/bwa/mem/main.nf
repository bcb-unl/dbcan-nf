process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf7890f8d4e38a7586581cb7fa13401b7af1582f21d94eef969df4cea852b6da/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4' }"

    input:
    tuple val(meta) , path(read1), path(read2)
    tuple val(meta), path(index)
    tuple val(meta2), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam")  , emit: bam,    optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.csi")  , emit: csi,    optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    
    // In coassembly mode, create symbolic link from sample index dir to coassembly index dir
    def index_prefix = meta.coassembly_index_id ? meta.coassembly_index_id : meta.id
    def index_setup = meta.coassembly_index_id ? """
    # Create symbolic link from sample index dir to coassembly index dir
    if [ ! -d "bwa_${meta.id}" ] && [ -d "bwa_${index_prefix}" ]; then
        ln -s bwa_${index_prefix} bwa_${meta.id}
    fi
    """ : ""
    
    """
    ${index_setup}
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        bwa_${meta.id}/bwa_${index_prefix} \\
        $read1 $read2 | samtools $samtools_command $args2 ${reference} --threads $task.cpus -o ${prefix}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.csi
    touch ${prefix}.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
