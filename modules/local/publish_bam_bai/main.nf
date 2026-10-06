process PUBLISH_BAM_BAI {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai)

    script:
    """
    touch ${bam} ${bai}
    """


}