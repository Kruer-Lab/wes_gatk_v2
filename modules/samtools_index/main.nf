process SAMTOOLS_INDEX {
    label 'short'
    
    publishDir "${params.outDataDir}/${sample}", mode: 'copy'

    input:
        tuple val(family), val(sample), path(bam)

    output:
        tuple val(family), val(sample), path("${bam}.bai"), emit: index

    script:
    """
    samtools index ${bam} "${bam}.bai"
    """
}