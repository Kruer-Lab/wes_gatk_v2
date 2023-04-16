process FASTQC {
    label 'short'

    publishDir "${params.qcDir}", mode: 'copy'

    input:
        path files

    output:
        path "*"

    script:
    """
    fastqc $files --threads $task.cpus
    """
}

process MULTIQC {
    label 'short'

    publishDir "${params.qcDir}", mode: 'copy'

    input:
        path qcFiles

    output:
        path "*"

    script:
    """
    multiqc $qcFiles
    """
}