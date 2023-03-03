process INDEX_REFERENCE {

    publishDir "${referenceDirectory}", mode: 'copy'

    input:
        val referenceDirectory
        path ref_index

    output:
        path("*.{amb,ann,bwt,pac,sa,fai,dict}"), emit: index

    script:
    def dictName = ref_index.name.split("\\.")[0] + ".dict"
    """
    bwa index "$ref_index"
    samtools faidx "$ref_index"
    samtools dict "$ref_index" -o "$dictName"
    """
}

process BWA_MEM {
    label 'alignment'
 
    input:
        path refGenome
        path index
        tuple val(family), val(sample), path(read1), path(read2)  

    output:
        tuple val(family), val(sample), path("${sample}.sorted.bam"), emit: sortedBam

    script:
    // Samtools memory input is memory per thread
    def samtoolsMemory = (task.memory.toBytes() / task.cpus).toString().split()[0]
    """
    bwa mem -M -t ${task.cpus} $refGenome $read1 $read2 | samtools view -bS | samtools sort -@ ${task.cpus} -m ${samtoolsMemory} > ${sample}.sorted.bam
    """
}