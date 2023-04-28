process PICARD_ADD_REPLACE_READ_GROUPS {
    label 'picard'

    input:
        tuple val(family), val(sample), path(sortedBam), path(samtoolsIndex)

    output:
        tuple val(family), val(sample), path("${sample}.grp.bam"), path("${sample}.grp.bam.bai")

    script:
    // Replace all read groups with a single new read group and assign all reads to this read group
    // Make index of fixed group file

    """
    java -jar \$PICARD AddOrReplaceReadGroups I=${sortedBam} O=${sample}.grp.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${sample}
    samtools index ${sample}.grp.bam ${sample}.grp.bam.bai
    """
}

process PICARD_MARK_DUPLICATES {
    label 'picard'
    
    input:
        tuple val(family), val(sample), path(grpBam), path(index)
    
    output:
        tuple val(family), val(sample), path("${sample}.metrics"), emit: metrics
        tuple val(family), val(sample), path("${sample}.dedupBam"), path("${sample}.dedupBam.bai"), emit: dedupBam

    script:
    """
    java -jar \$PICARD MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE="${sample}.metrics" REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT INPUT="${grpBam}" OUTPUT="${sample}.dedupBam"
    samtools index ${sample}.dedupBam ${sample}.dedupBam.bai
    """
}