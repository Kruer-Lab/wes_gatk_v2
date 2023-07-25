process BQSR {
    label 'medium'

    publishDir "${params.outDataDir}/${sample}", mode: 'copy'

    input:
        tuple val(family), val(sample), path(dedupBam), path(dedupIndex)
        path refGenome
        path refIndex
        path dbSNP
        path dbSNP_index

    output:
        tuple val(family), val(sample), path("${sample}.bqrecal.table"), emit: bqrecalTable
        tuple val(family), val(sample), path("${sample}.bqrecal.bam"), emit: bqrecalBam
        path("*.bai"), emit: index

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" BaseRecalibrator -R "${refGenome}" -I "${dedupBam}" --known-sites "${dbSNP}" -O "${sample}.bqrecal.table"
    gatk --java-options "-Xmx${gatkMemory}g" ApplyBQSR -R "${refGenome}" -I "${dedupBam}" --bqsr-recal-file "${sample}.bqrecal.table" -O "${sample}.bqrecal.bam"
    """
}

process HAPLOTYPE_CALLER {
    label 'medium'
    
    publishDir "${params.outDataDir}/${sample}", mode: 'copy'

    input:
        tuple val(family), val(sample), path(bqrecalBam)
        path refGenome
        path refIndex
        path intervalList

    output:
        tuple val(family), val(sample), path("${sample}.raw.g.vcf"), emit: rawGVCF
        tuple val(family), val(sample), path("${sample}.raw.g.vcf.idx"), emit: rawGVFIdx

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" HaplotypeCaller -R "${refGenome}" -I "${bqrecalBam}" -ERC BP_RESOLUTION --do-not-run-physical-phasing true -L "${intervalList}" -O "${sample}.raw.g.vcf"
    """
}

process COMBINE_AND_GENOTYPE_GVCFS {
    label = 'medium'

    publishDir "${params.outTrioDir}/${family}", mode: 'copy'

    input:
        tuple val(family), path(rawGVCF_child), path(rawGVCF_mother), path(rawGVCF_father), path(pedigree)
        path refGenome
        path refIndex

    output:
        tuple val(family), path(pedigree), path("${family}.trio.raw.vcf"), emit: genotype

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" CombineGVCFs -R $refGenome --variant $rawGVCF_child --variant $rawGVCF_mother --variant $rawGVCF_father -O "${family}.comb.trio.raw.vcf"
    gatk --java-options "-Xmx${gatkMemory}g" GenotypeGVCFs -R $refGenome --variant "${family}.comb.trio.raw.vcf" -O "${family}.trio.raw.vcf"
    """
}

process COMBINE_AND_GENOTYPE_GVCFS_SINGLE {
    label = 'medium'

    publishDir "${params.outTrioDir}/${family}", mode: 'copy'

    input:
        tuple val(family), path(rawGVCF_child), path(pedigree)
        path refGenome
        path refIndex

    output:
        tuple val(family), path(pedigree), path("${family}.trio.raw.vcf"), emit: genotype

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" CombineGVCFs -R $refGenome --variant $rawGVCF_child -O "${family}.comb.trio.raw.vcf"
    gatk --java-options "-Xmx${gatkMemory}g" GenotypeGVCFs -R $refGenome --variant "${family}.comb.trio.raw.vcf" -O "${family}.trio.raw.vcf"
    """

}
