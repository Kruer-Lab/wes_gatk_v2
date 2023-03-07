process KINSHIP {
    label = 'short'

    publishDir "${params.outTrioDir}/Kinship/${family}", pattern: "*.kinship.relatedness2",mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(rawVCF)
        
    output:
        path "{family}.kinship.relatedness2"

    script:
    """
    vcftools --vcf $rawVCF --relatedness2 --out ${family}.kinship
    """
}
