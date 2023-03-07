process KINSHIP {
    label = 'short'

    publishDir "${params.outTrioDir}/Kinship/${family}", pattern: "*.kinship",mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(rawVCF)
        
    output:
        path kinshipFile

    script:
    """
    vcftools --vcf $rawVCF --relatedness2 --out ${family}.kinship
    """
}
