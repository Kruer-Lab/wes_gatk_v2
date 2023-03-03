// GENE_FILTERING
// Filters genes by inheritance mode

nextflow.enable.dsl=2

include { COMBINE_AND_GENOTYPE_GVCFS } from '../../modules/gatk'
include { DENOVO_SNPS } from '../../modules/inheritanceModes'
include { DENOVO_INDELS } from '../../modules/inheritanceModes'
include { HOMOZYGOUS_RECESSIVE } from '../../modules/inheritanceModes'
include { COMPOUND_HETEROZYGOUS } from '../../modules/inheritanceModes'
include { X_LINKED_RECESSIVE } from '../../modules/inheritanceModes'
include { DOMINANT } from '../../modules/inheritanceModes'
include { VARIANTS_TO_TABLE } from '../../modules/inheritanceModes'
include { KINSHIP } from '../../modules/kinship'

workflow GENE_FILTERING {
    take:
        rawGVCF_mother
        rawGVCF_father
        rawGVCF_child
        pedigrees

    main:

        // Get reference genome and index
        refGenome = Channel.fromPath(params.referenceGenome, checkIfExists: true).collect()
        refIndex = Channel.fromPath(["${params.referenceGenome}.fai",
                                        "${params.referenceGenome}".split("\\.")[0] + ".dict"], checkIfExists: true).collect()

        // Combine family raw GVCFs and pedigree
        family_rawGVCF_pedigree = rawGVCF_child.join(rawGVCF_mother)
            .join(rawGVCF_father)
            .join(pedigrees)

        // Combine GVCFs into one GVCF file for the trio and genotype
        COMBINE_AND_GENOTYPE_GVCFS(family_rawGVCF_pedigree, refGenome, refIndex)

        // Call denovo SNPs
        DENOVO_SNPS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex)

        // Call denovo indels
        DENOVO_INDELS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex)

        // Call homozygous recessive variants
        HOMOZYGOUS_RECESSIVE(DENOVO_SNPS.out.snpSelectVCF, DENOVO_INDELS.out.indelSelectVCF, refGenome, refIndex)

        // Call compound heterozygous variants
        COMPOUND_HETEROZYGOUS(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex)

        // Call X-Linked recessive variants
        X_LINKED_RECESSIVE(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex)

        // Call dominant variants
        DOMINANT(rawGVCF_child, refGenome, refIndex)

        // Format outputs as .csv files
        VARIANTS_TO_TABLE(
                        DENOVO_SNPS.out.denovoIntersectVCF,
                        DENOVO_INDELS.out.denovoIntersectVCF,
                        HOMOZYGOUS_RECESSIVE.out.homCaddMetaSVMRare,
                        COMPOUND_HETEROZYGOUS.out.comphetCaddMetaSVMRareVCF,
                        X_LINKED_RECESSIVE.out.xLinkCaddMetaSVMRareVCF,
                        DOMINANT.out.domCaddMetaSVMRareVCF)

        // Get kinship
        KINSHIP(COMBINE_AND_GENOTYPE_GVCFS.out.genotype)
}