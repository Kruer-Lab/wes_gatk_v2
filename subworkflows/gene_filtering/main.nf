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
        (denovo_snpIntersectVCF, denovo_snpSelectVCF) = DENOVO_SNPS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex)

        // Call denovo indels
        (denovo_indelIntersectVCF, denovo_indelSelectVCF) = DENOVO_INDELS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex)

        // Merge DENOVO_SNPS and DENOVO_INDELS outputs
        denovo_variants = denovo_snpSelectVCF.join(denovo_indelSelectVCF)

        // Call homozygous recessive variants
        HOMOZYGOUS_RECESSIVE(denovo_variants, refGenome, refIndex)

        // Call compound heterozygous variants
        COMPOUND_HETEROZYGOUS(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex)

        // Call X-Linked recessive variants
        X_LINKED_RECESSIVE(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex)

        // Call dominant variants
        DOMINANT(rawGVCF_child, refGenome, refIndex)

        // Merge variant outputs
        merged_variants = denovo_snpIntersectVCF.join(denovo_indelIntersectVCF)
            .join(HOMOZYGOUS_RECESSIVE.out.homCaddMetaSVMRare)
            .join(COMPOUND_HETEROZYGOUS.out.comphetCaddMetaSVMRareVCF)
            .join(X_LINKED_RECESSIVE.out.xLinkCaddMetaSVMRareVCF)
            .join(DOMINANT.out.domCaddMetaSVMRareVCF)

        // Format outputs as .csv files
        VARIANTS_TO_TABLE(merged_variants)

        // Get kinship
        KINSHIP(COMBINE_AND_GENOTYPE_GVCFS.out.genotype)
}