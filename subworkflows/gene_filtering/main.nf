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

include { COMBINE_AND_GENOTYPE_GVCFS_SINGLE } from '../../modules/gatk'
include { HOMOZYGOUS_RECESSIVE as HOMOZYGOUS_RECESSIVE_INCOMPLETE } from '../../modules/inheritanceModesIncompleteTrios'
include { X_LINKED_RECESSIVE as X_LINKED_RECESSIVE_INCOMPLETE } from '../../modules/inheritanceModesIncompleteTrios'
include { DOMINANT as DOMINANT_INCOMPLETE } from '../../modules/inheritanceModesIncompleteTrios'
include { VARIANTS_TO_TABLE as VARIANTS_TO_TABLE_INCOMPLETE } from '../../modules/inheritanceModesIncompleteTrios'



workflow GENE_FILTERING {
    take:
        rawGVCF_mother
        rawGVCF_father
        rawGVCF_child
        pedigrees

    main:

        // Get reference genome, index, and Annovar
        refGenome = Channel.fromPath(params.referenceGenome, checkIfExists: true).collect()
        refIndex = Channel.fromPath(["${params.referenceGenome}.fai",
                                        "${params.referenceGenome}".split("\\.")[0] + ".dict"], checkIfExists: true).collect()
        annovarRef = Channel.fromPath("${params.resourcesDir}/Annovar/humandb/", checkIfExists: true).collect()


        // Incomplete trio
        if(rawGVCF_mother == null && rawGVCF_father == null) {
            subsetChild = rawGVCF_child.map{fam, sample, gvcf ->
                                            tuple(fam, gvcf)}
            
            // Join gvcf with pedigree
            family_rawGVCF_pedigree = subsetChild.join(pedigrees)

            COMBINE_AND_GENOTYPE_GVCFS_SINGLE(family_rawGVCF_pedigree, refGenome, refIndex)

            HOMOZYGOUS_RECESSIVE_INCOMPLETE(COMBINE_AND_GENOTYPE_GVCFS_SINGLE.out, refGenome, refIndex, annovarRef)

            X_LINKED_RECESSIVE_INCOMPLETE(COMBINE_AND_GENOTYPE_GVCFS_SINGLE.out, refGenome, refIndex, annovarRef)

            DOMINANT_INCOMPLETE(COMBINE_AND_GENOTYPE_GVCFS_SINGLE.out, refGenome, refIndex, annovarRef)

            merged_variants_incomplete = HOMOZYGOUS_RECESSIVE_INCOMPLETE.out.homCaddMetaSVMRare
                .join(X_LINKED_RECESSIVE_INCOMPLETE.out.xLinkCaddMetaSVMRareVCF)
                .join(DOMINANT_INCOMPLETE.out.domCaddMetaSVMRareVCF)

            VARIANTS_TO_TABLE_INCOMPLETE(merged_variants_incomplete)

        }

        // Complete trio
        else {
            subsetMother = rawGVCF_mother.map{fam, sample, gvcf -> 
                                            tuple(fam, gvcf)}

            subsetfather = rawGVCF_father.map{fam, sample, gvcf -> 
                                            tuple(fam, gvcf)}

            subsetChild = rawGVCF_child.map{fam, sample, gvcf -> 
                                            tuple(fam, gvcf)}                                

            // Combine family raw GVCFs and pedigree
            family_rawGVCF_pedigree = subsetChild.join(subsetMother)
                .join(subsetfather)
                .join(pedigrees)

            // Combine GVCFs into one GVCF file for the trio and genotype
            COMBINE_AND_GENOTYPE_GVCFS(family_rawGVCF_pedigree, refGenome, refIndex)

            // Call denovo SNPs
            (denovo_snpIntersectVCF, denovo_snpSelectVCF) = DENOVO_SNPS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex, annovarRef)

            // Call denovo indels
            (denovo_indelIntersectVCF, denovo_indelSelectVCF) = DENOVO_INDELS(COMBINE_AND_GENOTYPE_GVCFS.out.genotype, refGenome, refIndex, annovarRef)

            // Merge DENOVO_SNPS and DENOVO_INDELS outputs
            denovo_variants = denovo_snpSelectVCF.join(denovo_indelSelectVCF)

            // Call homozygous recessive variants
            HOMOZYGOUS_RECESSIVE(denovo_variants, refGenome, refIndex, annovarRef)

            // Call compound heterozygous variants
            COMPOUND_HETEROZYGOUS(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex, annovarRef)

            // Call X-Linked recessive variants
            X_LINKED_RECESSIVE(HOMOZYGOUS_RECESSIVE.out.mergedVariantsVCF, refGenome, refIndex, annovarRef)

            // Call dominant variants
            DOMINANT(rawGVCF_child, refGenome, refIndex, annovarRef)

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
}