// CALL_VARIANTS
// Performs GATK best practices up until HaplotypeCaller

nextflow.enable.dsl=2

include { INDEX_REFERENCE; BWA_MEM } from '../../modules/alignment'
include { SAMTOOLS_INDEX } from '../../modules/samtools_index'
include { PICARD_ADD_REPLACE_READ_GROUPS; PICARD_MARK_DUPLICATES } from '../../modules/picard'
include { BQSR; HAPLOTYPE_CALLER } from '../../modules/gatk'

workflow CALL_VARIANTS {
    take:
        read_pairs

    main:

        // Check for reference genome
        refGenome = Channel.fromPath(params.referenceGenome, checkIfExists: true).collect()
        referenceDirectory = params.referenceGenome.substring(0, params.referenceGenome.lastIndexOf("/")) + "/"

        // Check for dbSNP and index
        dbSNP = Channel.fromPath(params.dbSNP, checkIfExists: true).collect()
        dbSNP_index = Channel.fromPath("${params.dbSNP}.idx", checkIfExists: true).collect()

        // Check for interval list
        intervalList = Channel.fromPath(params.intervalList, checkIfExists: true).collect()

        // Index reference genome
        if(params.index_reference) {
            refIndex = INDEX_REFERENCE(referenceDirectory, refGenome).collect()
        }

        // Skip indexing
        else {
            refIndex = Channel.fromPath(["${params.referenceGenome}.{amb,ann,bwt,pac,sa,fai}",
                                        "${params.referenceGenome}".split("\\.")[0] + ".dict"], checkIfExists: true).collect()
        }

        // Align to reference genome
        BWA_MEM(refGenome, refIndex, read_pairs)

        // Make index of BAM file
        SAMTOOLS_INDEX(BWA_MEM.out.sortedBam)

        // Replace all read groups with a single new read group and assign all reads to this read group
        PICARD_ADD_REPLACE_READ_GROUPS(BWA_MEM.out.sortedBam, SAMTOOLS_INDEX.out.index)

        // Mark duplicates
        PICARD_MARK_DUPLICATES(PICARD_ADD_REPLACE_READ_GROUPS.out.grpBam, PICARD_ADD_REPLACE_READ_GROUPS.out.index)

        // Base Quality Score Recalibration
        BQSR(PICARD_MARK_DUPLICATES.out.dedupBam, PICARD_MARK_DUPLICATES.out.index, refGenome, refIndex, dbSNP, dbSNP_index)

        // Germline variant calling by HaplotypeCaller
        HAPLOTYPE_CALLER(BQSR.out.bqrecalBam, refGenome, refIndex, intervalList)

    emit:
        rawGVCF = HAPLOTYPE_CALLER.out.rawGVCF
}