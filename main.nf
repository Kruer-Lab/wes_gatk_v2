nextflow.enable.dsl=2

include { CALL_VARIANTS as CALL_VAR_M } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_F } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_C1 } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_C2 } from './subworkflows/call_variants'
include { GENE_FILTERING } from './subworkflows/gene_filtering'

workflow {

    // Regular expression patterns for getting sample and family IDs from file names
    def sampleRegexPattern = ~/F\d{3,}.*-\d{3}/
    def familyRegexPattern = ~/F\d{3,}/

    // Ensure sample sheet has complete trios
    def sheet = file(params.samplesheet, checkIfExists: true).readLines()*.split('\t') // Import sample sheet as .tsv
    sheet.remove(0) // Remove header

    def sheetFamilies = []
    for(line : sheet) {
        if(line.size() < 3) {
            error("ERROR: Sample sheet has incomplete trios.")
        }
        sheetFamilies.add((line[0] =~ familyRegexPattern).findAll()[0]) // Get family IDs
    }

    // Get family IDs for files in pedigree directory
    def pedDir = file(params.pedigreeDir).list()
    def pedDirFams = []
    for(item : pedDir) {
        if(item.endsWith(".ped")) {
            pedDirFams.add((item =~ familyRegexPattern).findAll()[0])
        }
    }

    // Error if no pedigrees are found
    if(pedDirFams.size() == 0) {
        error("ERROR: No pedigree families found.")
    }
    
    // Check pedigree directory for duplicates
    if(pedDirFams.toSet().size() != pedDirFams.size()) {
        error("ERROR: Duplicate sample pedigrees in pedigree directory.")
    }

    // Check for sample sheet families in pedigree families
    for(fam : sheetFamilies) {
        if(!pedDirFams.contains(fam)) {
            error("ERROR: Pedigree not found for " + fam)
        }
    }

    // Run variant calling (everything up to HaplotypeCaller)
    if(params.run_variant_calling == true) {

        // Import sample sheet and get reads
        // Get mother samples
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Mother =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row.Mother}/*_R1_*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row.Mother}/*_R2_*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of reads in ${row.Mother}"
                }

            return tuple(familyID, sampleID, r1, r2)
        }
            .set{read_pairs_mother}

        // Get father samples
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Father =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row.Father}/*_R1_*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row.Father}/*_R2_*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of reads in ${row.Father}"
                }

                return tuple(familyID, sampleID, r1, r2)
            }
            .set{read_pairs_father}

        // Get affected child samples
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Child_Affected =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row.Child_Affected}/*_R1_*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row.Child_Affected}/*_R2_*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of reads in ${row.Child_Affected}"
                }

                return tuple(familyID, sampleID, r1, r2)
            }
            .set{read_pairs_child_affected}

        // Get second child samples for quads (if any)
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                if(row.Child_Other) {
                    // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                    def sampleID = (row.Child_Other =~ sampleRegexPattern).findAll()[0]

                    // Get family ID ("F" number)
                    def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                    def r1 = file("$params.inDataDir/${row.Mother}/*_R1_*.fastq.gz", checkIfExists: true)
                    def r2 = file("$params.inDataDir/${row.Mother}/*_R2_*.fastq.gz", checkIfExists: true)

                    // Ensure only one file is used for each read
                    if(r1.size() != 1 || r2.size != 1) {
                        error "Incorrect number of reads in ${row.Mother}"
                    }

                    return tuple(familyID, sampleID, r1, r2)
                }
            }
            .set{read_pairs_child_other}

        CALL_VAR_M(read_pairs_mother)
        CALL_VAR_F(read_pairs_father)
        CALL_VAR_C1(read_pairs_child_affected)
        CALL_VAR_C2(read_pairs_child_other)
    }

    // Skipping variant calling
    else if(params.run_variant_calling == false) {

        // Get mother gvcf
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Mother =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                return tuple(familyID, gvcf)
            }
            .set{mother_calls}

        // Get father gvcf
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Father =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                return tuple(familyID, gvcf)
            }
            .set{father_calls}

        // Get affected child gvcf
        Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .ifEmpty{exit 1, "No sample sheet found at $params.samplesheet"}
            .splitCsv(header: true, sep: "\t", strip: true)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row.Child_Affected =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                return tuple(familyID, gvcf)
            }
            .set{proband_calls}

    }

    // Variant filtration
    if(params.run_variant_filtering == true) {

        // Get pedigrees
        Channel
        .fromPath("$params.pedigreeDir/*.ped", checkIfExists: true)
        .ifEmpty{exit 1, "No pedigrees found in $params.pedigreeDir"}
        .map{file ->
            // Get family ID (letter, 3 or more digits, -, 3 digits, -, letter (capture group gets just the 'F' number))
            // Alternatively will get family names with just the F### number
            def familyID = (file =~ familyRegexPattern).findAll()[0]

            return tuple(familyID, file)
        }
        .set{pedigrees}

        if(params.run_variant_calling == true) {
            GENE_FILTERING(CALL_VAR_M.out.rawGVCF,
                            CALL_VAR_F.out.rawGVCF,
                            CALL_VAR_C1.out.rawGVCF,
                            pedigrees)
        }
        else if(params.run_variant_calling == false) {
            GENE_FILTERING(mother_calls,
                father_calls,
                proband_calls,
                pedigrees)
        }
    }
}