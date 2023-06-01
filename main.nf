nextflow.enable.dsl=2

include { FASTQC; MULTIQC } from './modules/qc'
include { CALL_VARIANTS as CALL_VAR_M } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_F } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_C } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_OM } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_IP } from './subworkflows/call_variants'
include { CALL_VARIANTS as CALL_VAR_IO } from './subworkflows/call_variants'
include { GENE_FILTERING as GENE_FILTERING_COMPLETE } from './subworkflows/gene_filtering'
include { GENE_FILTERING as GENE_FILTERING_INCOMPLETE } from './subworkflows/gene_filtering'

workflow {

    // Regular expression patterns for getting sample and family IDs from file names
    def sampleRegexPattern = ~/F\d{3,}.*-\d{3}.*-[AU]/
    def familyRegexPattern = ~/F\d{3,}/

    // Validate and import samplesheet
    def sheet = file(params.samplesheet, checkIfExists: true).readLines()*.split('\t') // Import sample sheet as .tsv
    sheet.remove(0) // Remove header

    def completeTrios = []
    def incompleteTrios = []
    def sheetFamilies = []
    for(line : sheet) {

        // Get family IDs
        for(element : line) {
            if(!element.isEmpty()) {
                sheetFamilies.add((element =~ familyRegexPattern).findAll()[0])
                break
            }
        }

        // Add complete trios
        if(line.size() >= 3 && line[0].contains("-001") && line[1].contains("-002") && line[2].contains("-003")){
            completeTrios.add(line)
        }
        
        // Add incomplete trios
        else {
            incompleteTrios.add(line)
        }
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

    // Run FastQC and MultiQC
    if(params.run_qc == true) {
        
        // Get all fastq files in inDataDir
        Channel.fromPath("${params.inDataDir}/**.fastq.gz") | FASTQC | collect | MULTIQC
    }

    // Run variant calling (everything up to HaplotypeCaller)
    if(params.run_variant_calling == true) {

        // Import sample sheet and get reads
        // Get mother samples for complete trios
        Channel
            .fromList(completeTrios)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row[0] =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row[0]}/*_R1*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row[0]}/*_R2*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of file pairs in ${row.Mother}"
                }

            return tuple(familyID, sampleID, r1, r2)
        }
            .set{read_pairs_mother}

        // Get father samples for complete trios
        Channel
            .fromList(completeTrios)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row[1] =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row[1]}/*_R1*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row[1]}/*_R2*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of file pairs in ${row.Father}"
                }

                return tuple(familyID, sampleID, r1, r2)
            }
            .set{read_pairs_father}

        // Get affected child samples for complete trios
        Channel
            .fromList(completeTrios)
            .map{row ->
                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                def sampleID = (row[2] =~ sampleRegexPattern).findAll()[0]

                // Get family ID ("F" number)
                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                def r1 = file("$params.inDataDir/${row[2]}/*_R1*.fastq.gz", checkIfExists: true)
                def r2 = file("$params.inDataDir/${row[2]}/*_R2*.fastq.gz", checkIfExists: true)

                // Ensure only one file is used for each read
                if(r1.size() != 1 || r2.size != 1) {
                    error "Incorrect number of file pairs in ${row.Child_Affected}"
                }

                return tuple(familyID, sampleID, r1, r2)
            }
            .set{read_pairs_child_affected}

        // Get other samples in family (if any) for complete trios
        Channel
            .fromList(completeTrios)
            .flatMap{row ->
                if(row.size() == 4) {

                    def splitRow = row[3].split(',')*.trim()

                    def sampleList = []
                    for(member : splitRow) {

                        // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                        def sampleID = (member =~ sampleRegexPattern).findAll()[0]

                        // Get family ID ("F" number)
                        def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                        def r1 = file("$params.inDataDir/${member}/*_R1*.fastq.gz", checkIfExists: true)
                        def r2 = file("$params.inDataDir/${member}/*_R2*.fastq.gz", checkIfExists: true)

                        // Ensure only one file is used for each read
                        if(r1.size() != 1 || r2.size != 1) {
                            error "Incorrect number of file pairs in ${row[3]}"
                        }

                        sampleList.add(tuple(familyID, sampleID, r1, r2))
                    }
                    return sampleList
                }
            }
            .set{read_pairs_other_members}

        // Get incomplete trio samples
        Channel
            .fromList(incompleteTrios)
            .flatMap{row ->

                    def sampleList = []
                    for(member : row) {
                        // Check for multiple samples in element
                        if(member.contains(',')) {
                            def splitRow = member.split(',')*.trim()

                            for(splitMember : splitRow) {
                                
                                // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                                def sampleID = (member =~ sampleRegexPattern).findAll()[0]

                                // Get family ID ("F" number)
                                def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                                def r1 = file("$params.inDataDir/${splitMember}/*_R1*.fastq.gz", checkIfExists: true)
                                def r2 = file("$params.inDataDir/${splitMember}/*_R2*.fastq.gz", checkIfExists: true)

                                // Ensure only one file is used for each read
                                if(r1.size() != 1 || r2.size != 1) {
                                    error "Incorrect number of file pairs in ${splitMember}"
                                }

                                sampleList.add(tuple(familyID, sampleID, r1, r2))
                            }
                        }
                        else if(!member.isEmpty()){
                            // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                            def sampleID = (member =~ sampleRegexPattern).findAll()[0]

                            // Get family ID ("F" number)
                            def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                            def r1 = file("$params.inDataDir/${member}/*_R1*.fastq.gz", checkIfExists: true)
                            def r2 = file("$params.inDataDir/${member}/*_R2*.fastq.gz", checkIfExists: true)

                            // Ensure only one file is used for each read
                            if(r1.size() != 1 || r2.size != 1) {
                                error "Incorrect number of file pairs in ${member}"
                            }

                            sampleList.add(tuple(familyID, sampleID, r1, r2))
                        }
                    }
                    return sampleList
                }
            .branch{family, sample, read1, read2 ->
                proband: sample.contains("-003")
                other_members: true
            }
            .set{incompleteMembers}

        CALL_VAR_M(read_pairs_mother)
        CALL_VAR_F(read_pairs_father)
        CALL_VAR_C(read_pairs_child_affected)
        CALL_VAR_OM(read_pairs_other_members)
        CALL_VAR_IP(incompleteMembers.proband)
        CALL_VAR_IO(incompleteMembers.other_members)
    }

    // Skipping variant calling
    else if(params.run_variant_calling == false && params.run_variant_filtering == true) {

        // Get complete trio gvcf
        Channel
            .fromList(completeTrios)
            .flatMap{row ->

                def gvcfList = []
                for(member : row) {
                // Check for multiple samples in element
                    if(member.contains(',')) {
                        def splitRow = member.split(',')*.trim()

                        for(splitMember : splitRow) {
                                
                            // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                            def sampleID = (member =~ sampleRegexPattern).findAll()[0]

                            // Get family ID ("F" number)
                            def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                            def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                            gvcfList.add(tuple(familyID, sampleID, gvcf))
                        }
                    }
                    else if(!member.isEmpty()){
                        // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                        def sampleID = (member =~ sampleRegexPattern).findAll()[0]

                        // Get family ID ("F" number)
                        def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                        def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                        gvcfList.add(tuple(familyID, sampleID, gvcf))
                    }
                }
                return gvcfList
            }
            .branch{family, sample, gvcf ->
                mother: sample.contains("-001")
                father: sample.contains("-002")
                proband: sample.contains("-003")
            }
            .set{completeGVCFs}

        // Get incomplete trio gvcf
        Channel
            .fromList(incompleteTrios)
            .map{row ->
                if(!row[2].isEmpty()) {
                    // Get sample ID (letter, 3 or more digits, -, 3 digits, -, letter)
                    def sampleID = (row[2] =~ sampleRegexPattern).findAll()[0]

                    // Get family ID ("F" number)
                    def familyID = (sampleID =~ familyRegexPattern).findAll()[0]

                    def gvcf = file("$params.outDataDir/$sampleID/${sampleID}.raw.g.vcf", checkIfExists: true)

                    return tuple(familyID, gvcf)
                }
            }
            .set{incompleteGVCFs}
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
            GENE_FILTERING_COMPLETE(CALL_VAR_M.out.rawGVCF,
                            CALL_VAR_F.out.rawGVCF,
                            CALL_VAR_C.out.rawGVCF,
                            pedigrees)

            GENE_FILTERING_INCOMPLETE(null,
                            null,
                            CALL_VAR_IP.out.rawGVCF,
                            pedigrees)
        }
        else if(params.run_variant_calling == false) {
            GENE_FILTERING_COMPLETE(completeGVCFs.mother,
                completeGVCFs.father,
                completeGVCFs.proband,
                pedigrees)

            GENE_FILTERING_INCOMPLETE(null,
                            null,
                            incompleteGVCFs,
                            pedigrees)
        }
    }
}
