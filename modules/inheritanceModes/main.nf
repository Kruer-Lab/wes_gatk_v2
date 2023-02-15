process DENOVO_SNPS {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.snp.denovo.intersect.vcf",mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(trioRawVCF)
        path refGenome
        path refIndex

    output:
        tuple val(family), path(pedigree), path("${family}.snp.denovo.intersect.vcf"), emit: denovoIntersectVCF
        tuple val(family), path(pedigree), path("${family}.snp.select.vcf"), emit: snpSelectVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Select SNPs from raw vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V $trioRawVCF --select-type-to-include SNP -O ${family}.trio.snp.vcf
    
    # Hard filter trio SNPs
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.trio.snp.vcf -O ${family}.snp.varfilt.vcf -filter "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Low_Qual"
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.snp.varfilt.vcf -select 'vc.isNotFiltered()' -O ${family}.snp.select.vcf

    # Hard filter SNPs manually
    less ${family}.snp.select.vcf | perl -lane 'print if /#/; next if ("\$F[$params.childIDX]" =~ /^\\./ || "\$F[$params.motherIDX]" =~ /^\\./ || "\$F[$params.fatherIDX]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_DEN); \$F[$params.childIDX]=~/,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_DEN) || (\$1 < $params.MAD_PROB_DEN) || (\$2!=0 && \$1 >= $params.MAD_LINE_DEN && \$1/\$2 < $params.MA_FRAC1_PROB_DEN) || (\$2!=0 && \$1 < $params.MAD_LINE_DEN && \$1/\$2 < $params.MA_FRAC2_PROB_DEN); "\$F[$params.motherIDX]"=~/,(\\d+):(\\d+):(\\d+)/;next if (\$2 < $params.COV_PAR_DEN) || (\$2 != 0 && \$1/\$2 > $params.MAD_FRAC_PAR_DEN); "\$F[$params.fatherIDX]"=~ /,(\\d+):(\\d+):(\\d+)/;next if (\$2 < $params.COV_PAR_DEN) || (\$2 != 0 && \$1/\$2 > $params.MAD_FRAC_PAR_DEN); print' > ${family}.snp.hardfilt.vcf 
    
    # Annotate SNPs by Annovar
    perl \$ANNOVAR/table_annovar.pl ${family}.snp.hardfilt.vcf ${params.resourcesDir}/Annovar/humandb/ -buildver hg38 -out ${family}.snp.hardfilt -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
    
    # Index VCF file of annotated deleterious SNPs
    igvtools index ${family}.snp.hardfilt.hg38_multianno.vcf

    # Find denovo SNPs
    less ${family}.snp.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; "\$F[$params.childIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1==\$2; "\$F[$params.motherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==1; "\$F[$params.fatherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==1; print' > ${family}.snp.denovo.vcf
    
    # Filter denovo SNPs using GATK VariantFiltration to find rare variants
    gatk --java-options "-Xmx${gatkMemory}g"  VariantFiltration -R $refGenome -V ${family}.snp.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_DEN) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_DEN)" --filter-name "Rare" -O  ${family}.snp.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.snp.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.snp.rare.vcf
 
    # Select coding variants
    less ${family}.snp.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${family}.snp.rare.coding.vcf

    # Intersect trio denovo output with annotated rare variants
    bedtools intersect -a ${family}.snp.rare.coding.vcf -b ${family}.snp.denovo.vcf -header > ${family}.snp.denovo.intersect.vcf
    """
}

process DENOVO_INDELS {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.indel.denovo.intersect.vcf", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(trioRawVCF)
        path refGenome
        path refIndex

    output:
        tuple val(family), path("${family}.indel.denovo.intersect.vcf"), emit: denovoIntersectVCF
        tuple val(family), path("${family}.indel.select.vcf"), emit: indelSelectVCF
        
    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Select indels from raw vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.trio.raw.vcf --select-type-to-include INDEL -O ${family}.trio.indel.vcf

    # Hard filter trio indels using GATK VariantFiltration
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.trio.indel.vcf -O ${family}.indel.varfilt.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.indel.varfilt.vcf -select 'vc.isNotFiltered()' -O ${family}.indel.select.vcf

    # Hard filter indels manually
    less ${family}.indel.select.vcf | perl -lane  'print if /#/; next if ("\$F[$params.childIDX]" =~ /^\\./ || "\$F[$params.motherIDX]" =~ /^\\./ || "\$F[$params.fatherIDX]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_DEN); \$F[$params.childIDX]=~/,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_DEN) || (\$1 < $params.MAD_PROB_DEN) || (\$2!=0 && \$1 >= $params.MAD_LINE_DEN && \$1/\$2 < $params.MA_FRAC1_PROB_DEN) || (\$2!=0 && \$1 < $params.MAD_LINE_DEN && \$1/\$2 < $params.MA_FRAC2_PROB_DEN); "\$F[$params.motherIDX]"=~/,(\\d+):(\\d+):(\\d+)/;next if (\$2 < $params.COV_PAR_DEN) || (\$2 != 0 && \$1/\$2 > $params.MAD_FRAC_PAR_DEN); "\$F[$params.fatherIDX]"=~ /,(\\d+):(\\d+):(\\d+)/;next if (\$2 < $params.COV_PAR_DEN) || (\$2 != 0 && \$1/\$2 > "$params.MAD_FRAC_PAR_DEN"); print'  > ${family}.indel.hardfilt.vcf

    # Annotate indels by Annovar
    perl \$ANNOVAR/table_annovar.pl  ${family}.indel.hardfilt.vcf  ${params.resourcesDir}/Annovar/humandb/  -buildver hg38 -out ${family}.indel.hardfilt -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Index the VCF file of annotated deleterious indels
    igvtools index ${family}.indel.hardfilt.hg38_multianno.vcf

    # Find denovo indels
    less ${family}.indel.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; "\$F[$params.childIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1==\$2; "\$F[$params.motherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==1; "\$F[$params.fatherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==1; print' > ${family}.indel.denovo.vcf

    # Filter denovo indels using GATK VariantFiltration to find rare variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.indel.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_DEN) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_DEN)" --filter-name "Rare" -O  ${family}.indel.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.indel.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.indel.rare.vcf

    # Select coding variants
    less ${family}.indel.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${family}.indel.rare.coding.vcf

    # Intersect trio denovo output with annotated rare variants
    bedtools intersect -a ${family}.indel.rare.coding.vcf -b ${family}.indel.denovo.vcf -header > ${family}.indel.denovo.intersect.vcf
    """
}

process HOMOZYGOUS_RECESSIVE {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.hom.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(snpSelect)
        tuple val(family), path(indelSelect)
        path refGenome
        path refIndex

    output:
        tuple val(family), path("${family}.hom.cadd.metasvm.rare.vcf"), emit: homCaddMetaSVMRare
        tuple val(family), path("${family}.merged.variants.vcf"), emit: mergedVariantsVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Merge SNPs and indels
    java -jar \$PICARD MergeVcfs I=$snpSelect I=$indelSelect O=${family}.merged.variants.vcf

    less ${family}.merged.variants.vcf | perl -lane 'print if /#/; "\$F[$params.childIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==0; "\$F[$params.motherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1==\$2 && \$1==1; "\$F[$params.fatherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1==\$2 && \$1==1; print' > ${family}.hom.vcf

    less ${family}.hom.vcf | perl -lane 'print if /#/; next if ("\$F[$params.childIDX]" =~ /^\\./ || "\$F[$params.motherIDX]" =~ /^\\./ || "\$F[$params.fatherIDX]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_REC); "\$F[$params.childIDX]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_REC) || (\$3 < $params.GT_QUAL_REC); print'  > ${family}.hom.filt.vcf

    # Annotate HOMs by Annovar
    perl \$ANNOVAR/table_annovar.pl  ${family}.hom.filt.vcf  ${params.resourcesDir}/Annovar/humandb/  -buildver hg38 -out ${family}.hom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Filter homozygous variants to find MetaSVM deleterious variants
    more ${family}.hom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${family}.hom.metasvm.lof.vcf

    # Filter homozygous variants to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g"  VariantFiltration -R $refGenome -V ${family}.hom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR" --filter-name "High_CADD" -O  ${family}.hom.cadd.varfilt.vcf

    # Exclude non-exonic variants (intronic, UTR, synonymous, etc)
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.hom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.hom.cadd.select.vcf

    less ${family}.hom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${family}.hom.cadd.coding.vcf

    # Merge CADD and MetaSVM deleterious homozygous variants and remove duplicate variants
    java -jar \$PICARD MergeVcfs I=${family}.hom.cadd.coding.vcf I=${family}.hom.metasvm.lof.vcf O=${family}.hom.cadd.metasvm.vcf

    # Filter homozygous recessive variants to find rare variants and remove duplicates
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.hom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_REC)" --filter-name "Rare" -O  ${family}.hom.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.hom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.hom.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.hom.cadd.metasvm.rare.0.vcf > ${family}.hom.cadd.metasvm.rare.00.vcf

    # Remove homozygous variants in chrX if the proband is male
    SEX=\$(cat $pedigree | awk '\$2 ~ /-003-A/ {print}' | cut -f5 )

    if [[ \$SEX == 1 ]]; then	
	    awk '\$1 != "chrX" {print}' ${family}.hom.cadd.metasvm.rare.00.vcf > temp && mv temp ${family}.hom.cadd.metasvm.rare.vcf
    else
        mv ${family}.hom.cadd.metasvm.rare.00.vcf ${family}.hom.cadd.metasvm.rare.vcf
    fi
    """
}

process COMPOUND_HETEROZYGOUS {
    label = 'long'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.comphet.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(mergedVariantsVCF)
        path refGenome
        path refIndex

    output:
        tuple val(family), path("${family}.comphet.cadd.metasvm.rare.vcf"), emit: comphetCaddMetaSVMRareVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Filter the merged VCF file for GT quality and coverage
    less $mergedVariantsVCF | perl -lane 'print if /#/; next if ("\$F[$params.childIDX]" =~ /^\\./ || "\$F[$params.motherIDX]" =~ /^\\./ || "\$F[$params.fatherIDX]" =~ /^\\./ || \$F[6] ne "PASS" || \$F[5] < $params.MAP_QUAL_REC); "\$F[$params.childIDX]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_REC) || (\$3 < $params.GT_QUAL_REC); print' > ${family}.merged.filt.vcf

    # Annotate the merged filtered VCF file
    perl \$ANNOVAR/table_annovar.pl  ${family}.merged.filt.vcf  ${params.resourcesDir}/Annovar/humandb/  -buildver hg38 -out ${family}.merged -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Select columns: CHR, POS, REF, and gene name
    vcf-query -f '%CHROM:%POS %REF %ALT %INFO/Gene.refGene\\n' ${family}.merged.hg38_multianno.vcf > ${family}.merged.genes.vcf
    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V ${family}.merged.hg38_multianno.vcf -F CHROM -F POS -F REF -F ALT -F Gene.refGene -F HET -F HOM-REF -GF GT  -O ${family}.merged.genes.table

    # Move gene name column to the last column
    awk '{print \$0"\\t"\$5}' ${family}.merged.genes.table | cut --complement -f5 > ${family}.merged.genes2.table

    # Select the rows with 2 hets and 1 homozygous reference genotype
    # Add line number and select the rows with duplicate gene name
    awk  '\$5 == 2 && \$6 == 1' ${family}.merged.genes2.table | nl | uniq -D -f10 > ${family}.merged.dup.table

    # Add 2 more columns indicating whether the variant is inherited from the father (F) or mother (M)
    # Count the number of consecutive duplicate gene+inheritance
    awk '{if (\$8== \$10) print \$0,"\\tM\\t",\$11"M"; else if (\$9 == \$10) print \$0,"\\tF\\t",\$11"F" }' ${family}.merged.dup.table | uniq -c -f12>  ${family}.dup.parent.uniq.table

    # Move gene name column to the last column
    # Count number of rows with consecutive duplicate gene names
    # Export the gene names with potentially compound het variants to a new file
    awk '{print \$0"\\t"\$12}' ${family}.dup.parent.uniq.table |cut --complement -f11 | uniq -c -f13 | awk '{if (\$1 >1) print \$0}' | rev | cut -f 1| rev | cut -d'\\' -f 1 > ${family}.comphet.genenames.table

    # Select the rows in the initial table of duplicate gene names that have the final gene names
    for item in \$(cat ${family}.comphet.genenames.table); do grep "\\s\$item\$" ${family}.merged.dup.table || true; done > ${family}.comphet.lines.table
    cat ${family}.comphet.lines.table | cut -f 2,3 > ${family}.positions.tsv

    # Select the rows in the original annotated VCF file that have the derived chromosomal positions derived for the comp het variants
    vcftools --vcf ${family}.merged.hg38_multianno.vcf --positions ${family}.positions.tsv --recode  --recode-INFO-all --out ${family}.comphet
    more ${family}.comphet.recode.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${family}.comphet.metasvm.lof.vcf

    # Filter compound het variants using GATK VariantFiltration to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.comphet.recode.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR " --filter-name "High_CADD" -O  ${family}.comphet.cadd.varfilt.vcf

    # Exclude the non-exonic variants (intronic, UTR, synonymous, etc)
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.comphet.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.comphet.cadd.varfilt.select.vcf
    cat ${family}.comphet.cadd.varfilt.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream'> ${family}.comphet.cadd.coding.vcf

    # Merge CADD and MetaSVM deleterious compound het variants
    java -jar \$PICARD MergeVcfs I=${family}.comphet.cadd.coding.vcf I=${family}.comphet.metasvm.lof.vcf O=${family}.comphet.cadd.metasvm.vcf

    # Filter compound heterozygous variants to find rare variants
    # Remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.comphet.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_REC)" --filter-name "Rare" -O  ${family}.comphet.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.comphet.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.comphet.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.comphet.cadd.metasvm.rare.0.vcf > ${family}.comphet.cadd.metasvm.rare.vcf
    """
}

process X_LINKED_RECESSIVE {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.xlink.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(mergedVariantsVCF)
        path refGenome
        path refIndex

    output:
        tuple val(family), path("${family}.xlink.cadd.metasvm.rare.vcf"), emit: xLinkCaddMetaSVMRareVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    less $mergedVariantsVCF | perl -lane 'print if /#/; next if \$F[0] ne "chrX";"\$F[$params.childIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==0; "\$F[$params.motherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1==\$2 && \$1==1; "\$F[$params.fatherIDX]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2; print' > ${family}.xlink.vcf

    less ${family}.xlink.vcf | perl -lane 'print if /#/; next if ("\$F[$params.childIDX]" =~ /^\\./ || "\$F[$params.motherIDX]" =~ /^\\./ || "\$F[$params.fatherIDX]" =~ /^\\./  || \$F[6] ne "PASS" || \$F[5] < $params.MAP_QUAL_REC); "\$F[$params.childIDX]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_REC) || (\$3 < $params.GT_QUAL_REC); print' > ${family}.xlink.filt.vcf

    # Annotate X-Linked recessive variants by Annovar
    perl \$ANNOVAR/table_annovar.pl  ${family}.xlink.filt.vcf  ${params.resourcesDir}/Annovar/humandb/  -buildver hg38 -out ${family}.xlink -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Filter X-Linked variants to find MetaSVM deleterious variants
    more ${family}.xlink.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${family}.xlink.metasvm.lof.vcf

    # Filter X-Linked variants to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.xlink.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR " --filter-name "High_CADD" -O  ${family}.xlink.cadd.varfilt.vcf

    # Exclude the non-exonic X-Linked variants (intronic, UTR, synonymous, etc)
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.xlink.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.xlink.cadd.select.vcf

    less ${family}.xlink.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${family}.xlink.cadd.coding.vcf

    # Merge CADD and MetaSVM deleterious X-Linked variants
    java -jar \$PICARD MergeVcfs I=${family}.xlink.cadd.coding.vcf I=${family}.xlink.metasvm.lof.vcf O=${family}.xlink.cadd.metasvm.vcf 

    # Filter X-Linked variants to find rare variants
    # Remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.xlink.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_REC)" --filter-name "Rare" -O  ${family}.xlink.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.xlink.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.xlink.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.xlink.cadd.metasvm.rare.0.vcf > ${family}.xlink.cadd.metasvm.rare.vcf
    """
}

process DOMINANT {
    label = 'long'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.dom.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(childRawGVCF)
        path refGenome
        path refIndex

    output:
        tuple val(family), path("${family}.dom.cadd.metasvm.rare.vcf"), emit: domCaddMetaSVMRareVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" GenotypeGVCFs -R $refGenome --variant $childRawGVCF -O ${family}.het.raw.vcf

    cat ${family}.het.raw.vcf | perl -lane 'print if /#/;next if \$F[5]<30 ; \$F[9]=~/(\\d+)\\/(\\d+):/; next if \$1==\$2; print' >  ${family}.dom.vcf

    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.dom.vcf -filter "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Low_Qual" -O ${family}.dom.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.varfilt.vcf -select 'vc.isNotFiltered()' -O ${family}.dom.select.vcf

    less ${family}.dom.select.vcf | perl -lane 'print if /#/; next if ("\$F[9]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_DOM); "\$F[9]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_DOM) || (\$3 < $params.GT_QUAL_DOM) || (\$2!=0 && \$1 >= $params.MAD_LINE_DOM && \$1/\$2 < $params.MA_FRAC1_PROB_DOM) || (\$2!=0 && \$1 < $params.MAD_LINE_DOM && \$1/\$2 < $params.MA_FRAC2_PROB_DOM); print'  > ${family}.dom.filt.vcf

    # Annotate dominant variants by Annovar
    perl \$ANNOVAR/table_annovar.pl ${family}.dom.filt.vcf  ${params.resourcesDir}/Annovar/humandb/  -buildver hg38 -out ${family}.dom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    more ${family}.dom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${family}.dom.metasvm.lof.vcf

    # Filter dominant variants using GATK VariantFiltration to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.dom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR " --filter-name "High_CADD" -O  ${family}.dom.cadd.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.dom.cadd.select.vcf

    less ${family}.dom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${family}.dom.cadd.coding.vcf

    # Echo Merge CADD and MetaSVM dominant variants
    java -jar \$PICARD MergeVcfs I=${family}.dom.cadd.coding.vcf I=${family}.dom.metasvm.lof.vcf O=${family}.dom.cadd.metasvm.vcf 

    # Filter dominant variants to find rare variants
    # Remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.dom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_DOM) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_DOM)" --filter-name "Rare" -O  ${family}.dom.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.dom.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.dom.cadd.metasvm.rare.0.vcf > ${family}.dom.cadd.metasvm.rare.vcf
    """
}

process VARIANTS_TO_TABLE {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}/Excel_Tables/", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(snp_in)
        tuple val(family), path(indel_in)
        tuple val(family), path(hom_recessive_in)
        tuple val(family), path(comp_het_in)
        tuple val(family), path(xlinked_in)
        tuple val(family), path(dominant_in)

    output:
        path "${family}.snp.denovo.intersect.csv"
        path "${family}.indel.denovo.intersect.csv"
        path "${family}.hom.cadd.metasvm.rare.csv"
        path "${family}.comphet.cadd.metasvm.rare.dup.csv"
        path "${family}.xlink.cadd.metasvm.rare.csv"
        path "${family}.dom.cadd.metasvm.rare.csv"

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $snp_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.snp.denovo.intersect.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $indel_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.indel.denovo.intersect.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $hom_recessive_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.hom.cadd.metasvm.rare.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $comp_het_in  --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.comphet.cadd.metasvm.rare.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $xlinked_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.xlink.cadd.metasvm.rare.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable  -V $dominant_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.dom.cadd.metasvm.rare.csv

    # Remove the single compound hets in file
    awk '{print \$3}' ${family}.comphet.cadd.metasvm.rare.csv | sort | uniq -d | grep -F -f - ${family}.comphet.cadd.metasvm.rare.csv > ${family}.comphet.cadd.metasvm.rare.dup.noheader.csv
    head -n1 ${family}.comphet.cadd.metasvm.rare.csv > ${family}.comphet.cadd.metasvm.rare.dup.csv
    cat ${family}.comphet.cadd.metasvm.rare.dup.noheader.csv >> ${family}.comphet.cadd.metasvm.rare.dup.csv
    """
}