process HOMOZYGOUS_RECESSIVE_INCOMPLETE {
    label = 'medium'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.hom.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(rawVCF)
        path refGenome
        path refIndex
        path annovarRef

    output:
        tuple val(family), path("${family}.hom.cadd.metasvm.rare.vcf"), emit: homCaddMetaSVMRare

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Check the genotype of the proband for homozygous alternate
    less ${family}.trio.raw.vcf | perl -lane 'print if /#/; "\$F[9]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==0; print' > ${family}.hom.vcf

    # Apply filters on GT quality and coverage and make sure genotype is valid
    less ${family}.hom.vcf | perl -lane 'print if /#/; next if ("\$F[9]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_REC); "\$F[9]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_REC) || (\$3 < $params.GT_QUAL_REC); print'  > ${family}.hom.filt.vcf

    # Annotate HOMs by Annovar
    perl \$ANNOVAR/table_annovar.pl  ${family}.hom.filt.vcf  $annovarRef  -buildver hg38 -out ${family}.hom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Filter homozygous variants to find MetaSVM deleterious variants
    more ${family}.hom.hg38_multianno.vcf | { egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' || true; } | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' || true; } > ${family}.hom.metasvm.lof.vcf

    # Filter homozygous variants to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.hom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR" --filter-name "High_CADD" -O  ${family}.hom.cadd.varfilt.vcf

    # Exclude the non-exonic variants (intronic, UTR, synonymous, etc)
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.hom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.hom.cadd.select.vcf
    less ${family}.hom.cadd.select.vcf | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous' || true; } > ${family}.hom.cadd.coding.vcf

    # Merge CADD and MetaSVM deleterious homozygous variants and remove duplicate variants
    java -jar \$PICARD MergeVcfs I=${family}.hom.cadd.coding.vcf I=${family}.hom.metasvm.lof.vcf O=${family}.hom.cadd.metasvm.vcf 

    # Filter homozygous recessive variants to find rare variants and remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.hom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_REC)" --filter-name "Rare" -O  ${family}.hom.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.hom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.hom.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.hom.cadd.metasvm.rare.0.vcf > ${family}.hom.cadd.metasvm.rare.00.vcf

    # Remove homozygous variants in chrX if proband is male
    SEX=\$(cat $pedigree | awk '\$2 ~ /-003.*-A/ {print}' | cut -f5 )

    if [[ \$SEX == 1 ]]; then	
        awk '\$1 != "chrX" {print}' ${family}.hom.cadd.metasvm.rare.00.vcf > temp && mv temp ${family}.hom.cadd.metasvm.rare.vcf
    else
        mv ${family}.hom.cadd.metasvm.rare.00.vcf ${family}.hom.cadd.metasvm.rare.vcf
    fi
    """
}

process X_LINKED_RECESSIVE_INCOMPLETE {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.xlink.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(rawVCF)
        path refGenome
        path refIndex
        path annovarRef

    output:
        tuple val(family), path("${family}.xlink.cadd.metasvm.rare.vcf"), emit: xLinkCaddMetaSVMRareVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    # Check the genotype of the proband for xlink
    less ${family}.trio.raw.vcf | perl -lane 'print if /#/; next if \$F[0] ne "chrX";"\$F[9]" =~ /(\\d+)\\/(\\d+):/; next if \$1!=\$2 || \$1==0; print' > ${family}.xlink.vcf

    # Apply filters on mapping and GT quality, and coverage and make sure the genotype is valid
    less ${family}.xlink.vcf | perl -lane 'print if /#/; next if ("\$F[9]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_REC); "\$F[9]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_REC) || (\$3 < $params.GT_QUAL_REC); print' > ${family}.xlink.filt.vcf

    # Annotate x-linked recessive variants by Annovar
    perl \$ANNOVAR/table_annovar.pl  ${family}.xlink.filt.vcf  $annovarRef  -buildver hg38 -out ${family}.xlink -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

    # Filter x-linked variants to find MetaSVM deleterious variants
    more ${family}.xlink.hg38_multianno.vcf | { egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' || true; } | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' || true; } > ${family}.xlink.metasvm.lof.vcf

    # Filter x-linked variants to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.xlink.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR " --filter-name "High_CADD" -O  ${family}.xlink.cadd.varfilt.vcf

    # Exclude the non-exonic x-linked variants (intronic, UTR, synonymous, etc)
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.xlink.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.xlink.cadd.select.vcf
    less ${family}.xlink.cadd.select.vcf | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous' || true; } > ${family}.xlink.cadd.coding.vcf

    # Merge CADD and MetaSVM deleterious x-linked variants
    java -jar \$PICARD MergeVcfs I=${family}.xlink.cadd.coding.vcf I=${family}.xlink.metasvm.lof.vcf O=${family}.xlink.cadd.metasvm.vcf 

    # Filter x-linked variants to find rare variants and remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.xlink.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_REC)" --filter-name "Rare" -O  ${family}.xlink.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.xlink.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.xlink.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.xlink.cadd.metasvm.rare.0.vcf > ${family}.xlink.cadd.metasvm.rare.vcf
    """
}

process DOMINANT_INCOMPLETE {
    label = 'long'

    publishDir "${params.outTrioDir}/${family}", pattern: "*.dom.cadd.metasvm.rare.vcf", mode: 'copy'

    input:
        tuple val(family), path(pedigree), path(rawVCF)
        path refGenome
        path refIndex
        path annovarRef

    output:
        tuple val(family), path("${family}.dom.cadd.metasvm.rare.vcf"), emit: domCaddMetaSVMRareVCF

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """

    cat ${family}.trio.raw.vcf | perl -lane 'print if /#/;next if \$F[5]<30 ; "\$F[9]"=~/(\\d+)\\/(\\d+):/; next if \$1==\$2; print' >  ${family}.dom.vcf

    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.dom.vcf -filter "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Low_Qual" -O ${family}.dom.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.varfilt.vcf -select 'vc.isNotFiltered()' -O ${family}.dom.select.vcf

    less ${family}.dom.select.vcf | perl -lane 'print if /#/; next if ("\$F[9]" =~ /^\\./ || \$F[5] < $params.MAP_QUAL_DOM); "\$F[9]" =~ /,(\\d+):(\\d+):(\\d+)/; next if (\$2 < $params.COV_PROB_DOM) || (\$3 < $params.GT_QUAL_DOM) || (\$2!=0 && \$1 >= $params.MAD_LINE_DOM && \$1/\$2 < $params.MA_FRAC1_PROB_DOM) || (\$2!=0 && \$1 < $params.MAD_LINE_DOM && \$1/\$2 < $params.MA_FRAC2_PROB_DOM); print'  > ${family}.dom.filt.vcf

    # Annotate dominant variants by Annovar
    perl \$ANNOVAR/table_annovar.pl ${family}.dom.filt.vcf  $annovarRef  -buildver hg38 -out ${family}.dom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
    more ${family}.dom.hg38_multianno.vcf | { egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' || true; } | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous' || true; } > ${family}.dom.metasvm.lof.vcf

    # Filter dominant variants using GATK VariantFiltration to find CADD deleterious variants
    gatk --java-options "-Xmx${gatkMemory}g" VariantFiltration -R $refGenome -V ${family}.dom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $params.CADD_THR " --filter-name "High_CADD" -O  ${family}.dom.cadd.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${family}.dom.cadd.select.vcf
    less ${family}.dom.cadd.select.vcf | { egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous' || true; } > ${family}.dom.cadd.coding.vcf

    # Merge CADD and MetaSVM dominant variants
    java -jar \$PICARD MergeVcfs I=${family}.dom.cadd.coding.vcf I=${family}.dom.metasvm.lof.vcf O=${family}.dom.cadd.metasvm.vcf 

    # Filter dominant variants to find rare variants and remove duplicate variants
    gatk --java-options "-Xmx${gatkMemory}g"  VariantFiltration -R $refGenome -V ${family}.dom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $params.MAF_DOM) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $params.MAF_DOM)" --filter-name "Rare" -O  ${family}.dom.cadd.metasvm.varfilt.vcf
    gatk --java-options "-Xmx${gatkMemory}g" SelectVariants -R $refGenome -V ${family}.dom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${family}.dom.cadd.metasvm.rare.0.vcf
    uniq -f 7 ${family}.dom.cadd.metasvm.rare.0.vcf > ${family}.dom.cadd.metasvm.rare.vcf
    """
}

process VARIANTS_TO_TABLE_INCOMPLETE {
    label = 'short'

    publishDir "${params.outTrioDir}/${family}/Excel_Tables/", mode: 'copy'

    input:
        tuple val(family), path(hom_recessive_in), path(xlinked_in), path(dominant_in)

    output:
        path "${family}.hom.cadd.metasvm.rare.csv"
        path "${family}.xlink.cadd.metasvm.rare.csv"
        path "${family}.dom.cadd.metasvm.rare.csv"

    script:
    def gatkMemory = task.memory.toString().split()[0]

    """
    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $hom_recessive_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.hom.cadd.metasvm.rare.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable -V $xlinked_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.xlink.cadd.metasvm.rare.csv

    gatk --java-options "-Xmx${gatkMemory}g" VariantsToTable  -V $dominant_in --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O ${family}.dom.cadd.metasvm.rare.csv
    """
}
