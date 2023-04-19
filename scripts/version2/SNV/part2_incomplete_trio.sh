#!/bin/bash
source settings_snv.sh
#echo -e "Trio Child ID (e.g M_F071-003-A):"
CHILD_DIR=$1
TRIO_ID=$(echo $CHILD_DIR | awk -F "-" '{print $1"-"$2}')

echo ${ALL_SAMPLES_OUTPUT_DIR}

OUTPUTS_CHILD="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}"
SAMPLE_ID_CHILD=$(ls ${OUTPUTS_CHILD}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}'| awk -F"/" '{print $(NF)}')
CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}"
OUTPUTS_TRIO="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}_Trio"
TRIO="${OUTPUTS_TRIO}/${SAMPLE_ID_CHILD}"
BAM_CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}.bqrecal.reads.bam"
CHILD_SAMPLE_NAME=$CHILD_DIR
PEDIGREE="pedigree.${CHILD_SAMPLE_NAME}.p2.ped"

mkdir -p ${OUTPUTS_TRIO}



echo -e "\n\n*********************** Combine GVCs into one GVC file ***********************\n\n"
$GATK CombineGVCFs -R $REF --variant ${CHILD}.raw.g.vcf -O ${TRIO}.comb.trio.raw.vcf

echo -e "\n\n*********************** GenotypeGVCFs on trio *********************** \n\n"
$GATK --java-options "-Xmx4g" GenotypeGVCFs -R $REF --variant ${TRIO}.comb.trio.raw.vcf -O ${TRIO}.trio.raw.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                 HOMOZYGOUS RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n *********************** Check the genotype of the proband for homozygous alternate ***********************\n\n"
less ${TRIO}.trio.raw.vcf | perl -lane 'print if /#/; "$F[9]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; print' > ${TRIO}.hom.vcf

echo -e "\n\n *********************** Apply filters on GT quality and coverage and make sure the genotype is valid (not .) ***********************\n\n"
less ${TRIO}.hom.vcf | perl -lane 'print if /#/; next if ("$F[9]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_REC}); "$F[9]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print'  > ${TRIO}.hom.filt.vcf


echo -e "\n\n *********************** Annotate HOMs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.hom.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.hom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n*********************** Filter homozygous variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.hom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.hom.metasvm.lof.vcf


echo -e "\n\n*********************** Filter homozygous variants to find CADD deletrious variants *********************** \n\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR" --filter-name "High_CADD" -O  ${TRIO}.hom.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.hom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.hom.cadd.select.vcf

less ${TRIO}.hom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.hom.cadd.coding.vcf


echo -e "\n\n*********************** Hom: Merge CADD and MetaSVM deleterious homozygous variants, and remove duplicate variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.hom.cadd.coding.vcf I=${TRIO}.hom.metasvm.lof.vcf O=${TRIO}.hom.cadd.metasvm.vcf 
 

echo -e "\n ***********************  filter homozygous recessive variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.hom.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.hom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.hom.cadd.metasvm.rare.0.vcf
uniq -f 7 ${TRIO}.hom.cadd.metasvm.rare.0.vcf > ${TRIO}.hom.cadd.metasvm.rare.00.vcf

echo -e "\n *********************** Hom: remove homozygous variants in chrX if the proband is a male  ***********************\n\n" 
SEX=$(cat ${PEDIGREES_DIR}/$PEDIGREE | awk -v samp="$CHILD_SAMPLE_NAME" '$2 == samp {print}' | cut -f5 )

echo -e "proband sex is: $SEX "
if [[ $SEX == 1 ]]; then	
	awk '$1 != "chrX" {print}' ${TRIO}.hom.cadd.metasvm.rare.00.vcf > temp && mv temp ${TRIO}.hom.cadd.metasvm.rare.vcf
fi



echo -e "\n=============================================================================================\n"
echo -e " |                                    X_LINKED RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n *********************** Check the genotype of the proband for xlink  ***********************\n\n"
less ${TRIO}.trio.raw.vcf | perl -lane 'print if /#/; next if $F[0] ne "chrX";"$F[9]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; print' > ${TRIO}.xlink.vcf

echo -e "\n\n *********************** Apply filters on mapping and GT quality, and coverage and make sure the genotype is valid (not .) ***********************\n\n"
less ${TRIO}.xlink.vcf | perl -lane 'print if /#/; next if ("$F[9]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_REC}); "$F[9]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print' > ${TRIO}.xlink.filt.vcf


less ${TRIO}.hom.vcf | perl -lane 'print if /#/; next if ("$F[9]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_REC}); "$F[9]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print'  > ${TRIO}.hom.filt.vcf


echo -e "\n\n*********************** Annotate x-linked recessive variants by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.xlink.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.xlink -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

echo -e "\n\n*********************** Filter x-linked variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.xlink.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.xlink.metasvm.lof.vcf


echo -e "\n\n*********************** Filter x-linked variants to find CADD deletrious variants *********************** \n\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR " --filter-name "High_CADD" -O  ${TRIO}.xlink.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic x-linked variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.xlink.cadd.select.vcf

less ${TRIO}.xlink.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.xlink.cadd.coding.vcf


echo -e "\n\n*********************** Xlink: Merge CADD and MetaSVM deleterious x-linked variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.xlink.cadd.coding.vcf I=${TRIO}.xlink.metasvm.lof.vcf O=${TRIO}.xlink.cadd.metasvm.vcf 


echo -e "\n *********************** Filter x-linked variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.xlink.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.xlink.cadd.metasvm.rare.0.vcf
uniq -f 7 ${TRIO}.xlink.cadd.metasvm.rare.0.vcf > ${TRIO}.xlink.cadd.metasvm.rare.vcf





echo -e "\n=============================================================================================\n"
echo -e " |                                         DOMINANT                                            |\n"         
echo -e "\n=============================================================================================\n\n"

cat ${TRIO}.trio.raw.vcf | perl -lane 'print if /#/;next if $F[5]<30 ; "$F[9]"=~/(\d+)\/(\d+):/; next if $1==$2; print' >  ${TRIO}.dom.vcf


$GATK VariantFiltration -R $REF -V ${TRIO}.dom.vcf -filter "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Low_Qual" -O ${TRIO}.dom.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.dom.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.dom.select.vcf

less ${TRIO}.dom.select.vcf | perl -lane 'print if /#/; next if ("$F[9]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_DOM}); "$F[9]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_DOM}) || ($3 < $ENV{GT_QUAL_DOM}) || ($2!=0 && $1 >= $ENV{MAD_LINE_DOM} && $1/$2 < $ENV{MA_FRAC1_PROB_DOM}) || ($2!=0 && $1 < $ENV{MAD_LINE_DOM} && $1/$2 < $ENV{MA_FRAC2_PROB_DOM}); print'  > ${TRIO}.dom.filt.vcf

echo -e "\n\n*********************** Annotate dominant variants by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl ${TRIO}.dom.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.dom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

more ${TRIO}.dom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.dom.metasvm.lof.vcf

echo -e "\n\n*********************** Filter dominant variants using GATK VariantFiltration to find CADD deletrious variants *********************** \n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.dom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR " --filter-name "High_CADD" -O  ${TRIO}.dom.cadd.varfilt.vcf

$GATK SelectVariants -R $REF -V ${TRIO}.dom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.dom.cadd.select.vcf

less ${TRIO}.dom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.dom.cadd.coding.vcf

echo -e "\n\n*********************** Dominant: Merge CADD and MetaSVM deleterious x-linked variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.dom.cadd.coding.vcf I=${TRIO}.dom.metasvm.lof.vcf O=${TRIO}.dom.cadd.metasvm.vcf 

echo -e "\n *********************** Filter dominant variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.dom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DOM) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_DOM)" --filter-name "Rare" -O  ${TRIO}.dom.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.dom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.dom.cadd.metasvm.rare.0.vcf

uniq -f 7 ${TRIO}.dom.cadd.metasvm.rare.0.vcf > ${TRIO}.dom.cadd.metasvm.rare.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                         WRAP UP                                            |\n"         
echo -e "\n=============================================================================================\n\n"
echo -e  "\n\n*********************** Move final results to a folder ***********************\n\n"

mkdir -p ${OUTPUTS_TRIO}/Final_Results
mkdir -p ${OUTPUTS_TRIO}/Final_Results/Excel_Tables

cp  ${TRIO}.hom.cadd.metasvm.rare.vcf ${TRIO}.dom.cadd.metasvm.rare.vcf ${TRIO}.xlink.cadd.metasvm.rare.vcf ${OUTPUTS_TRIO}/Final_Results

echo -e "\n=============================================================================================\n"
echo -e " |                                   Generate EXCEL Reports                                    |\n"         
echo -e "\n=============================================================================================\n\n"

 
HOM_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.hom.cadd.metasvm.rare.vcf
XLINK_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.xlink.cadd.metasvm.rare.vcf
DOMINANT_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.dom.cadd.metasvm.rare.vcf

HOM_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.hom.cadd.metasvm.rare.csv
XLINK_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.xlink.cadd.metasvm.rare.csv
DOMINANT_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.dom.cadd.metasvm.rare.csv



$GATK VariantsToTable -V $HOM_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $HOM_OUT

$GATK VariantsToTable -V $XLINK_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $XLINK_OUT

$GATK VariantsToTable  -V $DOMINANT_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $DOMINANT_OUT


