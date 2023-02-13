#!/bin/bash
source settings_snv.sh
# Section1: uncomment this section and comment the section2 if the file formats FOLLOW our trio file name convention
#echo -e "Trio Child ID (e.g M_F071-003-A):"
#read CHILD_DIR
CHILD_DIR=$1
#TRIO_ID=$(echo $CHILD_DIR | awk -F "-" '{print $1}')
TRIO_ID=$(echo $CHILD_DIR | awk -F "-" '{print $1"-"$2}')
MOTHER_DIR=${TRIO_ID}-001-U
FATHER_DIR=${TRIO_ID}-002-U
#MOTHER_DIR=${TRIO_ID}-001_G-U
#FATHER_DIR=${TRIO_ID}-002_G-U
#CHILD_DIR=${TRIO_ID}-003-A

# End of Section1

# Section2: uncomment this section and comment the section1 if the file formats DO NOT FOLLOW our trio file name convention
#echo -e "Mother directory name:"
#read MOTHER_DIR
#echo -e "Father directory name:"
#read FATHER_DIR
#echo -e "Child directory name:"
#read CHILD_DIR

# End of section2

echo ${ALL_SAMPLES_OUTPUT_DIR}

OUTPUTS_CHILD="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}"
SAMPLE_ID_CHILD=$(ls ${OUTPUTS_CHILD}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}'| awk -F"/" '{print $(NF)}')
CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}"


OUTPUTS_MOTHER="${ALL_SAMPLES_OUTPUT_DIR}/${MOTHER_DIR}"
SAMPLE_ID_MOTHER=$(ls ${OUTPUTS_MOTHER}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}'| awk -F"/" '{print $(NF)}')
MOTHER=${OUTPUTS_MOTHER}/${SAMPLE_ID_MOTHER}


OUTPUTS_FATHER="${ALL_SAMPLES_OUTPUT_DIR}/${FATHER_DIR}"
SAMPLE_ID_FATHER=$(ls ${OUTPUTS_FATHER}/*.raw.g.vcf | awk -F".raw.g.vcf" '{print $1}' | awk -F"/" '{print $(NF)}')
FATHER=${OUTPUTS_FATHER}/${SAMPLE_ID_FATHER}

OUTPUTS_TRIO="${ALL_SAMPLES_OUTPUT_DIR}/${CHILD_DIR}_Trio"
TRIO="${OUTPUTS_TRIO}/${SAMPLE_ID_CHILD}"

BAM_CHILD="${OUTPUTS_CHILD}/${SAMPLE_ID_CHILD}.bqrecal.reads.bam"
BAM_MOTHER="${OUTPUTS_MOTHER}/${SAMPLE_ID_MOTHER}.bqrecal.reads.bam"
BAM_FATHER="${OUTPUTS_FATHER}/${SAMPLE_ID_FATHER}.bqrecal.reads.bam"

#CHILD_SAMPLE_NAME=$(echo $CHILD_DIR | awk -F "_" '{print $1"_"$2}')
#MOTHER_SAMPLE_NAME=$(echo $MOTHER_DIR | awk -F "_" '{print $1"_"$2}')
#FATHER_SAMPLE_NAME=$(echo $FATHER_DIR | awk -F "_" '{print $1"_"$2}')

CHILD_SAMPLE_NAME=$CHILD_DIR
MOTHER_SAMPLE_NAME=$MOTHER_DIR 
FATHER_SAMPLE_NAME=$FATHER_DIR

PEDIGREE="${ALL_SAMPLES_PEDIGREE_DIR}/pedigree.${CHILD_SAMPLE_NAME}"
SAMPLE_LIST=($CHILD_SAMPLE_NAME $MOTHER_SAMPLE_NAME $FATHER_SAMPLE_NAME)
SAMPLE_LIST_SORTED=( $(
for SAMP in "${SAMPLE_LIST[@]}"
do
     echo "$SAMP"
done | sort) )

echo "Sorted trio: ${SAMPLE_LIST_SORTED[@]}"

IDX=0
for SAMP in "${SAMPLE_LIST_SORTED[@]}"	
do
	if [ "$CHILD_SAMPLE_NAME" == "$SAMP" ] 
	then		
		let "CHILD_IDX=$IDX+9"
	elif [ "$MOTHER_SAMPLE_NAME" == "$SAMP" ]
	then
		let "MOTHER_IDX=$IDX+9"
	elif [ "$FATHER_SAMPLE_NAME" == "$SAMP" ]
	then
		let "FATHER_IDX=$IDX+9"
	else
		echo "Something wrong with sample names!!"
		exit
	fi
	let "IDX++"
done

export CHILD_IDX
export MOTHER_IDX
export FATHER_IDX



mkdir -p ${OUTPUTS_TRIO}
echo $REF

echo -e "\n\n*********************** Combine GVCs into one GVC file ***********************\n\n"
#$GATK CombineGVCFs -R $REF --variant ${CHILD}.raw.g.vcf --variant ${MOTHER}.raw.g.vcf --variant ${FATHER}.raw.g.vcf -O ${TRIO}.comb.trio.raw.vcf

echo -e "\n\n*********************** GenotypeGVCFs on trio *********************** \n\n"
#$GATK --java-options "-Xmx4g" GenotypeGVCFs -R $REF --variant ${TRIO}.comb.trio.raw.vcf -O ${TRIO}.trio.raw.vcf


echo -e "\n=============================================================================================\n"
echo -e " |                                      DENOVO SNPs                                            |\n"         
echo -e "=============================================================================================\n\n"

echo -e "\n\n *********************** Select SNPs from the raw vcf *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.trio.raw.vcf --select-type-to-include SNP -O ${TRIO}.trio.snps.vcf


echo -e "\n\n *********************** Hard filter trio SNPs using GATK VariantFiltration ***********************\n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.trio.snps.vcf -O ${TRIO}.snps.gatk.varfilt.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP Filter" 

$GATK VariantFiltration -R $REF -V ${TRIO}.trio.snps.vcf -O ${TRIO}.snps.gatk.varfilt.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 


echo -e "\n\n*********************** Hard filter the SNPs manually ***********************\n\n"
less ${TRIO}.snps.gatk.varfilt.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_DEN}); $F[$ENV{CHILD_IDX}]=~/,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_DEN}) || ($1 < $ENV{MAD_PROB_DEN}) || ($2!=0 && $1 >= $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC1_PROB_DEN}) || ($2!=0 && $1 < $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC2_PROB_DEN}); "$F[$ENV{MOTHER_IDX}]"=~/,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > $ENV{MAD_FRAC_PAR_DEN}); "$F[$ENV{FATHER_IDX}]"=~ /,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > "$ENV{MAD_FRAC_PAR_DEN}"); print' > ${TRIO}.snps.hardfilt.fordenovo.vcf 


echo -e "\n\n *********************** Annotate SNPs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.snps.hardfilt.fordenovo.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.snps.hardfilt.fordenovo -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,clinvar_20220320,gnomad211_exome,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n *********************** Index the vcf file of annotated deleterious SNPs ***********************\n\n"
$IGVTOOLS index ${TRIO}.snps.hardfilt.fordenovo.hg38_multianno.vcf

echo -e "\n\n*********************** Find Denovo SNPs by TrioDenovo ***********************\n\n"
less ${TRIO}.snps.hardfilt.fordenovo.hg38_multianno.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; print' > ${TRIO}.snps.denovo.td.vcf


echo -e "\n\n*********************** Filter denovo SNPs using GATK VariantFiltration to find rare variants *********************** \n\n"
less ${TRIO}.snps.hardfilt.fordenovo.hg38_multianno.vcf | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_DEN; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_DEN; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_DEN; print' > ${TRIO}.snps.hardfilt.fordenovo.rare.vcf


echo -e "\n\n*********************** Denovo SNP: Pass only coding regions *********************** \n\n"
less ${TRIO}.snps.hardfilt.fordenovo.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.snps.hardfilt.fordenovo.rare.coding.vcf

echo -e "\n\n*********************** Denovo SNP: Insersect trio denovo output with annotated rare variants  *********************** \n\n"
bedtools intersect -a ${TRIO}.snps.hardfilt.fordenovo.rare.coding.vcf -b ${TRIO}.snps.denovo.td.vcf -header > ${TRIO}.snps.denovo.td.intersect.vcf 




echo -e "\n=============================================================================================\n"
echo -e " |                                      DENOVO INDELs                                          |\n"         
echo -e "=============================================================================================\n\n"

echo -e "\n\n *********************** Select INDELs from the raw vcf *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.trio.raw.vcf --select-type-to-include INDEL -O ${TRIO}.trio.indels.vcf


echo -e "\n\n *********************** Hard filter trio INDELs using GATK VariantFiltration ***********************\n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.trio.indels.vcf -O ${TRIO}.indels.gatk.varfilt.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

echo -e "\n\n*********************** Hard filter the INDELs manually ***********************\n\n"
less ${TRIO}.indels.gatk.varfilt.vcf | perl -lane  'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_DEN}); $F[$ENV{CHILD_IDX}]=~/,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_DEN}) || ($1 < $ENV{MAD_PROB_DEN}) || ($2!=0 && $1 >= $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC1_PROB_DEN}) || ($2!=0 && $1 < $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC2_PROB_DEN}); "$F[$ENV{MOTHER_IDX}]"=~/,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > $ENV{MAD_FRAC_PAR_DEN}); "$F[$ENV{FATHER_IDX}]"=~ /,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > "$ENV{MAD_FRAC_PAR_DEN}"); print'  > ${TRIO}.indels.hardfilt.fordenovo.vcf


echo -e "\n\n *********************** Annotate INDELs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.indels.hardfilt.fordenovo.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.indels.hardfilt.fordenovo -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,clinvar_20220320,gnomad211_exome,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n *********************** Index the vcf file of annotated deleterious INDELs ***********************\n\n"
$IGVTOOLS index ${TRIO}.snps.hardfilt.fordenovo.hg38_multianno.vcf

echo -e "\n\n*********************** Find Denovo INDELs by TrioDenovo ***********************\n\n"
less ${TRIO}.indels.hardfilt.fordenovo.hg38_multianno.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; print' > ${TRIO}.indels.denovo.td.vcf


echo -e "\n\n*********************** Filter denovo INDELs using GATK VariantFiltration to find rare variants *********************** \n\n"
less ${TRIO}.indels.hardfilt.fordenovo.hg38_multianno.vcf | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_DEN; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_DEN; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_DEN; print' > ${TRIO}.indels.hardfilt.fordenovo.rare.vcf

echo -e "\n\n*********************** Denovo INDELs: Pass only coding regions *********************** \n\n"
less ${TRIO}.indels.hardfilt.fordenovo.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.indels.hardfilt.fordenovo.rare.coding.vcf

echo -e "\n\n*********************** Denovo INDELs: Insersect trio denovo output with annotated rare variants  *********************** \n\n"
bedtools intersect -a ${TRIO}.indels.hardfilt.fordenovo.rare.coding.vcf -b ${TRIO}.indels.denovo.td.vcf -header > ${TRIO}.indels.denovo.td.intersect.vcf 



echo -e "\n=============================================================================================\n"
echo -e " |                                 HOMOZYGOUS RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n *********************** Merge SNPs and INDELs *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.snps.gatk.varfilt.vcf I=${TRIO}.indels.gatk.varfilt.vcf O=${TRIO}.merged.variants.vcf
less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; print' > ${TRIO}.hom.vcf

less ${TRIO}.hom.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./  || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print'  > ${TRIO}.hom.filt.vcf


echo -e "\n\n *********************** Annotate HOMs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.hom.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.hom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,clinvar_20220320,gnomad211_exome,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n*********************** Filter homozygous variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.hom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.hom.metasvm.lof.vcf


echo -e "\n\n*********************** Filter homozygous variants to find CADD deletrious variants *********************** \n\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.hg38_multianno.vcf  -filter "vc.getAttribute('CADD13_PHRED').equals('.')" --filter-name "CADD" -filter "vc.getAttribute('CADD13_PHRED') < $CADD_THR" --filter-name "CADD" -O  ${TRIO}.hom.cadd.hardfilt.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.hom.cadd.hardfilt.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.hom.cadd.select.vcf

less ${TRIO}.hom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.hom.cadd.vcf


echo -e "\n\n*********************** Hom: Merge CADD and MetaSVM deleterious homozygous variants, and remove duplicate variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.hom.cadd.vcf I=${TRIO}.hom.metasvm.lof.vcf O=${TRIO}.hom.cadd.metasvm.vcf 

echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
echo -e "\n ***********************  filter homozygous recessive variants to find rare variants ***********************\n\n"
less ${TRIO}.hom.cadd.metasvm.vcf | uniq | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_REC; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_REC; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_REC; print' > ${TRIO}.hom.cadd.metasvm.rare.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                 COMPOUND HETEROZYGOUS                                       |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n*********************** Filter the merged VCF file for GT quality and coverage ***********************\n\n"
less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print' > ${TRIO}.merged.filt.vcf


echo -e "\n\n*********************** Annotate the merged filtered VCF file ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.merged.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.merged -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,clinvar_20220320,gnomad211_exome,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

echo -e "\n\n*********************** Remove the header of the VCF ***********************\n\n"
egrep -v "^##" ${TRIO}.merged.hg38_multianno.vcf > ${TRIO}.merged.noheader.vcf


echo -e "\n\n*********************** Select columns: CHR, POS, REF, and gene name ***********************\n\n"
vcf-query -f '%CHROM:%POS %REF %ALT  %INFO/{Gene.refGene}\n' ${TRIO}.merged.hg38_multianno.vcf > ${TRIO}.merged.genes.vcf
$GATK VariantsToTable -V ${TRIO}.merged.hg38_multianno.vcf -F CHROM -F POS -F REF -F ALT -F Gene.refGene -F HET -F HOM-REF -GF GT  -O ${TRIO}.merged.genes.table


echo -e "\n\n*********************** Move gene name column to the last column ***********************\n\n"
awk '{print $0"\t"$5}' ${TRIO}.merged.genes.table |cut --complement -f5 > ${TRIO}.merged.genes2.table


echo -e "\n\n*********************** Select the rows with 2 hets and 1 homozygous reference genotype ***********************\n\n"
awk  '$5 == 2 && $6 == 1' ${TRIO}.merged.genes2.table > ${TRIO}.merged.2hets.table 


echo -e "\n\n*********************** Add line number and select the rows with duplicate gene name ***********************\n\n"
cat ${TRIO}.merged.2hets.table  | nl | uniq -D -f10 > ${TRIO}.merged.dup.table


echo -e "\n\n*********************** Add 2 more columns indicating whether the variant is inherited from the father (F) or mother (M)***********************\n\n"
awk '{if ($8== $10) print $0,"\tM\t",$11"M"; else if ($9 == $10) print $0,"\tF\t",$11"F" }' ${TRIO}.merged.dup.table > ${TRIO}.dup.parent.table

echo -e "\n\n*********************** Count the number of consencutive duplicate gene+inheritance ***********************\n\n"
cat  ${TRIO}.dup.parent.table | uniq -c -f12>  ${TRIO}.dup.parent.uniq.table


echo -e "\n\n*********************** Move gene name column to the last column ***********************\n\n"
awk '{print $0"\t"$12}' ${TRIO}.dup.parent.uniq.table |cut --complement -f11 > ${TRIO}.dup.parent.uniq2.table


echo -e "\n\n*********************** Count number of rows  with consective duplicate gene names ***********************\n\n"
cat  ${TRIO}.dup.parent.uniq2.table | uniq -c -f13 | awk '{if ($1 >1) print $0}' >  ${TRIO}.dup.parent.uniq3.table


echo -e "\n\n======= Export the gene names with potentially compound het variants to a new file =========\n\n"
cat ${TRIO}.dup.parent.uniq3.table | rev | cut -f 1| rev | cut -d'\' -f 1 > ${TRIO}.comphet.genenames.table


echo -e "\n\n*********************** Select the rows in the initial table of duplicate gene names that have the final gene names ***********************\n\n"
for item in $(cat ${TRIO}.comphet.genenames.table); do grep "\s$item$" ${TRIO}.merged.dup.table ; done > ${TRIO}.comphet.lines.table

cat ${TRIO}.comphet.lines.table | cut -f 2,3 > ${TRIO}.positions.tsv

echo -e "\n\n*********************** Select the rows in the original annotated vcf file that have the derived chromosomal positions derived for the comp het variants ***********************\n\n"
vcftools --vcf ${TRIO}.merged.hg38_multianno.vcf --positions ${TRIO}.positions.tsv --recode  --recode-INFO-all --out ${TRIO}.comphet


more ${TRIO}.comphet.recode.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.comphet.metasvm.lof.vcf


echo -e "\n\n*********************** Filter compound het variants using GATK VariantFiltration to find CADD deletrious variants *********************** \n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.comphet.recode.vcf  -filter "(vc.getAttribute('CADD13_PHRED').equals('.') || vc.getAttribute('CADD13_PHRED') < $CADD_THR )" --filter-name "CADD" -O ${TRIO}.comphet.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.comphet.cadd.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.comphet.cadd.varfilt.select.vcf

cat ${TRIO}.comphet.cadd.varfilt.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream'> ${TRIO}.comphet.cadd.vcf


echo -e "\n\n*********************** Comphet: Merge CADD and MetaSVM deleterious homozygous variants, and remove duplicate variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.comphet.cadd.vcf I=${TRIO}.comphet.metasvm.lof.vcf O=${TRIO}.comphet.cadd.metasvm.vcf  


echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
echo -e "\n*********************** Filter compound heterozygous variants to find rare variants ***********************\n\n"
less ${TRIO}.comphet.cadd.metasvm.vcf | uniq | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac)  > $MAF_REC; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g)  > $MAF_REC; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_REC; print' > ${TRIO}.comphet.cadd.metasvm.rare.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                    X_LINKED RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; next if $F[0] ne "chrX";"$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2; print' > ${TRIO}.xlink.vcf

less ${TRIO}.xlink.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./  || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print' > ${TRIO}.xlink.filt.vcf

echo -e "\n\n*********************** Annotate x-linked recessive variants by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.xlink.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.xlink -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,clinvar_20220320,gnomad211_exome,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

echo -e "\n\n*********************** Filter x-linked variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.xlink.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deleteion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.xlink.metasvm.lof.vcf


echo -e "\n\n*********************** Filter x-linked variants to find CADD deletrious variants *********************** \n\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.hg38_multianno.vcf  -filter "vc.getAttribute('CADD13_PHRED').equals('.')" --filter-name "CADD" -filter "vc.getAttribute('CADD13_PHRED') < $CADD_THR" --filter-name "CADD" -O  ${TRIO}.xlink.cadd.hardfilt.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic x-linked variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.hardfilt.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.xlink.cadd.select.vcf

less ${TRIO}.xlink.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.xlink.cadd.vcf


echo -e "\n\n*********************** Hom: Merge CADD and MetaSVM deleterious x-linked variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.xlink.cadd.vcf I=${TRIO}.xlink.metasvm.lof.vcf O=${TRIO}.xlink.cadd.metasvm.vcf 


echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
echo -e "\n *********************** Filter x-linked variants to find rare variants ***********************\n\n"
less ${TRIO}.xlink.cadd.metasvm.vcf | uniq | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_REC; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_REC; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_REC; print' > ${TRIO}.xlink.cadd.metasvm.rare.vcf
