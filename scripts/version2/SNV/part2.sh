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

#PEDIGREE="${PEDIGREES_DIR}/pedigree.${CHILD_SAMPLE_NAME}.p2.ped"
PEDIGREE="pedigree.${CHILD_SAMPLE_NAME}.p2.ped"
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
$GATK SelectVariants -R $REF -V ${TRIO}.trio.raw.vcf --select-type-to-include SNP -O ${TRIO}.trio.snp.vcf


echo -e "\n\n *********************** Hard filter trio SNPs using GATK VariantFiltration ***********************\n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.trio.snp.vcf -O ${TRIO}.snp.varfilt.vcf -filter "QD < 2.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Low_Qual" 
$GATK SelectVariants -R $REF -V ${TRIO}.snp.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.snp.select.vcf


echo -e "\n\n*********************** Hard filter the SNPs manually ***********************\n\n"
less ${TRIO}.snp.select.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_DEN}); $F[$ENV{CHILD_IDX}]=~/,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_DEN}) || ($1 < $ENV{MAD_PROB_DEN}) || ($2!=0 && $1 >= $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC1_PROB_DEN}) || ($2!=0 && $1 < $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC2_PROB_DEN}); "$F[$ENV{MOTHER_IDX}]"=~/,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > $ENV{MAD_FRAC_PAR_DEN}); "$F[$ENV{FATHER_IDX}]"=~ /,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > "$ENV{MAD_FRAC_PAR_DEN}"); print' > ${TRIO}.snp.hardfilt.vcf 


echo -e "\n\n *********************** Annotate SNPs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.snp.hardfilt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.snp.hardfilt -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n *********************** Index the vcf file of annotated deleterious SNPs ***********************\n\n"
$IGVTOOLS index ${TRIO}.snp.hardfilt.hg38_multianno.vcf

echo -e "\n\n*********************** Find Denovo SNPs ***********************\n\n"
#less ${TRIO}.snp.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; print' > ${TRIO}.snp.denovo.vcf
docker run -it -v ${PEDIGREES_DIR}:/Pedigree_Dir_Dock  -v ${OUTPUTS_TRIO}:/Output_Dir_Dock spashleyfu/ubuntu20_triodenovo:0.0.6  triodenovo --ped /Pedigree_Dir_Dock/$PEDIGREE --in_vcf /Output_Dir_Dock/${SAMPLE_ID_CHILD}.snp.hardfilt.vcf --out /Output_Dir_Dock/${CHILD_SAMPLE_NAME}.snp.denovo.vcf


echo -e "\n\n*********************** Filter denovo SNPs using GATK VariantFiltration to find rare variants *********************** \n\n"
#less ${TRIO}.snp.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; $F[7] =~ /bravo_freeze8=(\d*).(\d*);/; $maf_bravo=$1.".".$2; next if int($maf_bravo) > $MAF_DEN; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_DEN; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_DEN; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_DEN; print' > ${TRIO}.snp.hardfilt.rare.vcf

####$GATK  VariantFiltration -R $REF -V ${TRIO}.snp.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DEN) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_DEN) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_DEN) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_DEN)" --filter-name "Rare" -O  ${TRIO}.snp.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.snp.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DEN) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_DEN)" --filter-name "Rare" -O  ${TRIO}.snp.varfilt.vcf

$GATK SelectVariants -R $REF -V ${TRIO}.snp.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.snp.rare.vcf


echo -e "\n\n*********************** Denovo SNP: Select coding variants *********************** \n\n"
less ${TRIO}.snp.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.snp.rare.coding.vcf

echo -e "\n\n*********************** Denovo SNP: Insersect trio denovo output with annotated rare variants  *********************** \n\n"
bedtools intersect -a ${TRIO}.snp.rare.coding.vcf -b ${TRIO}.snp.denovo.vcf -header > ${TRIO}.snp.denovo.intersect.vcf 




echo -e "\n=============================================================================================\n"
echo -e " |                                      DENOVO INDELs                                          |\n"         
echo -e "=============================================================================================\n\n"

echo -e "\n\n *********************** Select INDELs from the raw vcf *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.trio.raw.vcf --select-type-to-include INDEL -O ${TRIO}.trio.indel.vcf


echo -e "\n\n *********************** Hard filter trio INDELs using GATK VariantFiltration ***********************\n\n"
$GATK VariantFiltration -R $REF -V ${TRIO}.trio.indel.vcf -O ${TRIO}.indel.varfilt.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
$GATK SelectVariants -R $REF -V ${TRIO}.indel.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.indel.select.vcf


echo -e "\n\n*********************** Hard filter the INDELs manually ***********************\n\n"
less ${TRIO}.indel.select.vcf | perl -lane  'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_DEN}); $F[$ENV{CHILD_IDX}]=~/,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_DEN}) || ($1 < $ENV{MAD_PROB_DEN}) || ($2!=0 && $1 >= $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC1_PROB_DEN}) || ($2!=0 && $1 < $ENV{MAD_LINE_DEN} && $1/$2 < $ENV{MA_FRAC2_PROB_DEN}); "$F[$ENV{MOTHER_IDX}]"=~/,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > $ENV{MAD_FRAC_PAR_DEN}); "$F[$ENV{FATHER_IDX}]"=~ /,(\d+):(\d+):(\d+)/;next if ($2 < $ENV{COV_PAR_DEN}) || ($2 != 0 && $1/$2 > "$ENV{MAD_FRAC_PAR_DEN}"); print'  > ${TRIO}.indel.hardfilt.vcf


echo -e "\n\n *********************** Annotate INDELs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.indel.hardfilt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.indel.hardfilt -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n *********************** Index the vcf file of annotated deleterious INDELs ***********************\n\n"
$IGVTOOLS index ${TRIO}.indel.hardfilt.hg38_multianno.vcf

echo -e "\n\n*********************** Find Denovo INDELs  ***********************\n\n"
#less ${TRIO}.indel.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==1; print' > ${TRIO}.indel.denovo.vcf
docker run -it -v ${PEDIGREES_DIR}:/Pedigree_Dir_Dock  -v ${OUTPUTS_TRIO}:/Output_Dir_Dock spashleyfu/ubuntu20_triodenovo:0.0.6  triodenovo --ped /Pedigree_Dir_Dock/$PEDIGREE --in_vcf /Output_Dir_Dock/${SAMPLE_ID_CHILD}.indel.hardfilt.vcf --out /Output_Dir_Dock/${CHILD_SAMPLE_NAME}.indel.denovo.vcf


echo -e "\n\n*********************** Filter denovo INDELs using GATK VariantFiltration to find rare variants *********************** \n\n"
#less ${TRIO}.indel.hardfilt.hg38_multianno.vcf | perl -lane 'print if /#/; $F[7] =~ /bravo_freeze8=(\d*).(\d*);/; $maf_bravo=$1.".".$2; next if int($maf_bravo) > $MAF_DEN; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_DEN; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_DEN; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_DEN; print' > ${TRIO}.indel.hardfilt.rare.vcf

###### $GATK  VariantFiltration -R $REF -V ${TRIO}.indel.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DEN) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_DEN) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_DEN) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_DEN)" --filter-name "Rare" -O  ${TRIO}.indel.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.indel.hardfilt.hg38_multianno.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DEN) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_DEN)" --filter-name "Rare" -O  ${TRIO}.indel.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.indel.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.indel.rare.vcf


echo -e "\n\n*********************** Denovo INDELs: Select coding variants *********************** \n\n"
less ${TRIO}.indel.rare.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.indel.rare.coding.vcf

echo -e "\n\n*********************** Denovo INDELs: Insersect trio denovo output with annotated rare variants  *********************** \n\n"
bedtools intersect -a ${TRIO}.indel.rare.coding.vcf -b ${TRIO}.indel.denovo.vcf -header > ${TRIO}.indel.denovo.intersect.vcf 



echo -e "\n=============================================================================================\n"
echo -e " |                                 HOMOZYGOUS RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n *********************** Merge SNPs and INDELs *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.snp.select.vcf I=${TRIO}.indel.select.vcf O=${TRIO}.merged.variants.vcf
less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; "$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; print' > ${TRIO}.hom.vcf

less ${TRIO}.hom.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print'  > ${TRIO}.hom.filt.vcf


echo -e "\n\n *********************** Annotate HOMs by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.hom.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.hom -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n*********************** Filter homozygous variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.hom.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.hom.metasvm.lof.vcf


echo -e "\n\n*********************** Filter homozygous variants to find CADD deletrious variants *********************** \n\n"
########$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.hg38_multianno.vcf  -filter "vc.getAttribute('cadd16_phred').equals('.')" --filter-name "CADD" -filter "vc.getAttribute('cadd16_phred') < $CADD_THR" --filter-name "CADD" -O  ${TRIO}.hom.cadd.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR" --filter-name "High_CADD" -O  ${TRIO}.hom.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.hom.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.hom.cadd.select.vcf

less ${TRIO}.hom.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.hom.cadd.coding.vcf


echo -e "\n\n*********************** Hom: Merge CADD and MetaSVM deleterious homozygous variants, and remove duplicate variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.hom.cadd.coding.vcf I=${TRIO}.hom.metasvm.lof.vcf O=${TRIO}.hom.cadd.metasvm.vcf 
 

echo -e "\n ***********************  filter homozygous recessive variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
#less ${TRIO}.hom.cadd.metasvm.vcf | uniq | perl -lane 'print if /#/; $F[7] =~ /ExAC_ALL=(\d*).(\d*);/; $maf_exac=$1.".".$2; next if int($maf_exac) > $MAF_REC; $F[7] =~ /ALL.sites.2015_08=(\d*).(\d*);/; $maf_1000g=$1.".".$2; next if int($maf_1000g) > $MAF_REC; $F[7] =~ /esp6500siv2_all=(\d*).(\d*);/; $maf_esp=$1.".".$2; next if int($maf_esp) > $MAF_REC; print' > ${TRIO}.hom.cadd.metasvm.rare.vcf
###### $GATK  VariantFiltration -R $REF -V ${TRIO}.hom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_REC) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_REC) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.hom.cadd.metasvm.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.hom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.hom.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.hom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.hom.cadd.metasvm.rare.0.vcf
uniq -f 7 ${TRIO}.hom.cadd.metasvm.rare.0.vcf > ${TRIO}.hom.cadd.metasvm.rare.00.vcf

echo -e "\n *********************** Hom: remove homozygous variants in chrX if the proband is a male  ***********************\n\n" 
#SEX=$(cat $PEDIGREE | awk -v samp="$CHILD_SAMPLE_NAME" '$2 == samp {print}' | cut -f5 )
SEX=$(cat ${PEDIGREES_DIR}/$PEDIGREE | awk -v samp="$CHILD_SAMPLE_NAME" '$2 == samp {print}' | cut -f5 )

echo -e "proband sex is: $SEX "
if [[ $SEX == 1 ]]; then	
	awk '$1 != "chrX" {print}' ${TRIO}.hom.cadd.metasvm.rare.00.vcf > temp && mv temp ${TRIO}.hom.cadd.metasvm.rare.vcf
fi


echo -e "\n=============================================================================================\n"
echo -e " |                                 COMPOUND HETEROZYGOUS                                       |\n"         
echo -e "\n=============================================================================================\n\n"

echo -e "\n\n*********************** Filter the merged VCF file for GT quality and coverage ***********************\n\n"
less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./ || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print' > ${TRIO}.merged.filt.vcf


echo -e "\n\n*********************** Annotate the merged filtered VCF file ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.merged.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.merged -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput


echo -e "\n\n*********************** Select columns: CHR, POS, REF, and gene name ***********************\n\n"
vcf-query -f '%CHROM:%POS %REF %ALT  %INFO/{Gene.refGene}\n' ${TRIO}.merged.hg38_multianno.vcf > ${TRIO}.merged.genes.vcf
$GATK VariantsToTable -V ${TRIO}.merged.hg38_multianno.vcf -F CHROM -F POS -F REF -F ALT -F Gene.refGene -F HET -F HOM-REF -GF GT  -O ${TRIO}.merged.genes.table


echo -e "\n\n*********************** Move gene name column to the last column ***********************\n\n"
awk '{print $0"\t"$5}' ${TRIO}.merged.genes.table | cut --complement -f5 > ${TRIO}.merged.genes2.table


echo -e "\n\n*********************** Select the rows with 2 hets and 1 homozygous reference genotype ***********************\n"
echo -e "\n\n*********************** Add line number and select the rows with duplicate gene name ***********************\n\n"
awk  '$5 == 2 && $6 == 1' ${TRIO}.merged.genes2.table | nl | uniq -D -f10 > ${TRIO}.merged.dup.table


echo -e "\n\n*********************** Add 2 more columns indicating whether the variant is inherited from the father (F) or mother (M)***********************\n"
echo -e "\n\n*********************** Count the number of consencutive duplicate gene+inheritance ***********************\n\n"
awk '{if ($8== $10) print $0,"\tM\t",$11"M"; else if ($9 == $10) print $0,"\tF\t",$11"F" }' ${TRIO}.merged.dup.table | uniq -c -f12>  ${TRIO}.dup.parent.uniq.table


echo -e "\n\n*********************** Move gene name column to the last column ***********************\n"
echo -e "\n\n*********************** Count number of rows  with consective duplicate gene names ***********************\n"
echo -e "\n\n======= Export the gene names with potentially compound het variants to a new file =========\n\n"

awk '{print $0"\t"$12}' ${TRIO}.dup.parent.uniq.table |cut --complement -f11 | uniq -c -f13 | awk '{if ($1 >1) print $0}' | rev | cut -f 1| rev | cut -d'\' -f 1 > ${TRIO}.comphet.genenames.table


echo -e "\n\n*********************** Select the rows in the initial table of duplicate gene names that have the final gene names ***********************\n\n"
for item in $(cat ${TRIO}.comphet.genenames.table); do grep "\s$item$" ${TRIO}.merged.dup.table ; done > ${TRIO}.comphet.lines.table

cat ${TRIO}.comphet.lines.table | cut -f 2,3 > ${TRIO}.positions.tsv

echo -e "\n\n*********************** Select the rows in the original annotated vcf file that have the derived chromosomal positions derived for the comp het variants ***********************\n\n"
vcftools --vcf ${TRIO}.merged.hg38_multianno.vcf --positions ${TRIO}.positions.tsv --recode  --recode-INFO-all --out ${TRIO}.comphet





more ${TRIO}.comphet.recode.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.comphet.metasvm.lof.vcf


echo -e "\n\n*********************** Filter compound het variants using GATK VariantFiltration to find CADD deletrious variants *********************** \n\n"
$GATK  VariantFiltration -R $REF -V ${TRIO}.comphet.recode.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR " --filter-name "High_CADD" -O  ${TRIO}.comphet.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
$GATK SelectVariants -R $REF -V ${TRIO}.comphet.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.comphet.cadd.varfilt.select.vcf
cat ${TRIO}.comphet.cadd.varfilt.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream'> ${TRIO}.comphet.cadd.coding.vcf


echo -e "\n\n*********************** Comphet: Merge CADD and MetaSVM deleterious comp het variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.comphet.cadd.coding.vcf I=${TRIO}.comphet.metasvm.lof.vcf O=${TRIO}.comphet.cadd.metasvm.vcf  



echo -e "\n*********************** Filter compound heterozygous variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
####### $GATK  VariantFiltration -R $REF -V ${TRIO}.comphet.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_REC) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_REC) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.comphet.cadd.metasvm.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.comphet.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.comphet.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.comphet.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.comphet.cadd.metasvm.rare.0.vcf
uniq -f 7 ${TRIO}.comphet.cadd.metasvm.rare.0.vcf > ${TRIO}.comphet.cadd.metasvm.rare.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                    X_LINKED RECESSIVE                                        |\n"         
echo -e "\n=============================================================================================\n\n"

less ${TRIO}.merged.variants.vcf | perl -lane 'print if /#/; next if $F[0] ne "chrX";"$F[$ENV{CHILD_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2 || $1==0; "$F[$ENV{MOTHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1==$2 && $1==1; "$F[$ENV{FATHER_IDX}]" =~ /(\d+)\/(\d+):/; next if $1!=$2; print' > ${TRIO}.xlink.vcf

less ${TRIO}.xlink.vcf | perl -lane 'print if /#/; next if ("$F[$ENV{CHILD_IDX}]" =~ /^\./ || "$F[$ENV{MOTHER_IDX}]" =~ /^\./ || "$F[$ENV{FATHER_IDX}]" =~ /^\./  || $F[6] ne "PASS" || $F[5] < $ENV{MAP_QUAL_REC}); "$F[$ENV{CHILD_IDX}]" =~ /,(\d+):(\d+):(\d+)/; next if ($2 < $ENV{COV_PROB_REC}) || ($3 < $ENV{GT_QUAL_REC}); print' > ${TRIO}.xlink.filt.vcf

echo -e "\n\n*********************** Annotate x-linked recessive variants by Annovar ***********************\n\n"
perl ${ANNOVAR_DIR}/table_annovar.pl  ${TRIO}.xlink.filt.vcf  ${ANNOVAR_DIR}/humandb/  -buildver hg38 -out ${TRIO}.xlink -remove -protocol refGene,cytoBand,esp6500siv2_all,ALL.sites.2015_08,ljb26_all,exac03,dbnsfp42c,revel,intervar_20180118,cadd16all,bravo_v8,clinvar_20220320,gnomad30_genome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

echo -e "\n\n*********************** Filter x-linked variants to find MetaSVM deletrious variants *********************** \n\n"
more ${TRIO}.xlink.hg38_multianno.vcf | egrep 'RadialSVM_pred=D|ExonicFunc.refGene=stop|ExonicFunc.refGene=start|frameshift|deletion|insertion|splicing|#' | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream' > ${TRIO}.xlink.metasvm.lof.vcf


echo -e "\n\n*********************** Filter x-linked variants to find CADD deletrious variants *********************** \n\n"
######$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.hg38_multianno.vcf  -filter "vc.getAttribute('cadd16_phred').equals('.')" --filter-name "CADD" -filter "vc.getAttribute('cadd16_phred') < $CADD_THR" --filter-name "CADD" -O  ${TRIO}.xlink.cadd.hardfilt.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.hg38_multianno.vcf  -filter "vc.getAttributeAsDouble('cadd16_phred', 0) > $CADD_THR " --filter-name "High_CADD" -O  ${TRIO}.xlink.cadd.varfilt.vcf


echo -e "\n\n*********************** Exclude the non-exonic x-linked variants (intoric, UTR, synonymous, etc) *********************** \n\n"
######$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.hardfilt.varfilt.vcf -select 'vc.isNotFiltered()' -O ${TRIO}.xlink.cadd.select.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.varfilt.vcf -select 'vc.isFiltered()' -O ${TRIO}.xlink.cadd.select.vcf

less ${TRIO}.xlink.cadd.select.vcf | egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.xlink.cadd.coding.vcf


echo -e "\n\n*********************** Xlink: Merge CADD and MetaSVM deleterious x-linked variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.xlink.cadd.coding.vcf I=${TRIO}.xlink.metasvm.lof.vcf O=${TRIO}.xlink.cadd.metasvm.vcf 



echo -e "\n *********************** Filter x-linked variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
##### $GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_REC) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_REC) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.xlink.cadd.metasvm.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.xlink.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_REC) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_REC)" --filter-name "Rare" -O  ${TRIO}.xlink.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.xlink.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.xlink.cadd.metasvm.rare.0.vcf
uniq -f 7 ${TRIO}.xlink.cadd.metasvm.rare.0.vcf > ${TRIO}.xlink.cadd.metasvm.rare.vcf




echo -e "\n=============================================================================================\n"
echo -e " |                                         DOMINANT                                            |\n"         
echo -e "\n=============================================================================================\n\n"

$GATK --java-options "-Xmx4g" GenotypeGVCFs -R $REF --variant ${CHILD}.raw.g.vcf -O ${TRIO}.het.raw.vcf

cat ${TRIO}.het.raw.vcf | perl -lane 'print if /#/;next if $F[5]<30 ; $F[9]=~/(\d+)\/(\d+):/; next if $1==$2; print' >  ${TRIO}.dom.vcf


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

echo -e "\n\n*********************** Hom: Merge CADD and MetaSVM deleterious x-linked variants  *********************** \n\n"
picard-tools MergeVcfs I=${TRIO}.dom.cadd.coding.vcf I=${TRIO}.dom.metasvm.lof.vcf O=${TRIO}.dom.cadd.metasvm.vcf 

echo -e "\n *********************** Filter x-linked variants to find rare variants ***********************\n\n"
echo -e "\n\n *********************** Remove duplicate variants and  ***********************\n"
##### $GATK  VariantFiltration -R $REF -V ${TRIO}.dom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DOM) && (vc.getAttributeAsString('ExAC_ALL', null).equals('.') || vc.getAttributeAsDouble('ExAC_ALL', 0) < $MAF_DOM) && (vc.getAttributeAsString('ALL.sites.2015_08', null).equals('.') || vc.getAttributeAsDouble('ALL.sites.2015_08', 0) < $MAF_DOM) && (vc.getAttributeAsString('esp6500siv2_all', null).equals('.') || vc.getAttributeAsDouble('esp6500siv2_all', 0) < $MAF_DOM)" --filter-name "Rare" -O  ${TRIO}.dom.cadd.metasvm.varfilt.vcf
$GATK  VariantFiltration -R $REF -V ${TRIO}.dom.cadd.metasvm.vcf  -filter "(vc.getAttributeAsString('bravo_freeze8', null).equals('.') || vc.getAttributeAsDouble('bravo_freeze8', 0) < $MAF_DOM) && (vc.getAttributeAsString('AF', null).equals('.') || vc.getAttributeAsDouble('AF', 0) < $MAF_DOM)" --filter-name "Rare" -O  ${TRIO}.dom.cadd.metasvm.varfilt.vcf
$GATK SelectVariants -R $REF -V ${TRIO}.dom.cadd.metasvm.varfilt.vcf -select 'FILTER =~ Rare' -O ${TRIO}.dom.cadd.metasvm.rare.0.vcf

uniq -f 7 ${TRIO}.dom.cadd.metasvm.rare.0.vcf > ${TRIO}.dom.cadd.metasvm.rare.vcf



echo -e "\n=============================================================================================\n"
echo -e " |                                         WRAP UP                                            |\n"         
echo -e "\n=============================================================================================\n\n"
echo -e  "\n\n*********************** Move final results to a folder ***********************\n\n"

mkdir -p ${OUTPUTS_TRIO}/Final_Results
mkdir -p ${OUTPUTS_TRIO}/Final_Results/Excel_Tables

cp ${TRIO}.snp.denovo.intersect.vcf ${TRIO}.indel.denovo.intersect.vcf  ${TRIO}.hom.cadd.metasvm.rare.vcf ${TRIO}.comphet.cadd.metasvm.rare.vcf ${TRIO}.dom.cadd.metasvm.rare.vcf ${TRIO}.xlink.cadd.metasvm.rare.vcf ${OUTPUTS_TRIO}/Final_Results

echo -e "\n=============================================================================================\n"
echo -e " |                                   Generate EXCEL Reports                                    |\n"         
echo -e "\n=============================================================================================\n\n"

SNP_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.snp.denovo.intersect.vcf
INDEL_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.indel.denovo.intersect.vcf 
HOM_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.hom.cadd.metasvm.rare.vcf
COMPHET_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.comphet.cadd.metasvm.rare.vcf 
XLINK_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.xlink.cadd.metasvm.rare.vcf
DOMINANT_IN=${OUTPUTS_TRIO}/Final_Results/${SAMPLE_ID_CHILD}.dom.cadd.metasvm.rare.vcf

SNP_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.snp.denovo.intersect.csv
INDEL_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.indel.denovo.intersect.csv
HOM_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.hom.cadd.metasvm.rare.csv
COMPHET_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.comphet.cadd.metasvm.rare.csv
TEMPFILE=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.comphet.cadd.metasvm.rare.dup.noheader.csv
COMPHET_OUT2=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.comphet.cadd.metasvm.rare.dup.csv
XLINK_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.xlink.cadd.metasvm.rare.csv
DOMINANT_OUT=${OUTPUTS_TRIO}/Final_Results/Excel_Tables/${SAMPLE_ID_CHILD}.dom.cadd.metasvm.rare.csv


$GATK VariantsToTable -V $SNP_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $SNP_OUT

$GATK VariantsToTable -V $INDEL_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $INDEL_OUT

$GATK VariantsToTable -V $HOM_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $HOM_OUT

$GATK VariantsToTable -V $COMPHET_IN  --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $COMPHET_OUT

$GATK VariantsToTable -V $XLINK_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $XLINK_OUT

$GATK VariantsToTable  -V $DOMINANT_IN --show-filtered -F CHROM -F POS -F Gene.refGene -F ID -F REF -F ALT -F QUAL -F Func.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AC -F ExAC_ALL -F esp6500siv2_all -F ALL.sites.2015_08 -F bravo_freeze8 -F AF -F RadialSVM_pred -F SIFT_pred  -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_pred -F MutationTaster_pred -F REVEL -F MetaLR_pred  -F cadd16_phred -F fathmm-MKL_coding_pred -F PROVEAN_pred -F GERP++_RS_rankscore -F M-CAP_pred -F DANN_rankscore -F VEST3_score -F MutationAssessor_pred -F LRT_pred -F SiPhy_29way_logOdds_rankscore -F integrated_fitCons_rankscore -O $DOMINANT_OUT

echo -e "\n\n*********************** Remove the single compund hets in the excel sheet ***********************\n\n"

awk '{print $3}' $COMPHET_OUT | sort | uniq -d | grep -F -f - $COMPHET_OUT > $TEMPFILE
head -n1 $COMPHET_OUT > $COMPHET_OUT2
cat $TEMPFILE >> $COMPHET_OUT2
rm $TEMPFILE

