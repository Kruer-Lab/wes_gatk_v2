#!/bin/bash
source settings_snv.sh

TRIO=$1
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
