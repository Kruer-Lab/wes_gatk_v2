

#export PERL5LIB=/usr/local/lib/x86_64-linux-gnu/perl



echo -e "\n\n*********************** Exclude the non-exonic variants (intoric, UTR, synonymous, etc) *********************** \n\n"
java -jar $GATK -T SelectVariants -nt 32 -R $REF -V ${TRIO}.xlink.cadd.hardfilt.varfilt.vcf -select 'vc.isNotFiltered()'| egrep -v 'Func.refGene=intronic|Func.refGene=UTR3|Func.refGene=UTR5|Func.refGene=intergenic|Func.refGene=ncRNA|Func.refGene=downstream|Func.refGene=upstream|ExonicFunc.refGene=synonymous'> ${TRIO}.xlink.cadd.vcf

echo -e "\n\n*********************** Comphet: Merge CADD and MetaSVM deleterious homozygous variants  *********************** \n\n"
java -jar $GATK -T CombineVariants -R $REF --variant ${TRIO}.xlink.cadd.vcf --variant ${TRIO}.xlink.metasvm.lof.vcf -o ${TRIO}.xlink.cadd.metasvm.vcf -genotypeMergeOptions UNIQUIFY

echo -e "\n\n*********************** Index the vcf file of CADD-MetaSVM deleterious compound het variants ***********************\n\n"
$IGVTOOLS index ${TRIO}.xlink.cadd.metasvm.vcf

echo -e "\n\n*********************** Filter compound heterozygous variants to find rare variants ***********************\n\n"
java -jar $GATK -T VariantFiltration -R $REF -V ${TRIO}.xlink.cadd.metasvm.vcf  --filterExpression "vc.getAttribute('ExAC_ALL') > $MAF_XLINK" --filterName "rare_ExAC" --filterExpression " vc.getAttribute('1000g2015aug_all') > $MAF_XLINK" --filterName "rare_1000G" --filterExpression "vc.getAttribute('esp6500siv2_all') > $MAF_XLINK" --filterName "rare_esp" -o ${TRIO}.xlink.cadd.metasvm.rare.varfilt.vcf

java -jar $GATK -T SelectVariants -nt 32 -R $REF -V ${TRIO}.xlink.cadd.metasvm.rare.varfilt.vcf -select 'vc.isNotFiltered()'  -o ${TRIO}.xlink.cadd.metasvm.rare.vcf













#slivar expr --vcf ${TRIO}.merged.variants.vcf --ped ${PEDIGREES}/pedigree_${CHILD_SAMPLE_NAME}_slivar.ped -o ${TRIO}.slivar.expr.bcftools csq -s - --ncsq 40 -g $GFF -l -f $REF ${TRIO}.slivar.expr.vcf -O u  > ${TRIO}.test.vcf
#| slivar compound-hets --sample-field denovo -p ${PEDIGREES}/pedigree_${CHILD_SAMPLE_NAME}_slivar.ped > ${TRIO}.slivar.comphet.vcf

#echo $GFF

#bcftools csq --phase s  -g $GFF -f $REF ${TRIO}.slivar.expr.vcf  -O u | slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped > $out
echo "pedigreeee: $ped"

#bcftools csq --phase s  -g $GFF -f $REF ${TRIO}.slivar.expr.vcf   | slivar compound-hets --sample-field comphet_side --sample-field denovo -p ${PEDIGREES}/pedigree_${CHILD_SAMPLE_NAME}_slivar.ped > ${TRIO}.slivar.comphet.vcf



#bcftools csq -g $GFF -f $REF ${TRIO}.slivar.expr.vcf > ${TRIO}.test2.vcf







#bcftools csq -s - --ncsq 40 -l  -g $GFF -f $REF ${TRIO}.slivar.expr.vcf > ${TRIO}.test2.vcf
#slivar compound-hets -a  -v  ${TRIO}.test2.vcf --skip NONE --sample-field comphet_side --sample-field denovo -p ${PEDIGREES}/pedigree_${CHILD_SAMPLE_NAME}_slivar.ped -o ${TRIO}.slivar.comphet.vcf

#slivar compound-hets -v  ${TRIO}.comphet.hg38_multianno.vcf  --sample-field comphet_side  --sample-field denovo -p ${PEDIGREES}/pedigree_${CHILD_SAMPLE_NAME}_slivar.ped > ${TRIO}.slivar.comphet.vcf



 

#echo $VEP
#$VEP -i ${TRIO}.slivar.expr.vcf --cache -o ${TRIO}.slivar.expr.ann.vcf


#$VEP -i ${TRIO}.slivar.expr.vcf --gff $GFF --fasta $REF
