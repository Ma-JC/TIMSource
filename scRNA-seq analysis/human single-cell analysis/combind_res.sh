#!/bin/bash
>all_sample.txt
for i in ./sample/*.hg19_multianno.xls; do
	sample=$(echo $i | awk -F '.' '{print $2}')
	echo $sample  
	cut_output=$(cut -f 1,2,3,4,5,18,14,15,21,17,13 $i)
	sed_output=$(echo "$cut_output" | sed '1d')
	final_output=$(echo "$sed_output" | awk -v sample="$sample" '{print $0 "\t" sample}')
	echo "$final_output" >> all_sample.txt
done
sed -i '1i\Chr\tStart\tEnd\tRef\tAlt\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.ensGene\tHugo_Symbol\tTumor_Sample_Barcode' all_sample.txt

