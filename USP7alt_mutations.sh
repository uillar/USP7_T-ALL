#!/bin/bash

# Create a minibam with the regions of the genes of interest for faster analysis
samtools view -b $1_TARGET_STAR_hg38_aligned.bam 16:8892097-8975328 4:152320544-152536092 9:136494433-136546048 > $1_TARGET_USP7alt.bam

# Add readGroups
PicardCommandLine AddOrReplaceReadGroups \
	--INPUT $1_TARGET_USP7alt.bam \
	--OUTPUT $1_TARGET_USP7alt.RG.bam \
	--RGID 1 \
	--RGLB S9_L001 \
	--RGPL ILLUMINA \
	--RGPU BARCODE \
        --RGSM $1

samtools index $1_TARGET_USP7alt.RG.bam

# varScan
samtools mpileup -q 15 -Q 15 -A -f Homo_sapiens.GRCh38.dna.primary_assembly.fa $1_TARGET_USP7alt.RG.bam | varscan mpileup2cns - --variants --strand-filter 0 --p-value 0.5 --min-var-freq 0.05 --output-vcf 1 > $1_TARGET_USP7alt.vcf

# ANNOVAR Annotation
	
perl table_annovar.pl $1_TARGET_USP7alt.RG.bam.vcf annovar/humandb/ \
	-protocol refGene,ensGene,dbnsfp47a,cosmic70,gnomad41_genome,ALL.sites.2015_08,avsnp151,clinvar_20221231 \
	-operation g,g,f,f,f,f,f,f \
	-buildver hg38 \
	-out $1_TARGET_USP7alt.ANNOVAR \
	-remove \
	-nastring . \
	-vcfinput
 
awk 'BEGIN {FS=OFS="\t"} 
NR==1 { 
    for (i=1; i<=203; i++) printf "%s\t", $i; 
    print "GT", "GQ", "SDP", "DP", "RD", "AD", "FREQ", "PVAL", "RBQ", "ABQ", "RDF", "RDR", "ADF", "ADR"; 
    next 
} 
{ 
    split($204, a, ":"); 
    for (i=1; i<=203; i++) printf "%s\t", $i; 
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14] 
}' $1_TARGET_USP7alt.ANNOVAR.hg38_multianno.txt > $1_TARGET_USP7alt.ANNOVAR.txt

rm $1_TARGET_USP7alt.ANNOVAR.hg38_multianno.vcf
rm $1_TARGET_USP7alt.ANNOVAR.avinput
rm $1_TARGET_USP7alt.ANNOVAR.hg38_multianno.txt

awk -F'\t' '$11 == "exonic" && $9 != "synonymous SNV" && $9 !~ /nonframeshift/' $1_TARGET_USP7alt.ANNOVAR.txt > $1_TARGET_USP7alt.ANNOVAR_filtered.txt

cut -f1-18,186-217 $1_TARGET_USP7alt.ANNOVAR_filtered.txt > temp_file && mv temp_file $1_TARGET_USP7alt.ANNOVAR_filtered.txt

awk -F'\t' '($18 != "B") && ($19 == ".") && ($42 >= 15)' $1_TARGET_USP7alt.ANNOVAR_filtered.txt > $1_TARGET_USP7alt.ANNOVAR_filtered.txt
