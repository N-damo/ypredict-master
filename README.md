# Ypredict   
**Ypredict** is a python based software package that predicts y chromosome haplogroup. Here, I use calculate rank method to automatically find the most likely y haplogroup. For each y haplogroup, I give two mark (T or F) according their snp calling state. For Example, if the haplogroup O2a1a1a2a1 in isogg (<https://isogg.org/tree/>) haplogroup tree has six snps, then I calculate the number of true happening snp in the tested sample. For each haplogroup, I add all of the  rank along the routine from the 'Y' haplogroup to this haplogroup(rank = n_T + ((n_T**2)/(n_T + n_F))). After that, the max rank and 'T' mark of the haplogroup will be the most likely haplogroup. 
* The current version is 0.0.1  
  
# Dependence  
* vcftools(<http://vcftools.sourceforge.net/>)
* biopython(<https://biopython.org/wiki/Download>)

# Getting Started
***
## Step1  
Download y haplogroup tree from isogg. Then, filter snp by snpfilter.py. In this step, hotspot and backmutate snp will be removed. Finally, two files map.json and ref_vcf.gz will be generated.  
`python snpfilter.py -snp snp14.3.csv`  
## Step2  
In this step, we will use the file ref_vcf.gz generated by the step1 to make snp calling using gatk3.8 UnifiedGenotyper module. Critically, we use hg38 refrence genome in this step.

`java -Xmx32g -jar GenomeAnalysisTK.jar -T UnifiedGenotyper 
-R hg38.fa -I *.bam -o y.vcf.gz 
--intervals chrY 
-ploidy 1 
--output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES 
--alleles ref_vcf.gz `  
## Step3  
Y chromosome haplogroup can be predicted by ypredict.py. In the this step, the scrip will automatically output the most likely haplogroup. The final result can be seen in ypredict.txt. More detail output write in ystatistics.csv.

`python ypredict.py -vcf y.vcf.gz -s hfspecial.xlsx -m map.json`

If you need to update y haplogroup tree file downloaded from isogg, you can redo step1 and get an updated map.json.  
# Install  
`git clone https://github.com/N-damo/ypredict-master.git`