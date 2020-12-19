# repolarize-vcf
This code repolarizes a vcf file using one sample's genotypes as the ancestral state. 
#### the code was wirtten for an ongoing project ####

In many cases, the outgroup reference genome is not available, and both the ingroup and outgroup samples are genotyped by mapping to the ingroup reference genome. This makes the vcf file not polarized correctly, because the ancestral allele is the allele of the ingroup's reference genome, rather than that of the outgroup. However, when the outgroup sample is available, their genotypes can also indicate the ancestral condition, and thus the vcf file can be repolarized so that the information of ancestral/derived alleles can be used to calculate the derived site frequency spectrum (DSFS) for some analyses such as demography. 

To repolarize the vcf file, the heterozygous sites in the outgroup cannot be used because the ancestral state would be unknown. If the outgroup genotype (GT) is 0/0, then the outgroup state is the same as the current ancestral state and no change is needed. If the outgroup genotype (GT) is 1/1, then the current ancestral state and derived state should be reversed, so that the outgroup genotypes are the ancestral state.

## input file
This code was created for the vcf file of diploid individuals and biallelic SNP loci. The input file should have no missing data. There should be only one outgroup sample and this outgroup sample should be the last one in the vcf file (i.e. the last column). 
# NOTE: The vcf file needs to be modified before running the codes. To include the header lines, delete the # in the beginning of the vcf header line (i.e. in line 15, change #CHROM to CHROM). This can be done using the Notepad ++ or the unix nano command etc. Save the modified file and name it as myvcf.vcf. Put this file and the R script in the same directory to run the script.  

## output file
The script will output a file called mydata.vcf. All ingroup smaples are output in the original order, and the outgroup sample will be excluded. The number of loci kept in this file will be printed on the screen. Note that this file is not in the vcf format because the first 14 lines are not included. To make it into the vcf format, first add the # back in front of the header line (e.g. change CHROM to #CHROM). Second, copy the first 14 lines (all start with ##) in the original vcf data into the repolarized data. Again this can easily be done in Notepad or unix. 
