# repolarize-vcf
This script repolarizes a vcf file using one outgroup sample's genotypes as the ancestral state. 

The script was wirtten for the project published as: Yi, X., & Latch, E. K. (2022) Nuclear phylogeography reveals strong impacts of gene flow in big brown bats. Journal of Biogeography, 49(6), 1061– 1074. https://doi.org/10.1111/jbi.14362

When the ingroup and outgroup samples are genotyped by mapping to the same reference genome, such as the one of the ingroup species, the identified genotypes may not be polarized correctly in the vcf file, because the ancestral allele is the allele of the ingroup's reference genome, rather than that of the outgroup. However, when the outgroup sample is available, their genotypes can also indicate the ancestral condition, and thus the vcf file can be repolarized so that the information of ancestral/derived alleles can be used to calculate the derived site frequency spectrum (DSFS) for some analyses such as demography. 

To repolarize the vcf file, the heterozygous sites in the outgroup cannot be used because the ancestral state would be unknown. If the outgroup genotype (GT) is 0/0, then the outgroup state is the same as the current ancestral state and no change is needed. If the outgroup genotype (GT) is 1/1, then the current ancestral state and derived state should be reversed, so that the outgroup genotypes are the ancestral state.

## input file
This code was created for the vcf file of diploid individuals and biallelic SNP loci. The input file has no missing data and only one outgroup sample. 

## run the script
The script is written as a function in R. Once codes are downloaded, go to R and source the script "vcf-repo.R". Then run the function repolarizeVCF(vcf_file="mydata.vcf", outgroup_indv="sample01", out_file="repolarized.vcf").

## output file
The script will output a vcf file of only the ingroup individuals and the retained loci. The number of loci and individuals kept in this file will be printed on the screen. Note that this file does not have the vcf header lines. To add the header lines, go to the original input file and copy the first 14 lines then paste them into the repolarized file (this can be easily done using text editors such as Notepad or nano). 
