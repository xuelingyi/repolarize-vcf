# repolarize the vcf file based on the genotypes of one outgroup sample in the vcf file
# the goal is to used the repolarized vcf file to generate the derived site frequency spectrum (DSFS) for demographic analyses
#### the code was wirtten for an ongoing project ####

# NOTE: to include the header lines, delete the # in the beginning of the vcf header line (i.e. in line 15, change #CHROM to CHROM), save this file and name it as myvcf.vcf

# run in R

myvcf = read.table("myvcf.vcf", header=T)  # a vcf file of diploid individuals without missing data; all samples are genotyped by mapping to the reference genome of the ingroup species

ingroup = myvcf[,1:(ncol(myvcf)-1)]  ## there is only one outgroup, and the outgroup is the last sample in the data
outgroup = myvcf[,c(1:9,ncol(myvcf))]

# to repolarize: 
## heterozygote sites in the outgroup cannot be used, as the ancestral state is unknown
## if the outgroup genotype (GT) is 0/0: the outgroup state is the same as the current ancestral state, no change needed
## if the outgroup genotype (GT) is 1/1: the current ancestral state and derived state should be reversed, so that the outgroup genotypes are the ancestral state

ingroup$keep = "yes"  # remove loci later if outgroup genotypes are heterozygous
site = NULL  ## sites where ancestral and derived states need to be reversed

## find the sites
for (i in 1:nrow(outgroup)) {
GT = unlist(strsplit(outgroup[i,10], split=":"))
GT = unlist(strsplit(GT[1], split="/"))

if (GT[1]==1 & GT[2]==1) { 
site = c(site,i) } 

if (GT[1] != GT[2]) { 
ingroup[i, "keep"] = "no" }  ## sites that need to be removed
}

# iterate through the rows (i.e. sites) that need to be reversed
for (i in site) {
## exchange the current alternative and reference alleles

old_ref = ingroup[i,"REF"]   
ingroup[i,"REF"] = ingroup[i,"ALT"]
ingroup[i,"ALT"] = old_ref

## iterate through all individuals
for (j in 10:(ncol(myvcf)-1)) { 
data = unlist(strsplit(ingroup[i, j], split=":"))

GT = unlist(strsplit(data[1], split="/"))  # change the genotypes from 1 to 0 or 0 to 1 
GT[1] = as.character(1-as.numeric(GT[1]))
GT[2] = as.character(1-as.numeric(GT[2]))
data[1] = paste(GT[1],GT[2],sep="/")

AD = unlist(strsplit(data[3], split=","))  # change the order of the allele depths
data[3] = paste(AD[2],AD[1],sep=",")

GL = unlist(strsplit(data[5], split=","))  # change the order of the genotype likelihoods (only two alleles, so three genotypes)
data[5] = paste(GL[3],GL[2],GL[1],sep=",")

ingroup[i, j] = paste(data, collapse=":")
}
}

## output the repolarized data
mydata = ingroup[ingroup$keep == "yes", 1:(ncol(myvcf)-1)]
nrow(mydata) # the number of loci kept in the data 
write.table(mydata, "repolarized-data.vcf", row.names=F, quote=F, sep="\t")

## NOTE: the header lines are not included in the output file; to make the real vcf file:
### add the # in front of the header line
### copy the first few lines (all starting with ##) in the original vcf data into the repolarized data
