# repolarize the vcf file based on the genotypes of one outgroup sample in the vcf file
# the goal is to used the repolarized vcf file to generate the derived site frequency spectrum (DSFS) for demographic analyses
#### the code was wirtten for an ongoing project ####

repolarizeVCF = function(vcf_file="mydata.vcf", vcf.FORMAT = "GT:DP:AD:GQ:GL", outgroup_indv="sample01", out_file="repolarized.vcf", ...) {

	myvcf = read.table(vcf_file)
	
	print(paste0("input number of loci: ", nrow(myvcf)))

	cols = scan(vcf_file, skip=14, nlines=1, what=character())
	names(myvcf) = cols
	ingroup = myvcf[, -match(outgroup_indv, cols)]
	outgroup = myvcf[, c(1:9, match(outgroup_indv, cols))]

	# to repolarize: 
	## heterozygote sites in the outgroup cannot be used, as the ancestral state is unknown
	## if the outgroup genotype (GT) is 0/0: the outgroup state is the same as the current ancestral state, no change needed
	## if the outgroup genotype (GT) is 1/1: the current ancestral state and derived state should be reversed, so that the outgroup genotypes are the ancestral state

	ingroup$keep = "yes"  # remove loci later if outgroup genotypes are heterozygous
	site = NULL  ## sites where ancestral and derived states need to be reversed

	# find the sites to remove / reverse
	for (i in 1:nrow(outgroup)) {
		GT = unlist(strsplit(outgroup[i,10], split=":"))
		GT = unlist(strsplit(GT[1], split="/"))

		if (GT[1]==1 & GT[2]==1) { site = c(site,i) }      ## sites that need to be relabeled 
		if (GT[1] != GT[2]) { ingroup[i, "keep"] = "no" }  ## sites that need to be removed
	}

	
	# iterate through the rows (i.e. sites) that need to be reversed
	for (i in site) {  
		
		# relabel: 

		## exchange the current labels alternative and reference alleles
		old_ref = ingroup[i,"REF"]   
		ingroup[i,"REF"] = ingroup[i,"ALT"]
		ingroup[i,"ALT"] = old_ref

		## iterate through all individuals
		for (j in 10:(ncol(ingroup)-1)) { ## the last column would be "keep" 
		
			data = unlist(strsplit(ingroup[i, j], split=":"))

			vcf.FORMAT = unlist(strsplit(vcf.FORMAT, split=":"))
			
			n.GT = match("GT", vcf.FORMAT)
      			GT = unlist(strsplit(data[n.GT], split="/"))  # change the genotypes from 1 to 0 or 0 to 1 
      			GT[1] = as.character(1-as.numeric(GT[1]))
      			GT[2] = as.character(1-as.numeric(GT[2]))
      			data[n.GT] = paste(GT[1],GT[2],sep="/")
			
			n.AD = match("AD", vcf.FORMAT)
   		        AD = unlist(strsplit(data[n.AD], split=","))  # change the order of the allele depths
      			data[n.AD] = paste(AD[2],AD[1],sep=",")
      
      			n.GL = match("GL", vcf.FORMAT)
      			GL = unlist(strsplit(data[n.GL], split=","))  # change the order of the genotype likelihoods (only two alleles, so three genotypes)
      			data[n.GL] = paste(GL[3],GL[2],GL[1],sep=",")

			ingroup[i, j] = paste(data, collapse=":")
		}
	}

## output the repolarized data
mydata = ingroup[ingroup$keep == "yes", -ncol(ingroup)] # remove the last column "keep"

print(paste0("the number of retained loci: ", nrow(mydata)))
print(paste0("the number of retained individuals: ", (ncol(mydata)-9)))

write.table(mydata, out_file, sep="\t", row.names=F, quote=F)

print("NOTE: copy the first 14 header lines in the original vcf data into the repolarized data")

}
