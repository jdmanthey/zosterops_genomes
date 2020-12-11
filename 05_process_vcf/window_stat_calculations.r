#functions used to calculate diversity statistics within populations
# make sure the output file is already written in format (header):
# population1   population2     statistic       chromosome      start   end     number_sites    number_variable_sites   calculated_stat
# e.g., write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=outname, sep="\t")

# input is simplified vcf (entire vcf) and 3-column popmap (ordered the same as the vcf) and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only calculate for invariant sites and biallelic snps
heterozygosity <- function(xxx, popmap, outname, filename) {
        options(scipen=999)
        xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
        for(a in 1:nrow(popmap)) {
                output_rep <- c(popmap[a,1], "none", "heterozygosity", strsplit(filename, ":")[[1]][1], 
                        as.numeric(strsplit(strsplit(filename, ":")[[1]][2], "-")[[1]][1]),
                        as.numeric(strsplit(strsplit(filename, "-")[[1]][2], ".simple")[[1]][1]))
                # select this individual
                a_rep <- xxx[,a+3]
                # remove missing genotypes
                a_rep <- a_rep[a_rep != "./."]
                # count number of sites
                a_total <- length(a_rep)
                # remove phasing information
                a_rep <- gsub("\\|", "/", a_rep)
                # find number of heterozygous sites
                a_het <- length(a_rep[a_rep == "0/1"])
                # add to output
                output_rep <- c(output_rep, a_total, a_het, a_het/a_total)
                write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
        }
}

# helper function to identify transition variants from biallelic alleles in rows (used in titv function)
transitions <- function(x1) {
	a_rep <- sort(x1)
	if(a_rep[1] == "A" & a_rep[2] == "G") {
		return(1)
	} else if(a_rep[1] == "C" & a_rep[2] == "T") {
		return(1)
	} else {
		return(0)
	}
}



# input is simplified vcf (entire vcf) and other files as in heterozygosity calcs, to calculate transition transversion ratio
titv <- function(xxx, popmap, outname, filename) {
        options(scipen=999)
        xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
        # filter to variable sites
        xxx <- xxx[xxx[,3] != ".", ]
        
        # sites used for titv 
        # total sites used same as het sites
        a_total <- nrow(xxx)
        a_het <- a_total
        
        # find number of transitions using helper function
        total_transitions <- sum(apply(xxx[,2:3], 1, transitions))
        # calculate transition / transversion ratio
        titv_ratio <- total_transitions / (a_total - total_transitions)
        
        # output info
        output_rep <- c("all_inds", "none", "titv", strsplit(filename, ":")[[1]][1], 
                        as.numeric(strsplit(strsplit(filename, ":")[[1]][2], "-")[[1]][1]),
                        as.numeric(strsplit(strsplit(filename, "-")[[1]][2], ".simple")[[1]][1]))
        
        output_rep <- output_rep <- c(output_rep, a_total, a_het, titv_ratio)
        write(output_rep, file=outname, append=T, ncolumns=9, sep="\t")
        
}




# input is one simplified vcfs,  
# the popmap, and output file name
# and the input simple vcf file name that contains the chr and start and end information
# only for invariant sites and biallelic snps
# lastly, input number of sites needed to write file

create_fasta_from_vcf <- function(xxx, popmap, outname, filename, num_sites_needed) {
        options(scipen=999)
        # remove sites that are not either invariant or bi-allelic SNPs
        xxx <- xxx[nchar(xxx[,2]) == 1 & nchar(xxx[,3]) == 1, ]
        
        # keep going only if enough sites
        if(num_sites_needed <= nrow(xxx)) {
                # remove phasing information
                for(a in 4:ncol(xxx)) {
                        xxx[,a] <- gsub("\\|", "/", xxx[,a])
                }
                
                # define names of individuals in output fasta
                output_names_fasta <- paste(">", popmap[,1], sep="")
                
                # subset the genotypes from the allele info
                allele_info <- xxx[,2:3]
                genotypes <- xxx[,4:ncol(xxx)]
                
                # convert all numbers in genotypes to actual bases and ambiguities
                for(a in 1:nrow(genotypes)) {
                        if(allele_info[a,2] == ".") { # if non-polymorphic
                                genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
                                genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
                        } else { # if polymorphic
                                both_genotypes <- sort(as.character(allele_info[a,]))
                                if(both_genotypes[1] == "A" & both_genotypes[2] == "C") { het = "M" }
                                if(both_genotypes[1] == "A" & both_genotypes[2] == "G") { het = "R" }
                                if(both_genotypes[1] == "A" & both_genotypes[2] == "T") { het = "W" }
                                if(both_genotypes[1] == "C" & both_genotypes[2] == "G") { het = "S" }
                                if(both_genotypes[1] == "C" & both_genotypes[2] == "T") { het = "Y" }
                                if(both_genotypes[1] == "G" & both_genotypes[2] == "T") { het = "K" }
                                genotypes[a,] <- gsub("0/0", allele_info[a,1], genotypes[a,])
                                genotypes[a,] <- gsub("\\./\\.", "?", genotypes[a,])
                                genotypes[a,] <- gsub("1/1", allele_info[a,2], genotypes[a,])
                                genotypes[a,] <- gsub("0/1", het, genotypes[a,])
                        }
                }
  
                  # write output
                for(a in 1:ncol(genotypes)) {
                        if(a == 1) {
                                write(output_names_fasta[a], file=outname, ncolumns=1)
                                write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
                        } else {
                                write(output_names_fasta[a], file=outname, ncolumns=1, append=T)
                                write(paste(genotypes[,a], collapse=""), file=outname, ncolumns=1, append=T)
                        }
                }
        }
}
