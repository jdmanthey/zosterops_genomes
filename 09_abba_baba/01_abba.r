# make sure output directory written already 
# in hpcc : mkdir output

abba_kulamb_murphyi <- function(xxx, output_number) {
	# read in vcf window
	x <- read.table(xxx, sep="\t", stringsAsFactors=F)
	
	# some basic info about the structure of the input vcf
	# order of individuals in vcf: 
	# 1	10	11	12	13	14	15	16	17	2	3	4	5	6	7	8	9

	# set up names of individuals
	ind_names <- c("Z_griseotinctus_2003067", "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b",
	"Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819",
	"Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056",
	"Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")
	
	# set up names of taxa for the tests
	abba_inds <- list()
	#Z_griseotinctus_2003067, Z_murphyi_CEF838, Z_kulambangrae_VGR592
	abba_inds[[1]] <- c(1, 2, 12)
	# Z_splendidus_CEF825b, Z_kulambangrae_VGR592, Z_murphyi_CEF838
	abba_inds[[2]] <- c(4, 12, 2)
	# Z_griseotinctus_2003067, Z_murphyi_CEF838, Z_kulambangrae_CEF484
	abba_inds[[3]] <- c(1, 2, 11)
	#Z_splendidus_CEF825b, Z_kulambangrae_CEF484, Z_murphyi_CEF838
	abba_inds[[4]] <- c(4, 11, 2)
	
	# write initial output file
	write(paste("ind1", "ind2", "ind3", "n_snps", "d", "z-score", sep="\t"), file=paste("output/abba_", output_number, ".txt", sep=""), ncolumns=6)
	
	# subsets for each group of individuals and remove sites with missing data
	# and keep only sites biallelic in the quartet of individuals
	abba_data_list <- list()
	for(b in 1:length(abba_inds)) {
		b_columns <- c(1,2,3,abba_inds[[b]]+3)
		abba_data_list[[b]] <- x[, b_columns]
		# remove rows with missing data
		abba_data_list[[b]] <- abba_data_list[[b]][(grepl("\\./\\.", abba_data_list[[b]][,4])) == FALSE,]
		abba_data_list[[b]] <- abba_data_list[[b]][(grepl("\\./\\.", abba_data_list[[b]][,5])) == FALSE,]
		abba_data_list[[b]] <- abba_data_list[[b]][(grepl("\\./\\.", abba_data_list[[b]][,6])) == FALSE,]
		# keep biallelic sites
		abba_data_list[[b]] <- abba_data_list[[b]][apply(abba_data_list[[b]], 1, find_biallelic), ]
	}	
	
	# find d statistic and run randomizations for each set of individuals
	for(b in 1:length(abba_data_list)) {
		# subset the data to keep only the genotypes for this rep
		b_rep <- abba_data_list[[b]][,4:6]
		# replace genotypes with derived allele frequencies
		b_rep[b_rep == "0/0"] <- 0
		b_rep[b_rep == "0/1"] <- 0.5
		b_rep[b_rep == "1/1"] <- 1
		# make sure they are numeric
		b_rep <- data.frame(col1 = as.numeric(b_rep[,1]), col2 = as.numeric(b_rep[,2]), col3 = as.numeric(b_rep[,3]))
		
		# find abba frequencies
		abba_rep <- apply(b_rep, 1, find_abba)
		# find baba frequencies
		baba_rep <- apply(b_rep, 1, find_baba)
		# find d statistic
		d_rep <- sum(abba_rep - baba_rep) / sum(abba_rep + baba_rep)
		
		# 100 randomizations to get a z-score
		d_random_total <- c()
		for(a in 1:100) {
			a_rep <- abba_data_list[[b]][,4:6]
			a_rep <- t(apply(a_rep, 1, rand_geno))
			# replace genotypes with derived allele frequencies
			a_rep[a_rep == "0/0"] <- 0
			a_rep[a_rep == "0/1"] <- 0.5
			a_rep[a_rep == "1/1"] <- 1
			# make sure they are numeric
			a_rep <- data.frame(col1 = as.numeric(a_rep[,1]), col2 = as.numeric(a_rep[,2]), col3 = as.numeric(a_rep[,3]))
		
			# find abba frequencies
			abba_random <- apply(a_rep, 1, find_abba)
			# find baba frequencies
			baba_random <- apply(a_rep, 1, find_baba)
			# find d statistic
			d_random <- sum(abba_random - baba_random) / sum(abba_random + baba_random)
			d_random_total <- c(d_random_total, d_random)
		}
		
		# find z-score
		z_rep <- (d_rep - 0) / sd(d_random_total)
		
		output_rep <- c(ind_names[abba_inds[[b]][1]], ind_names[abba_inds[[b]][2]], ind_names[abba_inds[[b]][3]],
			nrow(b_rep), d_rep, z_rep)
		write(output_rep, file=paste("output/abba_", output_number, ".txt", sep=""), sep="\t", ncolumns=6, append=T)		
	}
	

}

# find abba statistic
# outgroup genotype = "0/0" for each allele, so frequency of derived allele = 0 always
find_abba <- function(xxxx) {
	return((1 - xxxx[1]) * xxxx[2] * xxxx[3] * (1 - 0))
}
# find baba statistic
# outgroup genotype = "0/0" for each allele, so frequency of derived allele = 0 always
find_baba <- function(xxxx) {
	return(xxxx[1] * (1 - xxxx[2]) * xxxx[3] * (1 - 0))	
}

# randomize genotype order
rand_geno <- function(xxxx) {
	x_samp <- xxxx
	xxxx <- sample(x_samp, 3, replace=T)
}




# find bi-allelic function
# use with apply over rows with 3 individuals (6 columns)
find_biallelic <- function(xxxx) {
	return(length(unique(xxxx[4:6])) > 1)
}





