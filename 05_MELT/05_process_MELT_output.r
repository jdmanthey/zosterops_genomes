# output base name
out_name <- "MELT_zosterops"
# set maximum allowable divergence from consensus
max_div <- 0.15
# set minimum distance between insertions
min_dist <- 100
# set maximum missing data for a call
max_missing <- 0.25

# list melt vcf files and setup names and input headers
x_files <- list.files(pattern="*comp.vcf")
x_names <- substr(x_files, 1, nchar(x_files) - 15)
x_headers <- c("scaffold", "position", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b", "Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819", "Z_griseotinctus_2003067", "Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056", "Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")
sample_names <- c("Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b", "Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819", "Z_griseotinctus_2003067", "Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056", "Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")
info_headers <- c("target_site_dup", "assess_qual", "near_gene", "SV_type", "SV_length", "ME_info", "cov_diff_to_reference", "discordant_pairs_left_side", "discordant_pairs_right_side", "left_to_right_discordant_ratio", "based_on_prior_info", "number_SRs")

# read in consensus sequences to get lengths
x_fasta <- list.files(pattern="*fasta")
consensus_lengths <- c()
for(a in 1:length(x_fasta)) {
	a_rep <- scan(x_fasta[a], what="character")
	consensus_lengths <- c(consensus_lengths, nchar(a_rep[2]))
}

# read in vcf files and initial processing
vcf <- list()
vcf_summary_filter <- x_names
# read in
for(a in 1:length(x_files)) {
	vcf[[a]] <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
}
# number calls initially pre-filtering
vcf_summary_filter <- cbind(vcf_summary_filter, unlist(lapply(vcf, nrow)))

# remove reads not passing filter
for(a in 1:length(vcf)) {
	vcf[[a]] <- vcf[[a]][vcf[[a]][,7] == "PASS", ]
}
vcf_summary_filter <- cbind(vcf_summary_filter, unlist(lapply(vcf, nrow)))

# remove reads with ASSESS < 3 (i.e., not 3, 4, or 5)
for(a in 1:length(vcf)) {
	vcf[[a]] <- vcf[[a]][as.numeric(sapply(strsplit(sapply(strsplit(vcf[[a]][,8], ";"), "[[", 2), "="), "[[", 2)) > 2,]
}
vcf_summary_filter <- cbind(vcf_summary_filter, unlist(lapply(vcf, nrow)))

# remove reads with % divergence greater than expected and those with less than 10% sequence coverage
for(a in 1:length(vcf)) {
	a_length_covered <- as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(vcf[[a]][,8], ";"), "[[", 7), ":"), "[[", 1), "="), "[[", 2))
	a_mutations <- strsplit(sapply(strsplit(sapply(strsplit(vcf[[a]][,8], ";"), "[[", 7), ":"), "[[", 2), ",")
	b_mutations <- c()
	for(b in 1:length(a_mutations)) {
		b_rep <- a_mutations[[b]]
		b_rep <- length(b_rep) - length(grep("-", b_rep))
		b_mutations <- c(b_mutations, b_rep)
	}
	vcf[[a]] <- vcf[[a]][a_length_covered >= 0.1 & (b_mutations / consensus_lengths[a]) <= max_div, ]
}
vcf_summary_filter <- cbind(vcf_summary_filter, unlist(lapply(vcf, nrow)))

# concatenate all vcf files
for(a in 1:length(vcf)) {
	if(a == 1) {
		total_vcf <- cbind(rep(x_names[a], nrow(vcf[[a]])), vcf[[a]])
	} else {
		total_vcf <- rbind(total_vcf, cbind(rep(x_names[a], nrow(vcf[[a]])), vcf[[a]]))
	}
}

# remove calls close to one another (use consensus with greater coverage)
total_vcf_copy <- total_vcf
keep <- c()
for(a in 1:nrow(total_vcf)) {
	a_rep <- total_vcf[a,]
	test_rep <- total_vcf[-a, ]
	test_rep <- test_rep[test_rep[,2] == a_rep[1,2], ]
	if(nrow(test_rep) > 0) {
		test_rep2 <- test_rep[abs(test_rep[,3] - a_rep[1,3]) <= min_dist, ]
		if(nrow(test_rep2) > 0) {
			# choose element with greatest coverage of consensus
			a_length_covered <- as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(a_rep[,9], ";"), "[[", 7), ":"), "[[", 1), "="), "[[", 2))
			b_length_covered <- as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(test_rep2[,9], ";"), "[[", 7), ":"), "[[", 1), "="), "[[", 2))
			if(length(length(a_length_covered > b_length_covered)[(a_length_covered > b_length_covered) == FALSE]) == 0)  {
				keep <- c(keep, a)
			}
		} else {
			keep <- c(keep, a)	
		}
	} else {
		keep <- c(keep, a)
	}
}
total_vcf <- total_vcf[keep,]
# summarize previous step
no_dups <- c()
for(a in 1:length(x_names)) {
	a_rep <- nrow(total_vcf[total_vcf[,1] == x_names[a], ])
	no_dups <- c(no_dups, a_rep)
}
vcf_summary_filter <- cbind(vcf_summary_filter, no_dups)


# remove calls with too much missing data
keep <- c()
for(a in 1:nrow(total_vcf)) {
	a_rep <- as.character(as.vector(total_vcf[a, 11:ncol(total_vcf)]))
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 1)
	a_rep2 <- length(a_rep[a_rep == "./."]) / length(a_rep)
	if(a_rep2 <= max_missing) {
		keep <- c(keep, a)
	}
}
total_vcf <- total_vcf[keep,]
# summarize previous step
not_missing <- c()
for(a in 1:length(x_names)) {
	a_rep <- nrow(total_vcf[total_vcf[,1] == x_names[a], ])
	not_missing <- c(not_missing, a_rep)
}
vcf_summary_filter <- cbind(vcf_summary_filter, not_missing)

# filter based on SR total > 2
total_vcf <- total_vcf[as.numeric(sapply(strsplit(sapply(strsplit(total_vcf[,9], ";"), "[[", 12), "="), "[[", 2)) > 2, ]
# summarize previous step
SR <- c()
for(a in 1:length(x_names)) {
	a_rep <- nrow(total_vcf[total_vcf[,1] == x_names[a], ])
	SR <- c(SR, a_rep)
}
vcf_summary_filter <- cbind(vcf_summary_filter, SR)

# remove calls with no variability
keep <- c()
for(a in 1:nrow(total_vcf)) {
	a_rep <- as.character(as.vector(total_vcf[a, 11:ncol(total_vcf)]))
	a_rep <- sapply(strsplit(a_rep, ":"), "[[", 1)
	a_rep2 <- length(a_rep[a_rep == "0/0"]) / length(a_rep)
	if(a_rep2 < 1) {
		keep <- c(keep, a)
	}
}
total_vcf <- total_vcf[keep,]
# summarize previous step
non_variable <- c()
for(a in 1:length(x_names)) {
	a_rep <- nrow(total_vcf[total_vcf[,1] == x_names[a], ])
	non_variable <- c(non_variable, a_rep)
}
vcf_summary_filter <- cbind(vcf_summary_filter, non_variable)

# filter for SV length >= 100
total_vcf <- total_vcf[as.numeric(sapply(strsplit(sapply(strsplit(total_vcf[,9], ";"), "[[", 5), "="), "[[", 2)) >= 100, ]
# summarize previous step
SV_length <- c()
for(a in 1:length(x_names)) {
	a_rep <- nrow(total_vcf[total_vcf[,1] == x_names[a], ])
	SV_length <- c(SV_length, a_rep)
}
vcf_summary_filter <- cbind(vcf_summary_filter, SV_length)




# write outputs

# summary
summary_output <- data.frame(Name=as.character(vcf_summary_filter[,1]), Initial=as.numeric(vcf_summary_filter[,2]), Pass_Filter=as.numeric(vcf_summary_filter[,3]), Assess=as.numeric(vcf_summary_filter[,4]), Low_Cons_Cov=as.numeric(vcf_summary_filter[,5]), Non_duplicate=as.numeric(vcf_summary_filter[,6]), Missing_filter=as.numeric(vcf_summary_filter[,7]), SR_filter=as.numeric(vcf_summary_filter[,8]), Variable=as.numeric(vcf_summary_filter[,9]), SV_length_filter=as.numeric(vcf_summary_filter[,10]))
write.table(summary_output, file="zosterops_filter_summary.txt", sep="\t", quote=F, row.names=F)

# genotypes after filtering
# colnames = name, chromosome, location, tsd, sv_length, perc_divergence, raw_info, "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b", "Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819", "Z_griseotinctus_2003067", "Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056", "Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869"

# find % divergence from consensus
perc_divergence <- c()
for(a in 1:nrow(total_vcf)) {
	a_length_covered <- as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(total_vcf[a,9], ";"), "[[", 7), ":"), "[[", 1), "="), "[[", 2))
	a_mutations <- strsplit(sapply(strsplit(sapply(strsplit(total_vcf[a,9], ";"), "[[", 7), ":"), "[[", 2), ",")[[1]]
	b_mutations <- length(a_mutations) - length(grep("-", a_mutations))
	a_consensus <- consensus_lengths[match(total_vcf[a,1], x_names)]
	perc_divergence <- c(perc_divergence, b_mutations / a_consensus)
}

output <- data.frame(Name=as.character(total_vcf[,1]),
					 Chromosome=as.character(total_vcf[,2]),
					 Location=as.numeric(total_vcf[,3]),
					 TSD=as.character(sapply(strsplit(sapply(strsplit(total_vcf[,9], ";"), "[[", 1), "="), "[[", 2)),
					 SV_length=as.numeric(sapply(strsplit(sapply(strsplit(total_vcf[,9], ";"), "[[", 5), "="), "[[", 2)),
					 Perc_Divergence=as.numeric(perc_divergence),
					 raw_info=as.character(total_vcf[,9]), 
					 Z_murphyi_CEF838=as.character(sapply(strsplit(total_vcf[,11], ":"), "[[", 1)),
					 Z_rendovae_33262=as.character(sapply(strsplit(total_vcf[,12], ":"), "[[", 1)),
					 Z_splendidus_CEF825b=as.character(sapply(strsplit(total_vcf[,13], ":"), "[[", 1)),
					 Z_stresemanni_32763=as.character(sapply(strsplit(total_vcf[,14], ":"), "[[", 1)),
					 Z_teteparius_CES726=as.character(sapply(strsplit(total_vcf[,15], ":"), "[[", 1)),
					 Z_ugiensis_32795=as.character(sapply(strsplit(total_vcf[,16], ":"), "[[", 1)),
					 Z_ugiensis_35017=as.character(sapply(strsplit(total_vcf[,17], ":"), "[[", 1)),
					 Z_vellalavella_CEF819=as.character(sapply(strsplit(total_vcf[,18], ":"), "[[", 1)),
					 Z_griseotinctus_2003067=as.character(sapply(strsplit(total_vcf[,19], ":"), "[[", 1)),
					 Z_japonicus_11220=as.character(sapply(strsplit(total_vcf[,20], ":"), "[[", 1)),
					 Z_kulambangrae_CEF484=as.character(sapply(strsplit(total_vcf[,21], ":"), "[[", 1)),
					 Z_kulambangrae_VGR592=as.character(sapply(strsplit(total_vcf[,22], ":"), "[[", 1)),
					 Z_luteirostris_PRS2615=as.character(sapply(strsplit(total_vcf[,23], ":"), "[[", 1)),
					 Z_metcalfi_32056=as.character(sapply(strsplit(total_vcf[,24], ":"), "[[", 1)),
					 Z_metcalfi_32096=as.character(sapply(strsplit(total_vcf[,25], ":"), "[[", 1)),
					 Z_metcalfi_34483=as.character(sapply(strsplit(total_vcf[,26], ":"), "[[", 1)),
					 Z_metcalfi_34869=as.character(sapply(strsplit(total_vcf[,27], ":"), "[[", 1))
					 )
save.image("zosterops_melt.RData")


# DONE FILTERING
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################







