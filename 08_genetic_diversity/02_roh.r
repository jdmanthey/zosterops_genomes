# Runs of homozygosity

# set ROH cutoff (minimum number of SNPs to call an ROH)
ROH_cutoff <- 100

# order of individuals in vcf: 
# 1	10	11	12	13	14	15	16	17	2	3	4	5	6	7	8	9

# set up names of individuals
ind_names <- c("Z_griseotinctus_2003067", "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b",
	"Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819",
	"Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056",
	"Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")

# list vcf files
x_files <- list.files(pattern="*vcf")

# write initial output file
# output file for ROH
# columns = individual, scaffold, start, end, number of snps in ROH
write(paste("#Calculated with an ROH minimum size cutoff of ", ROH_cutoff, " SNPs", sep=""), file="runs_of_homozygosity.txt", sep="\t", ncolumns=1)
write(c("individual", "scaffold", "start", "end", "n_snps"), file="runs_of_homozygosity.txt", sep="\t", ncolumns=5, append=T)


# loop for each vcf file	
for(a in 1:length(x_files)) {
	a_rep <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
	a_length <- a_rep[nrow(a_rep), 1] - a_rep[1, 1] + 1
	
	# loop for each individual
	for(b in 1:length(ind_names)) {
		b_rep <- a_rep[,c(1, b+3)]
		# remove non genotyped sites
		b_rep <- b_rep[b_rep[,2] != "./.", ]
		
		# identify rows that are heterozygous sites
		het_sites1 <- seq(from=1, to=nrow(b_rep), by=1)
		het_sites2 <- het_sites1[b_rep[,2] == "0/1" | b_rep[,2] == "0|1" | b_rep[,2] == "0/2" | b_rep[,2] == "0|2" |
								b_rep[,2] == "1/2" | b_rep[,2] == "1|2"]
		# het sites shifted is everything moved over one spot to find ranges
		het_sites_shifted <- c(0,het_sites2[1:(length(het_sites2) - 1)])
		
		# measure distances between het sites
		het_diff <- het_sites2 - c(1, het_sites2[1:(length(het_sites2)-1)]) - 1
		# check how many SNPs are after the last site
		het_diff_end <- nrow(b_rep) - het_sites2[length(het_sites2)]
		
		# ROH start and end locations in the b_rep matrix
		ROH_start <- het_sites_shifted[het_diff >= ROH_cutoff] + 1
		ROH_end <- het_sites2[het_diff >= ROH_cutoff] - 1
		# add the end of the scaffold if necessary
		if(het_diff_end >= ROH_cutoff) {
			ROH_start <- c(ROH_start, het_sites2[length(het_sites2)] + 1)
			ROH_end <- c(ROH_end, nrow(b_rep))
		}
		
		# obtain output values
		individual_rep <- rep(ind_names[b], length(ROH_start))
		scaffold_rep <- rep(strsplit(x_files[a], "_")[[1]][1], length(ROH_start))
		start_rep <- b_rep[ROH_start, 1]
		end_rep <- b_rep[ROH_end, 1]
		nsnps_rep <- ROH_end - ROH_start + 1
		output_rep <- cbind(individual_rep, scaffold_rep, start_rep, end_rep, nsnps_rep)
		write.table(output_rep, file="runs_of_homozygosity.txt", sep="\t", quote=F, col.names=F, row.names=F, append=T)
		
	}
}
	


# summarize the results
x <- read.table("runs_of_homozygosity.txt", sep="\t", header=T, stringsAsFactors=F)

x_names <- unique(x[,1])

output <- c()
# summarize for each individual
for(a in 1:length(x_names)) {
	a_rep <- x[x[,1] == x_names[a], ]
	num_ROH <- nrow(a_rep)
	length_ROH <- sum(a_rep$end - a_rep$start + 1)
	num_SNPs <- sum(a_rep$n_snps)
	output <- rbind(output, c(x_names[a], num_ROH, length_ROH, num_SNPs))
}
output <- data.frame(individual=as.character(output[,1]), num_ROH=as.numeric(output[,2]),
					length_ROH=as.numeric(output[,3]), num_SNPs=as.numeric(output[,4]))
write.table(output, file="zosterops_ROH_summary.txt", sep="\t", row.names=F, quote=F)



