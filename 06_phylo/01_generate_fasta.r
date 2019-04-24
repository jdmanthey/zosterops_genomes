x_files <- list.files(pattern="*phylo.simple.vcf")
x_names <- sapply(strsplit(x_files, "_"), "[[", 1)

# window size
window_size <- 50000
# minimum content to keep the window (e.g., 0.5 = half the window size)
keep_window <- 0.5
output_directory <- "zosterops_fasta"
dir.create(output_directory)

# order of individuals in vcf: 
# 1	10	11	12	13	14	15	16	17	2	3	4	5	6	7	8	9

# set up names of individuals
ind_names <- c("Z_griseotinctus_2003067", "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b",
	"Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819",
	"Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056",
	"Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")
ind_fasta_names <- paste(">", ind_names, sep="")


# loop for each chromosome
for(a in 1:length(x_files)) {
	a_start <- 1
	a_end <- window_size
	a_rep <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
	# remove non bi-allelic SNPs (how to code so many possibilities and null alleles?)
	a_rep <- a_rep[nchar(a_rep[,3]) == 1, ]
	a_window_number <- ceiling(max(a_rep[,1]) / window_size)
	
	# loop for each window
	for(b in 1:a_window_number) {
		b_rep <- a_rep[a_rep[,1] >= a_start & a_rep[,1] <= a_end, ]
		if(nrow(b_rep) > (window_size * keep_window)) {
			allele1 <- b_rep[,2]
			allele2 <- b_rep[,3]
			# set up heterozygous codes
			heterozygous <- rep("?", length(allele1))
			heterozygous[allele1 == "A" & allele2 == "C"] <- "M"
			heterozygous[allele1 == "C" & allele2 == "A"] <- "M"
			heterozygous[allele1 == "A" & allele2 == "G"] <- "R"
			heterozygous[allele1 == "G" & allele2 == "A"] <- "R"
			heterozygous[allele1 == "A" & allele2 == "T"] <- "W"
			heterozygous[allele1 == "T" & allele2 == "A"] <- "W"
			heterozygous[allele1 == "C" & allele2 == "G"] <- "S"
			heterozygous[allele1 == "G" & allele2 == "C"] <- "S"
			heterozygous[allele1 == "C" & allele2 == "T"] <- "Y"
			heterozygous[allele1 == "T" & allele2 == "C"] <- "Y"
			heterozygous[allele1 == "G" & allele2 == "T"] <- "K"
			heterozygous[allele1 == "T" & allele2 == "G"] <- "K"
				
			# loop for each individual
			for(c in 1:length(ind_names)) {
				c_rep <- b_rep[,c+3]
				c_rep[c_rep == "./."] <- "?"
				c_rep[c_rep == "0/0"] <- allele1[c_rep == "0/0"]
				c_rep[c_rep == "0|0"] <- allele1[c_rep == "0|0"]
				c_rep[c_rep == "1/1"] <- allele2[c_rep == "1/1"]
				c_rep[c_rep == "1|1"] <- allele2[c_rep == "1|1"]
				c_rep[c_rep == "0/1"] <- heterozygous[c_rep == "0/1"]
				c_rep[c_rep == "0|1"] <- heterozygous[c_rep == "0|1"]
				c_rep[c_rep == "1|0"] <- heterozygous[c_rep == "1|0"]
				c_rep <- paste(c_rep, collapse="", sep="")
				if(c == 1) {
					write(ind_fasta_names[c], paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1)
					write(c_rep, paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
				} else {
					write(ind_fasta_names[c], paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
					write(c_rep, paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
				}
			}
		}
		a_start <- a_start + window_size
		a_end <- a_end + window_size
	}
	
	
	
}



