# measure observed heterozygosity for each individual


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
# columns = individual, scaffold, number heterozygous sites, total number sites
write(c("individual", "scaffold", "het_sites", "total_sites_genotyped"), file="zosterops_diversity_per_scaffold.txt", sep="\t", ncolumns=4)

# loop for each vcf file
for(a in 1:length(x_files)) {
	a_rep <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
	
	# loop for each individual
	for(b in 1:length(ind_names)) {
		b_rep <- a_rep[,b+3]
		# remove non genotyped sites
		b_rep <- b_rep[b_rep != "./."]
		b_total <- length(b_rep)
		# remove 0/0 (most common genotype)
		b_rep <- b_rep[b_rep != "0/0"]
		# remove the phase bar and replace with /
		b_rep <- gsub("\\|", "/", b_rep)
		b_het <- length(b_rep[sapply(strsplit(b_rep, "/"), "[[", 1) != sapply(strsplit(b_rep, "/"), "[[", 2)])
		write(c(ind_names[b], strsplit(x_files[a], "_")[[1]][1], b_het, b_total), file="zosterops_diversity_per_scaffold.txt", sep="\t", ncolumns=4, append=T)
	}
}

