# create msmc files for each individual, one per chromosome
library(data.table)

# order of individuals in vcf: 
# 1	10	11	12	13	14	15	16	17	2	3	4	5	6	7	8	9

# set up names of individuals
ind_names <- c("Z_griseotinctus_2003067", "Z_murphyi_CEF838", "Z_rendovae_33262", "Z_splendidus_CEF825b",
	"Z_stresemanni_32763", "Z_teteparius_CES726", "Z_ugiensis_32795", "Z_ugiensis_35017", "Z_vellalavella_CEF819",
	"Z_japonicus_11220", "Z_kulambangrae_CEF484", "Z_kulambangrae_VGR592", "Z_luteirostris_PRS2615", "Z_metcalfi_32056",
	"Z_metcalfi_32096", "Z_metcalfi_34483", "Z_metcalfi_34869")

# list vcf files
x_files <- list.files(pattern="*vcf")
x_names <- sapply(strsplit(x_files, "_"), "[[", 1)

out_dir <- "zosterops_demography"
dir.create(out_dir)


# loop for each chromosome
for(a in 1:length(x_files)) {
	a_rep <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
	for(b in 1:length(ind_names)) {
		if(b %in% c(2:9,11:17)) {
			if(a == 1) {
				dir.create(paste(out_dir, "/", ind_names[b], sep=""))
			}
			# subset locations and variants for this individual
			b_rep <- a_rep[,c(1,2,3,b+3)]
			b_rep <- b_rep[b_rep[,4] != "./.",]
			# replace all phased genotypes to keep it simple
			b_rep[,4] <- gsub("\\|", "/", b_rep[,4])
			
			# subset only the variants
			b_rep_matches <- b_rep[,4] == "0/0" | b_rep[,4] == "1/1" | b_rep[,4] == "2/2" | b_rep[,4] == "3/3" | b_rep[,4] == "4/4" | 
								b_rep[,4] == "5/5" | b_rep[,4] == "6/6" | b_rep[,4] == "7/7" | b_rep[,4] == "8/8" | b_rep[,4] == "9/9"
			b_rep_matches <- b_rep_matches == FALSE
			b_rep_variant <- b_rep[b_rep_matches, ]
			# remove the null alleles asterisks
			b_rep_variant[,3] <- gsub("\\*,", "", b_rep_variant[,3])
			b_rep_variant[,3] <- gsub(",\\*", "", b_rep_variant[,3])
			
			# determine number of sites before each match
			b_rep_sites <- seq(from=1, to=nrow(b_rep), by=1)
			b_rep_sites <- c(0, b_rep_sites[b_rep_matches])
			b_rep_sites <- diff(b_rep_sites)
			
			# output
			output <- cbind(rep(x_names[a], length(b_rep_sites)),
							b_rep_variant[,1],
							b_rep_sites,
							paste(b_rep_variant[,2], substr(b_rep_variant[,3], 1, 1), sep=""))
			write.table(output, file=paste(out_dir, "/", ind_names[b], "/", x_names[a], ".txt", sep=""), sep="\t",
							quote=F, row.names=F, col.names=F)
			}
	}
}




