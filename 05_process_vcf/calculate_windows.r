# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# file name for parsing chromosome and window bounds
filename_simple <- strsplit(args[1], "windows/")[[1]][2]

# add functions for calculations
source("window_stat_calculations.r")

# define minimum number of sites to keep a fasta file 
min_sites <- 10000

# no scientific notation
options(scipen=999)

# read in input file
input_file <- read.table(args[1], sep="\t", stringsAsFactors=F)

# subset input file 
input_file_genotypes <- input_file[,4:ncol(input_file)]

# read in populations
populations <- read.table(args[2], sep="\t", stringsAsFactors=F, header=T)

# define output name
output_name <- paste(strsplit(args[1], ".simple")[[1]][1], "__stats.txt", sep="")
# define output fasta name
output_fasta <- paste(strsplit(args[1], ".simple")[[1]][1], ".fasta", sep="")

# write output file
write(c("pop1", "pop2", "stat", "chr", "start", "end", "number_sites", "number_variable_sites", "calculated_stat"), ncolumns=9, file=output_name, sep="\t")

# calculate heterozygosity for each individual
heterozygosity(input_file, populations, output_name, filename_simple)

# calculate transition / transversion ratio
titv(input_file, populations, output_name, filename_simple)

# create fasta sequence alignments for all files with a minimum number of sites
create_fasta_from_vcf(input_file, populations, output_fasta, filename_simple, min_sites)
