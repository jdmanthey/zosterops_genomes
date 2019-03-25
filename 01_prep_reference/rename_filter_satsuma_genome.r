require(Biostrings)
genome <- readDNAStringSet("pseudochromosomes.fasta")

# filter scaffolds >= 1 Mbp
keep <- genome@ranges@width >= 1000000
genome <- genome[keep]

# rename the chromosomes to simpler names
genome_names <- genome@ranges@NAMES
genome_names <- paste("CHR", sapply(strsplit(sapply(strsplit(as.character(genome_names), "_"), "[[", 8), ","), "[[", 1), sep="")
genome@ranges@NAMES <- genome_names

# write output
writeXStringSet(genome, file="ref.fa")
