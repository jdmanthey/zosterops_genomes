require(ape)


x_files <- list.files(pattern="*tre$")
x_chrom <- sapply(strsplit(x_files, "_"), "[[", 1)
x_locations <- as.numeric(sapply(strsplit(sapply(strsplit(x_files, "_"), "[[", 2), ".tre"), "[[", 1))

# read in and root the trees
for(a in 1:length(x_files)) {
	a_rep <- read.tree(file=x_files[a])
	a_rep <- root(a_rep, outgroup="Z_japonicus_11220", resolve.root=T)
	if(a == 1) {
		output_trees <- a_rep
	} else {
		output_trees <- c(output_trees, a_rep)
	}
}

# write the trees to a new file
write.tree(output_trees, file="zosterops_total.trees")
write.table(cbind(seq(from=1, to=length(x_files), by=1), x_chrom, x_locations), file="zosterops_total_notes.txt", quote=F, row.names=F, col.names=F)












