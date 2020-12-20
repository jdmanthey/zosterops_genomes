module load intel R

R

options(scipen=999)

# list all the files in the trees directory
x_files <- list.files(pattern="*tre")

# find the chromosome, start, and end for each tree
x_names <- list.files(pattern="*tre")
x_chrom <- sapply(strsplit(sapply(strsplit(x_names, "RAxML_bipartitions."), "[[", 2), ":"), "[[", 1)
x_start <- sapply(strsplit(sapply(strsplit(x_names, ":"), "[[", 2), "-"), "[[", 1)
x_end <- sapply(strsplit(sapply(strsplit(x_names, "-"), "[[", 2), ".tre"), "[[", 1)

# write tree info
write.table(cbind(x_chrom, x_start, x_end), file="tree_info.txt", sep="\t", quote=F, row.names=F, col.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="zost_trees_25kbp.trees", ncolumns=1)

