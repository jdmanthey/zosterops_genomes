ref <- read.table("ref.fa.fai", sep="\t", stringsAsFactors=F)

output <- c()
for(a in 1:nrow(ref)) {
	pos1 <- sample(1:10000, 1)
	pos2 <- sample(31000:40000, 1)
	strand <- sample(c("+", "-"), 1)
	a_rep <- c(ref[a,1], pos1, pos2, paste("fake_gene_", a, sep=""), 0, strand, pos1, pos2, 0, 1, paste(pos2-pos1, ",", sep=""), "0,")
	output <- rbind(output, a_rep)
}

write.table(output, file="fake_genes.bed", sep="\t", quote=F, row.names=F, col.names=F)

