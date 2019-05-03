# compare the main summary tree from the overall signal as well as ASTRAL, and compare all of the
# sliding window gene trees to that summary tree

require(ape)

species_tree <- read.nexus("zost_summed.tre")
plot(species_tree)
species_tree <- root(species_tree, outgroup="Z_japonicus_11220", resolve.root=T)

# read in all gene trees and notes about those trees and alignments
total_trees <- read.tree("zosterops_total.trees")
# order of trees
tree_notes <- read.table("zosterops_total_notes.txt", sep="", stringsAsFactors=F)
# summaries of fastas
fasta_summary <- read.table("zosterops_fasta_summary.txt", sep="", header=T, stringsAsFactors=F)
	
	
	
# loop for each tree to measure distance and summarize
output2 <- c()
output3 <- c()
for(a in 1:nrow(fasta_summary)) {
	# read in summary data
	a_output <- fasta_summary[a,]
	# read in gene tree
	a_tree <- total_trees[tree_notes[tree_notes[,2] == a_output[1,1] & tree_notes[,3] == a_output[1,2],1]]
	# make sure it is rooted
	a_tree <- root(a_tree, outgroup="Z_japonicus_11220", resolve.root=T)
	# measure the Kuhner and Felsenstein (1994) branch length score
	a_kf94 <- dist.topo(species_tree, a_tree, method="score")[1]
	a_ph85 <- dist.topo(species_tree, a_tree, method="PH85")[1]
	# cat the output
	output2 <- c(output2, a_kf94)
	output3 <- c(output3, a_ph85)
}	
kf94 <- output2
ph85 <- output3

# combine tree distance with other summary info
output <- cbind(fasta_summary, kf94, ph85)

# rearrange the output based on chromosome names
output <- rbind(output[output[,1] == "CHR1", ],
				output[output[,1] == "CHR1A", ],
				output[output[,1] == "CHR1B", ],
				output[output[,1] == "CHR2", ],
				output[output[,1] == "CHR3", ],
				output[output[,1] == "CHR4", ],
				output[output[,1] == "CHR4A", ],
				output[output[,1] == "CHR5", ],
				output[output[,1] == "CHR6", ],
				output[output[,1] == "CHR7", ],
				output[output[,1] == "CHR8", ],
				output[output[,1] == "CHR9", ],
				output[output[,1] == "CHR10", ],
				output[output[,1] == "CHR11", ],
				output[output[,1] == "CHR12", ],
				output[output[,1] == "CHR13", ],
				output[output[,1] == "CHR14", ],
				output[output[,1] == "CHR15", ],
				output[output[,1] == "CHR17", ],
				output[output[,1] == "CHR18", ],
				output[output[,1] == "CHR19", ],
				output[output[,1] == "CHR20", ],
				output[output[,1] == "CHR21", ],
				output[output[,1] == "CHR22", ],
				output[output[,1] == "CHR23", ],
				output[output[,1] == "CHR24", ],
				output[output[,1] == "CHR25", ],
				output[output[,1] == "CHR26", ],
				output[output[,1] == "CHR27", ],
				output[output[,1] == "CHR28", ],
				output[output[,1] == "CHRLGE22", ],
				output[output[,1] == "CHRZ", ])

# rename row names
rownames(output) <- seq(from=1, to=nrow(output), by=1)

# find sliding window means of the distance metrics
# 500 kbp sliding windows
sliding_kf94 <- c()
sliding_ph85 <- c()
for(a in 1:nrow(output)) {
	a_range <- seq(from=(a - 4), to=(a + 5), by = 1)
	# don't go over the limits of the genome!
	a_range <- a_range[a_range > 0 & a_range <= nrow(output)]
	# only windows on same chromosome
	a_range <- a_range[output[a_range,1] %in% output[a,1]]
	sliding_kf94 <- c(sliding_kf94, mean(output$kf94[a_range]))
	sliding_ph85 <- c(sliding_ph85, mean(output$ph85[a_range]))
}



# what are the unique chromosomes
chr <- unique(output[,1])
chr_polygons_ph <- list()
chr_polygons_kf <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- rownames(output[output[,1] == chr[a],])
	a2 <- as.numeric(a1[length(a1)])
	a1 <- as.numeric(a1[1])
	chr_polygons_ph[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 30), c(a1, 30), c(a1, 0))
	chr_polygons_kf[[a]] <- rbind(c(a1, 0.002), c(a2, 0.002), c(a2, 0.008), c(a1, 0.008), c(a1, 0.002))	
}

par(mar=c(0.5,5,0.5,0))
par(mfrow=c(2,1))
plot(c(-1,-1), ylim=c(0,30), xlim=c(1,length(output$ph85)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="PH85 Distance")
odd <- 0
for(a in 1:length(chr_polygons_ph)) {
	if(odd == 1) {
	polygon(chr_polygons_ph[[a]], col="lightsteelblue1", border="white")
	odd <- 0	
	} else {
		odd <- 1
	}
}
points(output$ph85, pch=19, cex=0.02)
lines(sliding_ph85, lwd=0.7)


plot(c(-1,-1), ylim=c(0.002,0.008), xlim=c(1,length(output$kf94)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="KF94 Distance")
odd <- 0
for(a in 1:length(chr_polygons_kf)) {
	if(odd == 1) {
	polygon(chr_polygons_kf[[a]], col="lightsteelblue1", border="white")
	odd <- 0	
	} else {
		odd <- 1
	}
}
points(output$kf94, pch=19, cex=0.02)
lines(sliding_kf94, lwd=0.7)


