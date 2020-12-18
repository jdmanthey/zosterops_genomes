options(scipen=999)
x <- read.table("window_titv.txt", sep="\t", stringsAsFactors=F, header=T)
# remove sex chromosomes
x <- x[x$chr != "CHR4A" & x$chr != "CHRZ",]
# calculate mean transition / transversion rate across all variable sites
sum(x$number_variable_sites * x$calculated_stat) / sum(x$number_variable_sites)
# value = 2.2351

# calculate mean heterozygosity per individual
x <- read.table("window_heterozygosity.txt", sep="\t", stringsAsFactors=F, header=T)
x <- x[x$number_sites > 0, ]
output <- c()
for(a in 1:length(unique(x[,1]))) {
	a_rep <- x[x$pop1 == unique(x[,1])[a], ]
	# remove chr 4A (neo sex chromosome, with some individuals female having excess het.)
	a_rep <- a_rep[a_rep$chr != "CHR4A",]
	a_div <- sum(a_rep$number_sites * a_rep$calculated_stat) / sum(a_rep$number_sites)
	
	a_sd <- sd(a_rep$calculated_stat)
	output <- rbind(output, c(a_rep[1,1], a_div, a_sd))
}
# write to table
output <- data.frame(ID=as.character(output[,1]), obs_het=as.numeric(output[,2]), het_sd=as.numeric(output[,3]))
write.table(output, file="zost_het.txt", sep="\t", quote=F, row.names=F)



# plot diversity across the genome for each individual (those from Solomon Islands)
x <- read.table("window_heterozygosity.txt", sep="\t", stringsAsFactors=F, header=T)

inds <- unique(x[,1])
inds <- inds[inds != "Z_griseotinctus_2003067" & inds != "Z_japonicus_11220"]

# list of chromosomes
chr <- unique(x$chr)
# reorder (w/o 4A again)
chr <- chr[c(10,11,12,22,23,24,26,27,28,29,1,2,3,4,5,6,7,8,9,13,14,15,16,17,18,19,20,21,31,32)]


# define window size from file
window_size <- x[1,6] - x[1,5] + 1
# define number of windows for line plots
num_windows <- 10

par(mfrow=c(3,5))
par(mar=c(0,0,0,0))


# loop for each individual
for(a in 1:length(inds)) {
	a_rep <- x[x[,1] == inds[a],]
	
	# reorder by chromosome and position
	new_a_rep <- c()
	for(b in 1:length(chr)) {
		b_rep <- a_rep[a_rep$chr == chr[b], ]
		b_rep <- b_rep[order(b_rep$start), ]
		# add to new object
		new_a_rep <- rbind(new_a_rep, b_rep)
	}
	# rename rows 
	rownames(new_a_rep) <- seq(from=1, to=nrow(new_a_rep), by=1)
	
	# what are the unique chromosomes and their bounding areas for plotting?
	total_windows <- nrow(new_a_rep)
	chr_polygons <- list()
	# make the plotting polygons
	for(d in 1:length(chr)) {
		a1 <- rownames(new_a_rep)[new_a_rep[,4] == chr[d]]
		a2 <- a1[length(a1)]
		a1 <- a1[1]
		chr_polygons[[d]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 0.03), c(a1, 0.03), c(a1, 0))
	}
	
	plot(c(-1,-1), ylim=c(0,0.03), xlim=c(1, nrow(new_a_rep)), xaxt="n", yaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="")
	odd <- 0
	for(d in 1:length(chr_polygons)) {
		if(odd == 1) {
			polygon(chr_polygons[[d]], col="snow2", border="white")
			odd <- 0	
		} else {
			odd <- 1
		}
	}
	# can't plot points for each window (too many points)
	
	# sliding windows
	x_axis_rep <- seq(from=1, to=total_windows, by=num_windows)
	y_axis_rep <- c()
	scaffold <- c()
	for(d in 1:length(x_axis_rep)) {
		# pull out the rows for the window
		b_rows <- seq(from=(x_axis_rep[d] - as.integer(num_windows / 2) + 1), to=(x_axis_rep[d] + as.integer(num_windows / 2)))
		b_rows <- new_a_rep[b_rows[b_rows > 0 & b_rows <= total_windows],]
		# keep only those on the same scaffold
		b_rows <- b_rows[b_rows[,1] == new_a_rep[x_axis_rep[d], 1],]
		# get the mean of the stat
		y_axis_rep <- c(y_axis_rep, mean(b_rows$calculated_stat))
		# record the scaffold
		scaffold <- c(scaffold, new_a_rep[x_axis_rep[d],4])
	}
	
	# plot each scaffold at a time (so lines don't connect between scaffolds)
	for(d in 1:length(unique(scaffold))) {
		a_rep <- cbind(x_axis_rep, y_axis_rep, scaffold)
		a_rep <- a_rep[a_rep[,3] == unique(scaffold)[d], ]
		lines(a_rep[,1:2], lwd=1.2, col="navyblue")
	}

}

















