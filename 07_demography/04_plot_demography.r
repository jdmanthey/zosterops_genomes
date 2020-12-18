# list all zosterops bootstrap 1 files (to get names of all files)
x_files <- list.files(pattern="*b1.final.txt")
x_files_base <- sapply(strsplit(x_files, "_b1"), "[[", 1)

# define parameters
gen <- 1
mu <- 3.16e-9 * gen # needs to be per generation
min_age <- 000
max_age <- 500000
plotting_age_range <- 1000 # years for x axis


# loop through each output directory
output <- list()
output_bootstraps <- list()
for(a in 1:length(x_files_base)) {
	
	# identify main output file
	a_main <- paste(x_files_base[a], ".final.txt", sep="")	
	
	# identify bootstrap outputs
	a_bootstraps <- paste(x_files_base[a], "_b", seq(from=1, to=10, by=1), ".final.txt", sep="")

	# read in main file
	a_rep <- read.table(a_main, sep="\t", header=T)
	# rearrange main file for plotting lines
	for(d in 1:nrow(a_rep)) {
		if(d == 1) {
			a_rep2 <- rbind(c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		} else {
			a_rep2 <- rbind(a_rep2, c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		}
	}
	a_rep <- a_rep2
	# scale by mutation rate
	a_rep[,1] <- a_rep[,1] / mu #don't multiply by generation because mu is in change per year
	a_rep[,2] <- (1 / a_rep[,2]) / (2 * mu)
	# remove very young and old time frames prone to error
	a_rep <- a_rep[a_rep[,1] >= min_age & a_rep[,1] <= max_age,]
	# scale by plotting age range and pop size range
	a_rep <- a_rep / plotting_age_range
	# add to output list
	output[[a]] <- a_rep
	
	# output for each bootstrap
	output_bootstraps[[a]] <- list()
	for(b in 1:length(a_bootstraps)) {
		b_rep <- read.table(a_bootstraps[b], sep="\t", header=T)
		# rearrange main file for plotting lines
		for(d in 1:nrow(b_rep)) {
			if(d == 1) {
				b_rep2 <- rbind(c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			} else {
				b_rep2 <- rbind(b_rep2, c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			}
		}
		b_rep <- b_rep2
		# scale by mutation rate
		b_rep[,1] <- b_rep[,1] / mu #don't multiply by generation because mu is in change per year
		b_rep[,2] <- (1 / b_rep[,2]) / (2 * mu)
		# remove very young and old time frames prone to error
		b_rep <- b_rep[b_rep[,1] >= min_age & b_rep[,1] <= max_age,]
		# scale by plotting age range and pop size range
		b_rep <- b_rep / plotting_age_range
		# add to output list
		output_bootstraps[[a]][[b]] <- b_rep		
	}
}

# set plot items
a_col <- "darkgreen"
par(mfrow=c(3,5))
par(mar=c(1,1,2,0.2))
for(a in 1:length(output)) {
	plot_name1 <- paste("Z.", sapply(strsplit(x_files[a], "_"), "[[", 2))
	plot_name2 <- sapply(strsplit(x_files[a], "_"), "[[", 3)
	plot(c(-1,1), xlim=c(20, 350), ylim=c(0,1300), pch=19, cex=0.01, log="x", xlab="", ylab="", main="", xaxt="n", yaxt="n")
	title(main=bquote(italic(.(plot_name1)) ~ .(plot_name2)), adj=0,line=0.5)
	if(a == 1 | a == 6 | a == 11) {
		axis(side=2, at=c(0, 200, 400, 600, 800, 1000, 1200), labels=FALSE)
	} else {
		axis(side=2, at=c(0, 200, 400, 600, 800, 1000, 1200), labels=FALSE)
	}
	if(a == 11 | a == 12 | a == 13 | a == 14 | a == 15) {
		axis(side=1, at=c(10, 20, 50, 100, 200), labels=F)		
	} else {
		axis(side=1, at=c(10, 20, 50, 100, 200), labels=F)
	}
	
	# plot bootstraps
	for(b in 1:length(output_bootstraps[[1]])) {
		lines(output_bootstraps[[a]][[b]][,1], output_bootstraps[[a]][[b]][,2], col=a_col, lwd=0.3)
	}
	lines(output[[a]][,1], output[[a]][,2], col=a_col, lwd=3)
}

current_pop <- c()
for(a in 1:length(output)) {
	current_pop <- c(current_pop, output[[a]][1,2])
}

# harmonic mean of pop sizes from most recent to 200k years ago
harmonic_pop <- c()
for(a in 1:length(output)) {
	out_rep <- output[[a]]
	
	# define time series
	time_series <- seq(from=as.integer(out_rep[2,1])+1, to=200, by=1)
	# time series pops
	time_pops <- c()
	for(b in 1:length(time_series)) {
		time_pops <- c(time_pops, out_rep[time_series[b] < out_rep[,1],][1,2])
	}
	# harmonic mean of this individual
	harm_rep <- length(time_pops) / sum((1 / time_pops))
	
	# add to output element
	harmonic_pop <- c(harmonic_pop, harm_rep)
}

harmonic_pop <- data.frame(id=as.character(x_files_base), current_pop=as.numeric(current_pop), harmonic_pop=as.numeric(harmonic_pop))
harmonic_pop



