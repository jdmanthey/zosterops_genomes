options(scipen=999)
# read in heterozygosity table
x <- read.table("window_heterozygosity.txt", sep="\t", stringsAsFactors=F, header=T)

# define individuals
inds <- unique(x[,1])

# identify window size
win_size <- x$end[1] - x$start[1] + 1

# keep only windows w/ no heterozygosity and >=20000 sites genotyped
x <- x[x$number_sites >= 20000 & x$calculated_stat == 0, ]

# loop for each individual
output <- c()
for(a in 1:length(inds)) {
	a_rep <- x[x[,1] == inds[a], ]
	if(nrow(a_rep) > 0) { # if there are some ROH
		# loop for each unique chromosome with ROH
		a_chr <- unique(a_rep$chr)
		# remove sex chromosomes
		a_chr <- a_chr[a_chr != "CHR4A" & a_chr != "CHRZ"]
		for(b in 1:length(a_chr)) {
			b_rep <- a_rep[a_rep$chr == a_chr[b], ]
			# order windows by position
			b_rep <- b_rep[order(b_rep$start),]
			# find differences between window starts
			b_diff <- c(1, diff(b_rep$start))
			# make a new data frame with merged ROH in neighboring windows
			for(d in 1:nrow(b_rep)) {
				if(d == 1) {
					d_rep <- t(as.matrix(c(inds[a], b_rep$chr[1], b_rep$start[1], b_rep$end[1])))
				} else {
					if(b_diff[d] == win_size) { # change end site of window if neighboring the previous one
						d_rep[nrow(d_rep),4] <- b_rep$end[d]
					} else { # or add to matrix if not neighboring previous ROH
						d_rep <- rbind(d_rep, c(inds[a], b_rep$chr[d], b_rep$start[d], b_rep$end[d]))
					}
				}
			}
			output <- rbind(output, d_rep)
		}
	}
}

output <- data.frame(ind=as.character(output[,1]), chr=as.character(output[,2]), start=as.numeric(output[,3]), end=as.numeric(output[,4]))


# summarize
#loop for each individual
sum_output <- c()
for(a in 1:length(inds)) {
	a_rep <- output[output$ind == inds[a],]
	if(nrow(a_rep) > 0) {
		length_roh <- sum(a_rep$end - a_rep$start + 1)
		sum_output <- rbind(sum_output, c(inds[a], nrow(a_rep), length_roh))		
	} else {
		sum_output <- rbind(sum_output, c(inds[a], 0, 0))
	}
}
sum_output <- data.frame(ind=as.character(sum_output[,1]), n_ROH=as.numeric(sum_output[,2]), len_ROH=as.numeric(sum_output[,3]))
sum_output











