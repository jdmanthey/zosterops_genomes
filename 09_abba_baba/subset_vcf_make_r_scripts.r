# subset the vcf files into sliding windows and keep a record of which file is which

options(scipen=999)

file_number <- 1

# list the files and get the chromosome names
x_files <- list.files(pattern="*abba.simple.vcf")
x_names <- sapply(strsplit(x_files, "_"), "[[", 1)

# make window directory
dir.create("window_vcf")

# define window size
window_size <- 100000

# output tracking file
output_tracker <- c()

for(a in 1:length(x_files)) {
	print(a)
	a_start <- 1
	a_end <- window_size
	a_rep <- read.table(x_files[a], stringsAsFactors=F, sep="\t")
	# replace all phased haplotypes with the regular /
	for(b in 4:ncol(a_rep)) {
		a_rep[,b] <- gsub("\\|", "/", a_rep[,b])
	}
	
	# define the number of windows
	a_windows <- floor(a_rep[nrow(a_rep), 1] / window_size)
	
	# loop for each window
	for(b in 1:a_windows) {
		b_rep <- a_rep[a_rep[,1] >= a_start & a_rep[,1] <= a_end, ]
		write.table(b_rep, file=paste("window_vcf/", file_number, ".txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
		
		# update the output tracking file
		b_output <- c(x_names[a], a_start, a_end, file_number)
		output_tracker <- rbind(output_tracker, b_output)
		
		# move the sliding window
		a_start <- a_start + window_size
		a_end <- a_end + window_size
		# update the file number
		file_number <- file_number + 1
	}
	
}

write.table(output_tracker, file="abba_windows_info.txt", sep="\t", quote=F, row.names=F, col.names=F)


# write r scripts for submitting jobs
for(a in 1:file_number) {
	a_name <- paste("window_vcf/abba_", a, ".r", sep="")
	a_rep <- 'source("01_abba.r")'
	write(a_rep, a_name, ncolumns=1)
	a_rep <- paste('abba_kulamb_murphyi("', a, '.txt", ', a, ")", sep="")
	write(a_rep, a_name, ncolumns=1, append=T)
}
















