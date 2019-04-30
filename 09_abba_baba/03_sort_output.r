# summarize the thousands of abba baba outputs

options(scipen=999)

# read in the file about which files are which
x_summary <- read.table("abba_windows_info.txt", stringsAsFactors=F)

# chromosome order
chromosomes <- c("CHR1", "CHR1A", "CHR1B", "CHR2", "CHR3", "CHR4", "CHR4A", "CHR5", "CHR6", "CHR7", "CHR8", "CHR9",
					"CHR10", "CHR11", "CHR12", "CHR13", "CHR14", "CHR15", "CHR17", "CHR18", "CHR19", "CHR20", "CHR21", 
					"CHR22", "CHR23", "CHR24", "CHR25", "CHR26", "CHR27", "CHR28", "CHRLGE22", "CHRZ")
					
# list output files
x_files <- list.files("output", full.names=T)

# check how many comparisons there were by reading the first file
x_test <- read.table(x_files[1], header=T)
num_tests <- nrow(x_test)
test_names <- c()
for(a in 1:num_tests) {
	a_temp <- paste(substr(sapply(strsplit(as.character(as.matrix(x_test[a,1:3])), "_"), "[[", 2), 1, 5),
			sapply(strsplit(as.character(as.matrix(x_test[a,1:3])), "_"), "[[", 3), sep="")
	a_temp <- paste(a_temp[1], "_", a_temp[2], "_", a_temp[3], "_abba.txt", sep="")
	test_names <- c(test_names, a_temp)
}


# loop for each chromosome
for(a in 1:length(chromosomes)) {
	a_windows <- x_summary[x_summary[,1] == chromosomes[a], ]
	# loop for each window
	for(b in 1:nrow(a_windows)) {
		b_rep <- read.table(paste("output/abba_", a_windows[b,4], ".txt", sep=""), header=T, stringsAsFactors=F)
		# loop for each test
		for(c in 1:length(test_names)) {
			if(a == 1 & b == 1) {
				write(c("chromosome", "start", "end", "window", "ind1", "ind2", "ind3", "n_snps", "d", "z_score"), file=test_names[c], ncolumns=10, sep="\t")
			} else {
				write(c(as.character(a_windows[b,]), as.character(b_rep[c,])), file=test_names[c], ncolumns=10, sep="\t", append=T)
			}
		}
	}
}
 

