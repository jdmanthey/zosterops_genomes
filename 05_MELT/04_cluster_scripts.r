# popmap = tab delimited file of all individual name (i.e. BAM file names), 2nd column = number, 3rd column = coverage

	popmap <- "popmap.txt"
	melt_location <- "/lustre/work/jmanthey/MELTv2.1.5/MELT.jar"
	project_directory <- "/lustre/scratch/jmanthey/10_zosterops/04_MELT/"
	reference_genome_location <- "/lustre/work/jmanthey/zosterops_genome/ref.fa"
	bam_file_location <- "/lustre/scratch/jmanthey/10_zosterops/01_bam_files/"
	gene_file_location <- "/lustre/work/jmanthey/MELTv2.1.5/reference/02_zost_me_refs/fake_genes.bed"
	TE_bed_file_location <- "/lustre/work/jmanthey/MELTv2.1.5/reference/02_zost_me_refs/MELT_ERV.bed"
	TE_zips <- c("TguERVK7LTR4zLat_MELT.zip","TguERVK8LTR1fzLat_MELT.zip","TguERVL2b5LTRazLat_MELT.zip","TguERVL2b5LTRbzLat_MELT.zip","TguLTR11gzLat_MELT.zip","TguLTR13czLat_MELT.zip","zLatLTR01_MELT.zip","zLatLTR02_MELT.zip","zLatLTR03_MELT.zip","zLatLTR09_MELT.zip")
	TE_zips_location <- "/lustre/work/jmanthey/MELTv2.1.5/reference/02_zost_me_refs/"
	TE_output <- substr(TE_zips, 1, nchar(TE_zips) - 4)
	queue <- "omni"
	cluster <- "quanah"
	runtime <- "48:00:00"
	memory <- "8G"
	
	min_scaffold_size <- as.character("1000000")
	
	individuals <- read.table(popmap, sep="\t", stringsAsFactors=F)
	dir.create("melt_scripts")
	dir.create("melt_scripts/step1")
	dir.create("melt_scripts/step2")
	dir.create("melt_scripts/step3")
	dir.create("melt_scripts/step4")
	#$ -t 1-13

	# set up scripts for individual analyses
	for(a in 1:length(TE_zips)) {
		for(b in 1:nrow(individuals)) {
			a.script <- paste("melt_scripts/step1/", individuals[b,1], "_", a, ".sh", sep="")
			indiv_short <- sapply(strsplit(individuals[b,1], "_"), "[[", 3)
			write("#!/bin/sh", file=a.script)
			write("#$ -V", file=a.script, append=T)
			write("#$ -cwd", file=a.script, append=T)
			write("#$ -S /bin/bash", file=a.script, append=T)
			write(paste("#$ -N a_", indiv_short, "_", a, sep=""), file=a.script, append=T)
			write("#$ -o $JOB_NAME.o$JOB_ID", file=a.script, append=T)
			write("#$ -e $JOB_NAME.e$JOB_ID", file=a.script, append=T)
			write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
			write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
			write("#$ -pe sm 1", file=a.script, append=T)
			write(paste("#$ -l h_rt=", runtime, sep=""), file=a.script, append=T)
			write(paste("#$ -l h_vmem=", memory, sep=""), file=a.script, append=T)
			write("", file=a.script, append=T)
			write("module load intel samtools java", file=a.script, append=T)
			write("", file=a.script, append=T)
			write(paste("mkdir ", project_directory, TE_output[a], sep=""), file=a.script, append=T)
			write("", file=a.script, append=T)
			write(paste("java -Xmx", memory, " -jar ", melt_location, " IndivAnalysis -bamfile ",
				bam_file_location, individuals[b,2], "_final.bam -c ", individuals[b,3], " -d ",
				min_scaffold_size, " -h ", reference_genome_location, " -t ", TE_zips_location, TE_zips[a],
				" -w ", project_directory, TE_output[a], "/", sep=""), file=a.script, append=T)
		}
	}
	
	# scripts for group analysis
	for(a in 1:length(TE_zips)) {
		a.script <- paste("melt_scripts/step2/", "group_", a, ".sh", sep="")
		write("#!/bin/sh", file=a.script)
		write("#$ -V", file=a.script, append=T)
		write("#$ -cwd", file=a.script, append=T)
		write("#$ -S /bin/bash", file=a.script, append=T)
		write(paste("#$ -N group_", a, sep=""), file=a.script, append=T)
		write("#$ -o $JOB_NAME.o$JOB_ID", file=a.script, append=T)
		write("#$ -e $JOB_NAME.e$JOB_ID", file=a.script, append=T)
		write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
		write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
		write("#$ -pe sm 1", file=a.script, append=T)
		write(paste("#$ -l h_rt=", runtime, sep=""), file=a.script, append=T)
		write(paste("#$ -l h_vmem=", memory, sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write("module load intel samtools java", file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("java -Xmx", memory, " -jar ", melt_location, " GroupAnalysis -discoverydir ",
			project_directory, TE_output[a], "/ -h ", reference_genome_location, " -n ", gene_file_location,
			" -t ", TE_zips_location, TE_zips[a], " -w ", project_directory, TE_output[a], "/", sep=""), 
			file=a.script, append=T)
	}
	
# set up scripts for genotype analyses
	for(a in 1:length(TE_zips)) {
		for(b in 1:nrow(individuals)) {
			a.script <- paste("melt_scripts/step3/", individuals[b,1], "_", a, ".sh", sep="")
			indiv_short <- sapply(strsplit(individuals[b,1], "_"), "[[", 3)
			write("#!/bin/sh", file=a.script)
			write("#$ -V", file=a.script, append=T)
			write("#$ -cwd", file=a.script, append=T)
			write("#$ -S /bin/bash", file=a.script, append=T)
			write(paste("#$ -N a_", indiv_short, "_", a, sep=""), file=a.script, append=T)
			write("#$ -o $JOB_NAME.o$JOB_ID", file=a.script, append=T)
			write("#$ -e $JOB_NAME.e$JOB_ID", file=a.script, append=T)
			write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
			write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
			write("#$ -pe sm 1", file=a.script, append=T)
			write(paste("#$ -l h_rt=", runtime, sep=""), file=a.script, append=T)
			write(paste("#$ -l h_vmem=", memory, sep=""), file=a.script, append=T)
			write("", file=a.script, append=T)
			write("module load intel samtools java", file=a.script, append=T)
			write("", file=a.script, append=T)
			write(paste("java -Xmx", memory, " -jar ", melt_location, " Genotype -bamfile ",
				bam_file_location, individuals[b,2], "_final.bam -h ", reference_genome_location, 
				" -p ", project_directory, TE_output[a], "/ -t ", TE_zips_location, TE_zips[a],
				" -w ", project_directory, TE_output[a], "/", sep=""), file=a.script, append=T)
		}
	}

# scripts for making VCFs
	for(a in 1:length(TE_zips)) {
		a.script <- paste("melt_scripts/step4/", "group_", a, ".sh", sep="")
		write("#!/bin/sh", file=a.script)
		write("#$ -V", file=a.script, append=T)
		write("#$ -cwd", file=a.script, append=T)
		write("#$ -S /bin/bash", file=a.script, append=T)
		write(paste("#$ -N group_", a, sep=""), file=a.script, append=T)
		write("#$ -o $JOB_NAME.o$JOB_ID", file=a.script, append=T)
		write("#$ -e $JOB_NAME.e$JOB_ID", file=a.script, append=T)
		write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
		write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
		write("#$ -pe sm 1", file=a.script, append=T)
		write(paste("#$ -l h_rt=", runtime, sep=""), file=a.script, append=T)
		write(paste("#$ -l h_vmem=", memory, sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write("module load intel samtools java", file=a.script, append=T)
		write("", file=a.script, append=T)
		write(paste("java -Xmx", memory, " -jar ", melt_location, " MakeVCF -genotypingdir ",
			project_directory, TE_output[a], "/ -h ", reference_genome_location, " -o ",
			project_directory, TE_output[a], "/ -p ", project_directory, TE_output[a], "/ -t ",
			TE_zips_location, TE_zips[a], " -w ", project_directory, TE_output[a], "/", sep=""), 
			file=a.script, append=T)
	}
	












