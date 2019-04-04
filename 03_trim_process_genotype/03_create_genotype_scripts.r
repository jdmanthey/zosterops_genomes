# popmap = individual base names of fastq files, one line per individual

# make sure reference is indexed with bwa and samtools before use, and use CreateSequenceDictionary in GATK 
	
	project_directory <- "/lustre/scratch/jmanthey/10_zosterops"
	directory_name <- "zost_scripts"
	reference_genome_location <- "/lustre/work/jmanthey/zosterops_genome/ref.fa"
	queue <- "omni"
	cluster <- "quanah"
	output_name <- "zost"
	popmap <- "zost_popmap.txt"
	individuals <- read.table(popmap, sep="\t")

	# make directories
	dir.create(directory_name)
	dir.create(paste(directory_name, "/01_gatk_split", sep=""))
	dir.create(paste(directory_name, "/02_gatk_combine", sep=""))
	dir.create(paste(directory_name, "/03_group_genotype", sep=""))

	chr_list_to_genotype <- c("CHR1", "CHR1A", "CHR1B", "CHR2", "CHR3", "CHR4", "CHR4A", "CHR5", "CHR6", "CHR7", "CHR8", "CHR9", "CHR10", "CHR11", "CHR12", "CHR13", "CHR14", "CHR15", "CHR17", "CHR18", "CHR19", "CHR20", "CHR21", "CHR22", "CHR23", "CHR24", "CHR25", "CHR26", "CHR27", "CHR28", "CHRZ", "CHRLGE22")
	
		
	# step 1
	# genotype all individuals using GATK, one array job per chromosome
	for(a in 1:length(chr_list_to_genotype)) {
		a.name <- paste(project_directory, "/01_bam_files/", "${SGE_TASK_ID}", "_final.bam", sep="")
		a.script <- paste(directory_name, "/01_gatk_split/", chr_list_to_genotype[a], ".sh", sep="")
		write("#!/bin/sh", file=a.script)
		write("#$ -V", file=a.script, append=T)
		write("#$ -cwd", file=a.script, append=T)
		write("#$ -S /bin/bash", file=a.script, append=T)
		write(paste("#$ -N ", chr_list_to_genotype[a], "_step1", sep=""), file=a.script, append=T)
		write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
		write("#$ -pe sm 4", file=a.script, append=T)
		write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
		write("#$ -l h_rt=48:00:00", file=a.script, append=T)
		write("#$ -l h_vmem=10G", file=a.script, append=T)
		write(paste("#$ -t 1:", nrow(individuals), sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write("module load intel java", file=a.script, append=T)
		write("", file=a.script, append=T)
		#gatk 4.0
		gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx40g" HaplotypeCaller -R ', reference_genome_location, " -I ", a.name, " -ERC GVCF -O ", project_directory, "/02_vcf/", chr_list_to_genotype[a], "._${SGE_TASK_ID}_.g.vcf", " --QUIET --intervals ", chr_list_to_genotype[a], sep="")
		write(gatk_command, file=a.script, append=T)
	}

	# step 2 
	# combine the vcfs for each chromosome and each individual
	for(a in 1:length(chr_list_to_genotype)) {
		a.script <- paste(directory_name, "/02_gatk_combine/", chr_list_to_genotype[a], ".sh", sep="")
		write("#!/bin/sh", file=a.script)
		write("#$ -V", file=a.script, append=T)
		write("#$ -cwd", file=a.script, append=T)
		write("#$ -S /bin/bash", file=a.script, append=T)
		write(paste("#$ -N ", chr_list_to_genotype[a], "_step2", sep=""), file=a.script, append=T)
		write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
		write("#$ -pe sm 4", file=a.script, append=T)
		write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
		write("#$ -l h_rt=48:00:00", file=a.script, append=T)
		write("#$ -l h_vmem=10G", file=a.script, append=T)
		write(paste("#$ -t 1:", nrow(individuals), sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write("module load intel java", file=a.script, append=T)
		write("", file=a.script, append=T)
		
		#make list of all vcfs
		for(b in 1:nrow(individuals)) {
			if(b == 1) {
				vcf_total <- paste("--variant ", project_directory, "/02_vcf/", chr_list_to_genotype[a], "._", b, "_.g.vcf", sep="")
			} else {
				vcf_total <- paste(vcf_total, " --variant ", project_directory, "/02_vcf/", chr_list_to_genotype[a], "._", b, "_.g.vcf", sep="")
			}
		}
		
		#gatk 4.0
		gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx40g" CombineGVCFs -R ', reference_genome_location, " ", vcf_total, " -O ", project_directory, "/02_vcf/", chr_list_to_genotype[a], ".g.vcf", sep="")
		write(gatk_command, file=a.script, append=T)
	}





	# step 3 
	# group genotype
	for(a in 1:length(chr_list_to_genotype)) {
		a.script <- paste(directory_name, "/03_group_genotype/", chr_list_to_genotype[a], ".sh", sep="")
		write("#!/bin/sh", file=a.script)
		write("#$ -V", file=a.script, append=T)
		write("#$ -cwd", file=a.script, append=T)
		write("#$ -S /bin/bash", file=a.script, append=T)
		write(paste("#$ -N ", chr_list_to_genotype[a], "_step3", sep=""), file=a.script, append=T)
		write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
		write("#$ -pe sm 10", file=a.script, append=T)
		write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
		write("#$ -l h_rt=48:00:00", file=a.script, append=T)
		write("#$ -l h_vmem=10G", file=a.script, append=T)
		write(paste("#$ -t 1:", nrow(individuals), sep=""), file=a.script, append=T)
		write("", file=a.script, append=T)
		write("module load intel java", file=a.script, append=T)
		write("", file=a.script, append=T)
		
		#gatk 4.0
		gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx100g" GenotypeGVCFs -R ', reference_genome_location, " -V ", project_directory, "/02_vcf/", chr_list_to_genotype[a], ".g.vcf --include-non-variant-sites -O ", project_directory, "/03_vcf/", chr_list_to_genotype[a], ".g.vcf", sep="")
		write(gatk_command, file=a.script, append=T)
	}
	
















