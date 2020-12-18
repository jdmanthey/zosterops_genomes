options(scipen=999)

individuals <- c("Z_kulambangrae_CEF484","Z_kulambangrae_VGR592","Z_luteirostris_PRS2615","Z_metcalfi_32056","Z_metcalfi_32096","Z_metcalfi_34483","Z_metcalfi_34869","Z_murphyi_CEF838","Z_rendovae_33262","Z_splendidus_CEF825b","Z_stresemanni_32763","Z_teteparius_CES726","Z_ugiensis_32795","Z_ugiensis_35017","Z_vellalavella_CEF819")

# read in diversity table
div <- read.table("zost_het.txt", header=T, stringsAsFactors=F)
# reorder diversity table
div <- div[match(individuals, div[,1]), ]

# half the heterozygosity for the MSMC r parameter
div$obs_het <- div$obs_het / 2

directory <- "/lustre/scratch/jmanthey/10_zosterops/zosterops_demography"

queue <- "omni"
cluster <- "quanah"

# determine all the input directories
directories <- c()
for(a in 1:length(individuals)) {
	directories <- c(directories, paste(directory, "/", individuals[a], "/*txt", sep=""))
	# bootstraps
	for(b in 1:10) {
		directories <- c(directories, paste(directory, "/", individuals[a], "/bootstrap_", b, "/*txt", sep=""))
	}
}

# determine the output names for each 
outputs <- substr(sapply(strsplit(directories, directory), "[[", 2), 2, nchar(sapply(strsplit(directories, directory), "[[", 2)) - 5)
outputs <- gsub("/bootstrap_", "_b", outputs)

# determine the heterozygosities to use for each 
heterozygosities <- div[match(sort(rep(individuals, 11)), div[,1]),2]


# write helper files
write(directories, file="helper1.txt", ncolumns=1)
write(outputs, file="helper2.txt", ncolumns=1)
write(heterozygosities, file="helper3.txt", ncolumns=1)

# write the array job
a_script <- "01_zost_msmc_array.sh"
write("#!/bin/sh", file=a_script)
write("#$ -V", file=a_script, append=T)
write("#$ -cwd", file=a_script, append=T)
write("#$ -S /bin/bash", file=a_script, append=T)
write(paste("#$ -N ", "zost_dem", sep=""), file=a_script, append=T)
write(paste("#$ -q ", queue, sep=""), file=a_script, append=T)
write("#$ -pe sm 2", file=a_script, append=T)
write(paste("#$ -P ", cluster, sep=""), file=a_script, append=T)
write("#$ -l h_rt=48:00:00", file=a_script, append=T)
write("#$ -l h_vmem=8G", file=a_script, append=T)
write(paste("#$ -t 1:", length(directories), sep=""), file=a_script, append=T)
write("", file=a_script, append=T)

write("input_array=$( head -n${SGE_TASK_ID} helper1.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("output_array=$( head -n${SGE_TASK_ID} helper2.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("heter_array=$( head -n${SGE_TASK_ID} helper3.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)

msmc_command <- "~/msmc2_linux64bit -o ${output_array} -i 20 -t 2 -m ${heter_array} -p 1*2+20*1+1*2+1*3 ${input_array}"
write(msmc_command, file=a_script, append=T)








