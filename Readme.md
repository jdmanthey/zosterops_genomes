
Run steps in directories to get from raw data to output data tables and some of the figures (or partial figures later edited in Illustrator).

01_prep_reference = modify the published Zosterops lateralis genome for use in this project, including annotation of TEs

02_qual_stats = run fastqc and summarize with multiqc on all raw data

03_trim_process_genotype = quality trim with bbduk, samtools to convert to bam, and then GATK for bam processing and genotyping

04_depth = get alignment depth statistics and plots

05_MELT = MELT workflow to call polymorphic transposable elements

05_process_vcf = filter all vcf files, make fasta files for phylogenomics, window calculations, summarize diversity and ROH

06_phylo = phylogenomics of windowed fasta files

07_demography = demographic analyses in MSMC, plotting, and calculating pop. sizes

08_island_size = summarize and plot relationships of island size and demography and genomic diversity

09_plotting_correlations_2024.r = updated plotting and analyses in 2024 (by Ethan Gyllenhaal)

