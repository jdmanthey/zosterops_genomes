
Run steps in directories to get from raw data to output data tables and some of the figures (or partial figures later edited in Illustrator).

01 directory = modify the published Zosterops lateralis genome for use in this project, including annotation of TEs

02 directory = run fastqc and summarize with multiqc on all raw data

03 directory = quality trim with bbduk, samtools to convert to bam, and then GATK for bam processing and genotyping

04 directory = get alignment depth statistics and plots

05 directory = MELT workflow to call polymorphic transposable elements

05 directory = genotype all individuals using GATK

06 directory = subset genome into 50kbp segments and estimate gene trees from each

07 directory = demographic analyses in MSMC

08 directory = estimate observed heterozygosity and runs of homozygosity

09 directory = abba/baba (d statistic) tests
