Workflow to run a MELT analysis to call polymorphic transposable elements in the Zosterops resequencing data.

Input necessary for all TEs that are going to be analyzed:
1) Consensus sequence
2) .bed file showing location of all TEs closely-related to the TE to be checked for polymorphisms (e.g., all ERVs)
3) a unique format gene annotation file, a psuedo-gene file can be created with the r script in this directory 
if you do not have real annotations
4) the reference genome that was used for alignment (must have associated fai file)

In Solomon Islands' Zosterops, I looked at histograms of divergence to consensus for CR1s and ERVs, two groups of 
retrotransposons with recent activity in birds. ERVs showed a clear peak in recent time, so I looked at what families 
were <= 10% divergent from the consensus in ERVs. Of these, ten ERVs had more than 500 copies in the reference 
genome at low divergence, and these were chosen as candidates for MELT polymorphism analyses.
