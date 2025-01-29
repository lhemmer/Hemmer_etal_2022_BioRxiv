<h1>GATK Variant Analysis<h1>

<h2>Note<h2> 

The VCF files as the output of GATK and used for analysis are found on Dryad.


<h2>Scripts to run GATK and process<h2> 

1. Run GATK for all files separately: GATK.SNP.analysis_new_calling.sh
	- The intervals.select.list file is to run GATK on the regions in the genome, either specifically or all
	- Run with the fastq file per sample as ```GATK.SNP.analysis_new_calling.sh samples.fastq```
	- The DGRP Samples are from BioProject PRJNA36679

2. Combine all vcfs into one and filter: GATK.SNP.analysis_new_calling_2half.sh
	- Just run in the same directory as your outputs from above
	
3. Calculate Tajima's D: tajima.slurm
	- Runs ```SimplifyVCF_Basic.v2.pl``` followed by ```Windows_Basic.v2.pl``` to calculate Tajima's D in several windows
	
4. Calculate pi with Pixy
	- Example of running pixy for each vcf
	- ```pixy --stats pi --vcf All_filtered_snps.select.Contig119.vcf --populations DGRP.vcf.list.txt --window_size 500 --n_cores 20 --chromosomes 'Contig119' --output_prefix All_filtered_snps.select.Contig119```

5. Plots of pi and Tajima's D: Figure_S8-9_pi_tajimasD_DGRP3.R

6. Analyze Tajima's D output: Fig S9, S10, and S11 R files
	- VCF files were modified to only include the genome contig for each analysis
	- The variants or SNPs in the satellite DNA regions were also removed for some analysis

