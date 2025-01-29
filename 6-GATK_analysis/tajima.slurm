#!/bin/bash
#
#SBATCH -p charlesworths
#SBATCH --job-name="gatk.test"
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 4
#SBATCH --mem=100G
#SBATCH --time=0-48:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

###Load modules

module load perl

#zcat All_filtered_snps_highcov.vcf.gz > All_filtered_snps_highcov.vcf

perl SimplifyVCF_Basic.v2.pl All_filtered_snps_highcov.vcf

perl Windows_Basic.v2.pl All_filtered_snps_highcov.simple.vcf DGRP dummy.gff 30 30 1000 10000 DGRP.test.10000.csv

perl Windows_Basic.v2.pl All_filtered_snps_highcov.simple.vcf DGRP dummy.gff 30 30 1000 1000 DGRP.test.1000.csv

perl Windows_Basic.v2.pl All_filtered_snps_highcov.simple.vcf DGRP dummy.gff 30 30 1000 500 DGRP.test.500.csv
