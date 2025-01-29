#!/bin/bash
#
#SBATCH -p preempt
#SBATCH --job-name="mcc.91"
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH --time=0-48:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

## Load conda environment

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate mcclintock

module load git

line="RAL91" #An example of your sample name, for your fastq files

module load parallel

mkdir output/${line}_all
cd output/${line}_all

ref_dir="~/reference"
run_dir="~/mcclintock"

###Commands

python3 $run_dir/mcclintock.py \
    -r $ref_dir/dmel_scaffold2_plus0310_2.fasta \
    -c $ref_dir/specieslib_mcClintock_2020_mod.fasta \
    -g $ref_dir/dmel.chromosomes.fa.TE.mcClintock.2020.gff \
    -t $ref_dir/dmel.chromosomes.fa.TE.mcClintock.2020.tsv \
    -1 $run_dir/${line}_1.fastq \
    -2 $run_dir/${line}_2.fastq \
    -p 20 \
    -m trimgalore,TEMP,ngs_te_mapper,retroseq,popoolationte,popoolationte2,te-locate,relocate2
