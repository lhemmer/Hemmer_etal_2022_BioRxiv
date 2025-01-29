#!/bin/bash

#SBATCH -J cat
#SBATCH -e cat.errL%j
#SBATCH -o cat.outL%j
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH -c 1
#SBATCH -n 1
#SBATCH -p preempt

module load python

#set working directory
workingDir="/scratch/Lucas/work_dir"
#set subprograms in the McClintock Pipeline
subprograms=("ngs" "temp" "popoolationte" "telocate" "retroseq")
#set populations if you need to
#populations=("59" "91")


## When you have a list of files
populations=()
while IFS= read -r line; do
 populations+=("$line")
done < /scratch/Lucas/work_dir/mcc.list.txt

echo "populations"

echo "Working in dir: $workingDir"


#loop through populations
for population in ${populations[@]}; do
	echo "population: ${population}"
	#loop through McClintock components
	for subProgram in ${subprograms[@]}; do
		echo "subProgram: ${subProgram}"

		#popFile="${population}sequence_${subProgram}_nonredundant.bed"
		popFile="${population}_${subProgram}_nonredundant.bed"
		echo "Population file: $popFile"

		#categorize all subprograms' results
		python /scratch/Lucas/work_dir/categorize_mcc_output.py ${workingDir}/${popFile} /scratch/Lucas/work_dir/${population}.${subProgram}.ref.txt /scratch/Lucas/work_dir/${population}.${subProgram}.nonref.txt /scratch/Lucas/work_dir/TableS11.dmel_scaffold2_plus0310_context_without.csv

	done
done

