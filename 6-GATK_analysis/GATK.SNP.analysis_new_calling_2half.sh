############### 	Variant analysis ##############
# echo "begin variant analysis..."
# 
# while IFS='' read -r line || [[ -n "$line" ]]; do
# 	echo "creating all variants list, added: $line"
# 	echo "$line.all_variants.vcf ">>$PREFIX.vcf.list
# 
# done < "$INFILE"
# 
cpu=$SLURM_NTASKS_PER_NODE
let mem=$SLURM_MEM_PER_NODE/1024
function GATK {
        module load java
	module load R
	/scratch/alarracu_lab/GATK_2019/gatk-4.1.2.0/gatk --java-options "-Xmx"$mem"g -Xms"$mem"g" $@
	module unload java
	module unload R
#    module unload gatk/3.5
}
REFERENCE="dmel_scaffold2_plus0310_2.fasta"
# 
# ############### 	Variant analysis ##############
# echo "begin variant analysis..."
# 
# GATK GenomicsDBImport \
#        --genomicsdb-workspace-path my_database \
#        --sample-name-map /scratch/alarracu_lab/centromere_population_genomics/gatk_run_dir/GATK_pipeline_update/variants/DGRP.test.vcf.list \
#        --reader-threads 24 \
#        -L intervals.list
#        #-L 2L -L 2R -L X -L 3L -L 3R -L X -L 4
# ##Merge all
# echo "merge vcf files..."
# 
# GATK  GenotypeGVCFs \
# -R $REFERENCE \
# -V gendb://my_database \
# -O All_variants.vcf.gz 
# ##-all-sites    
##--sample-ploidy 1

if [ -f "All_variants.vcf.gz" ]; then  
GATK VariantRecalibrator \
   -R $REFERENCE \
   -V All_variants.vcf.gz \
   --resource:DGRP,known=false,training=true,truth=true,prior=15.0 DGRP_All_SNPs.var_only.multiple_filters.All_filtered_snps.dip.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode BOTH \
   -O output.recal \
   --tranches-file output.tranches \
   --rscript-file output.plots.R
GATK ApplyVQSR \
   -R $REFERENCE \
   -V All_variants.vcf.gz \
   -O All_recal.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file output.tranches \
   --recal-file output.recal \
   -mode BOTH

echo "exclude indels from vcf files..."
GATK  SelectVariants \
-V All_variants.vcf.gz \
--select-type-to-include SNP \
-R $REFERENCE \
-O All_SNPs.vcf.gz
 
 else 
	echo "$INFILE.All_variants.vcf not created"
	exit 1
fi	

echo "final SNP filtering"
expression="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
#expression="\"AB < 0.2 || MQ0 > 50\""
module load java
/scratch/alarracu_lab/GATK_2019/gatk-4.1.2.0/gatk  VariantFiltration \
-R $REFERENCE \
-V All_SNPs.vcf.gz \
-O All_filtered_snps.vcf.gz \
--filter-name "filter" \
--filter-expression "$expression" 
 	
echo "Done"
