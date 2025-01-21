The following repository contains datasets for analysis and generating figures in:

Centromere-associated retroelement evolution in *Drosophila melanogaster* reveals an underlying conflict. BioRxiv [https://doi.org/10.1101/2022.11.25.518008](https://doi.org/10.1101/2022.11.25.518008).
Lucas Hemmer<sup>1,2</sup>, Sherif Negm<sup>1,3</sup>, Xuewen Geng<sup>1</sup>, Cecile Courret<sup>1</sup>, Beatriz Navarro-Domínguez<sup>1,4</sup>, Iain Speece<sup>1</sup>, Eddyson Altidor<sup>1</sup>, James Chaffer<sup>1</sup>, John S. Sproul<sup>1,5</sup>, Amanda M. Larracuente<sup>1</sup>

<sup>1</sup> Department of Biology, University of Rochester, Rochester, NY, USA
<sup>2</sup> Invaio Sciences, Cambridge, MA, USA
<sup>3</sup> Department of Human Genetics, University of Chicago, Chicago, IL, USA
<sup>4</sup> Departamento de Genética, Facultad de Ciencias, Universidad de Granada, Granada, 18071, Spain
<sup>5</sup> Department of Biology, University of Nebraska Omaha, Omaha, NY, USA


# **ABSTRACT**


Centromeres are chromosomal regions essential for coordinating chromosome segregation during cell division. While centromeres are defined by the presence of a centromere-specific histone H3 variant called CENP-A rather than a particular DNA sequence, they are typically embedded in repeat-dense and heterochromatic chromosomal genome regions. In many species, centromeres are associated with transposable elements but it is unclear if these elements are selfish and target centromeres for insertion or if they play a role in centromere specification and function. Here we use Drosophila melanogaster as a model to understand the evolution of centromere-associated transposable elements. G2/Jockey-3 is a non-LTR retroelement in the Jockey clade and the only sequence shared by all centromeres. We study the evolution of G2/Jockey-3 using short and long read population genomic data to infer insertion polymorphisms across the genome. We combine estimates of the age, frequency, and location of insertions to infer the evolutionary processes shaping G2/Jockey-3 and its association with the centromeres. We find that G2/Jockey-3 is an active retroelement that is targeted by the piRNA pathway. Our results suggest that G2-Jockey-3 is highly enriched in centromeres at least in part due to an insertion bias. We do not detect any signature of positive selection on any G2/Jockey-3 insertions that would suggest than individual insertions are favored by natural selection. Instead, we infer that most insertions are neutral or weakly deleterious both inside and outside of the centromeres. Therefore, G2/Jockey-3 evolution is consistent with it being a selfish genetic element that targets centromeres. We suspect targeting centromeres for insertion helps active retroelements like G2/Jockey-3 escape host defenses, as the unique centromeric chromatin may prevent targeting by the host silencing machinery. On the other hand, centromeric TEs insertions may be tolerated or even beneficial if they also contribute to the right transcriptional and chromatin environment. Thus, we suspect centromere-associated retroelements like G2/Jockey-3 reflect a balance between conflict and cooperation at the centromeres.


Please direct questions to Lucas Hemmer at lhemmer08@gmail.com or Amanda Larracuente at alarracu@bio.rochester.edu.

The scripts used to analyze these data are found at [https://github.com/LarracuenteLab/Dmel_Jockey-3_Evolution/](https://github.com/LarracuenteLab/Dmel_Jockey-3_Evolution/)

We organized the data files and scripts in this repository according to the figures or supplementary files in Hemmer et al. 2023. 


## Datasets

The datasets are organized and used for the following:

### Alignments

These are primarily used for the Age of Allele Test. Those elements used in the Age-of-Allele test were detected by the McClintock pipeline and Figures 6B and S6

The files include:

alignment_dmel_scaffold2_plus0310_2.fasta.out.G2.fasta -- alignment of all *G2/Jockey-3 elements*

alignment_mel_BS_agest_DGRP3.fasta -- alignment of *BS* copies for age-of-allele test

alignment_mel_DOC_agest_DGRP3.fasta -- alignment of *Do*c copies for age-of-allele test

alignment_mel_DOC2_DM_agest_DGRP3.fasta -- alignment of *Doc2* copies for age-of-allele test

alignment_mel_G5_DM_agest_DGRP3.fasta -- alignment of *G5* copies for age-of-allele test

alignment_mel_jockey-3_agetest_DGRP3.fasta -- alignment of *G2/Jockey-3* copies for age-of-allele test


### Diversity

These are the outputs of analyzing the VCF files using SimplifyVCF_Basic.v2.pl and Windows_Basic.v2.pl Perl scripts from [https://github.com/bnavarrodominguez/SD-Population-genomics](https://github.com/bnavarrodominguez/SD-Population-genomics). Genetic diversity (pi) and Tajima's D were calculated in 500 bp windows along the entire length of each contig. These data were used to generate Figure S7 and S8.

The files include:DGRP3.all.2L_1.500.csvDGRP3.all.2R_21.500.csvDGRP3.all.3L_1.500.csvDGRP3.all.3R_28.500.csvDGRP3.all.other.500.csvDGRP3.all.X_1.500.csv

### GATK

These filtered SNP vcf files were generated by the GATK pipeline using best practices. They are divided by contig. The centromeric contigs have two vcf files, one with SNPs detected in the satellite DNA regions and one without SNPs from those regions. These were used to generate Figures S9-S11.

The files include:

All_filtered_snps.2L_1.vcf.gz

All_filtered_snps.2R_21.vcf.gz

All_filtered_snps.3L_1.vcf.gz

All_filtered_snps.3R_28.vcf.gz

All_filtered_snps.variant.select.3R_5.2.nosat.vcf

All_filtered_snps.variant.select.3R_5.2.vcf

All_filtered_snps.variant.select.Contig79.2.vcf

All_filtered_snps.variant.select.Contig79.nosat.vcf

All_filtered_snps.variant.select.Contig119.2.vcf

All_filtered_snps.variant.select.Contig119.nosat.vcf

All_filtered_snps.variant.select.tig00057289.2.nosat.vcf

All_filtered_snps.variant.select.tig00057289.2.vcf

All_filtered_snps.variant.select.vcf


### McClintock

These files were used with the McClintock pipeline found at [https://github.com/bergmanlab/mcclintock](https://github.com/bergmanlab/mcclintock). They include the following:

dmel_scaffold2_plus0310_2.fasta -- genome FASTA file

dmel_scaffold2_plus0310_2.fasta.out.gff -- genome transposon annotation file from Repeatmasker

dmel.chromosomes.fa.TE.mcClintock.2020.gff -- modified transposon annotation file for running McClintock

dmel.chromosomes.fa.TE.mcClintock.2020.tsv -- modified transposon annotation list for running McClintock

specieslib_mcClintock_2020_mod.fasta -- the transposon library


### Paper

This are files directly mentioned in the paper. They include the following:

File_S4_McClintock_output_DGRP.txt
File_S5_DSPR_calls_wFreq.txt
File_S10_specieslib_McClintock.fasta

Note that File_S10 is the same as specieslib_mcClintock_2020_mod.fasta under the McClintock directory.

