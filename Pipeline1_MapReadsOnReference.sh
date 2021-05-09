Location -
/storage/data/DATA4/analysis/27_Project

#################################################################################################################################

1.fastqc on Raw Reads 
######################

mkdir 1_Rawdata

cd 1_Rawdata

# Move the files to respective folders
awk -F/ '{$NF=""; print $0}' Paeruginosa_183_samples | tr ' ' '/' >Paeruginosa_183_samples_locations
for d in $(cat Paeruginosa_183_samples_locations); do ls "$d"*.fq*; ln -s "$d"*.fq.gz .; done #Never use force (-sf) when creating softlinks using ln

# Make sure the file extension is as needed (fastq,fq,fastq etc)

for d in $(ls *.gz | cut -f1 -d "_"  | sort -u); do mkdir $d; mv $d*.gz $d; done

# Run fastqc on the raw reads
time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 48 *.fq.gz; cd ..; done

# Samples with Topup sequencing are copied seperately into the folder "Samples_with_Topup" and the rest of the samples are in Samples_without_Topup

$ mkdir Samples_with_Topup
for d in $(ls */*.fq.gz | awk -F/ '{print $1}' | sort | uniq -c | sort -n | grep "^\s*4" | awk '{print $2}'); do echo $d; mv $d Samples_with_Topup/; done

#Combine the topup sequencing files

$ for d in $(ls -d */ | tr -d "/"); do cd $d; echo $d; base=`ls *DDM*HFKWVCCX2*L8_1.fq.gz | awk '{print $1}' | cut -f2 -d "_"`; cat *_L8_1.fq.gz >"$d"_"$base"_combined_topup_L8_1.fq.gz; cd ..; done
$ for d in $(ls -d */ | tr -d "/"); do cd $d; echo $d; base=`ls *DDM*HFKWVCCX2*L8_2.fq.gz | awk '{print $1}' | cut -f2 -d "_"`; cat *_L8_2.fq.gz >"$d"_"$base"_combined_topup_L8_2.fq.gz; cd ..; done

# Unzip the gz files 
time for d in $(ls -d */ | tr -d "/"); do echo $d; cd $d; gunzip *.gz; cd ..; done

cd /storage/data/DATA4/analysis/27_Project/1_Rawdata/Samples_without_Topup

# Unzip the gz files - because these are softlinks we write to  an other file
time for d in $(ls *.fq.gz); do base=`echo $d | sed 's/.fq.gz//'`; gunzip -c $d >"$base".fq; done

#combine all fastqc reports using multiqc
multiqc .

Observations: adapters are present in many sequences. So, trimmed them using BBMAP.

But base quality of all the samples looking high (above Q30).

Q30 Stats
#########

The perl script calculates the Total_Reads, Total_Bases, Total_Q20_Bases, Total_Q30_Bases, Mean Read Length from the fastq file

time for d in $(ls */*fq); do echo $d; perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d & done 

mv Q30_Q20_readstats.txt rawData_Q30_Q20_readstats.txt

Let us plot the numbers generated from the above perl script

mkdir rawDataStats_RPlots
mv rawData_Q30_Q20_readstats.txt rawDataStats_RPlots
cd rawDataStats_RPlots

# Change variables in the Rscript beginning portion and run the R script from terminal: ---->had some issues in making R run on server so running Rscript steps on my desktop

# The Rscript generates the figures and the results of the above quality stats can be found in file:
 /data02/Analysis/for_Colloborators/for_Wenjian/Feb2021/Wenjian_12_Samples_DataAnalysis/rawDataStats_RPlots/Wenjian_samples_Sequencing Company_Raw data_aggregated_stats.txt


2. Adapter trimming BBDuk - Done!
####################################

#Samples_with_Topup
$ time for d in $(ls -d */); do echo $d; subdir=`echo $d`; cd $subdir; R1=`ls *8_1.fq | sed 's/.fq//g'`; R2=`ls *8_2.fq | sed 's/.fq//g'`; echo "$R1 $R2"; /storage/apps/bbmap/bbduk.sh -Xmx6g in1=$R1.fq in2=$R2.fq out1=$R1\_bbmap_adaptertrimmed.fq out2=$R2\_bbmap_adaptertrimmed.fq ref=/storage/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=30 minavgquality=30; cd ..; done

#Samples_without_Topup
time for d in $(ls -d */); do echo $d; subdir=`echo $d`; cd $subdir; R1=`ls *8_1.fq | sed 's/.fq//g'`; R2=`ls *8_2.fq | sed 's/.fq//g'`; echo "$R1 $R2"; /storage/apps/bbmap/bbduk.sh -Xmx6g in1=$R1.fq in2=$R2.fq out1=$R1\_bbmap_adaptertrimmed.fq out2=$R2\_bbmap_adaptertrimmed.fq ref=/storage/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=30 minavgquality=30; cd ..; done

manually copied the bbmap command log into adapter_trimming.log

mkdir ../2_AdapterTrimmed_bbduk_Q30

# move the adapter trimmed files into another directory
mv */*bbmap_adaptertrimmed.fq /storage/data/DATA4/analysis/27_Project/2_AdapterTrimmed_bbduk_Q30/2_AdapterTrimmed_bbduk_Q30

# move to specific folder
$ for d in $(ls *.fq| awk -F_ '{print $1}' | sort -u); do echo $d; mkdir $d; mv "$d"_* $d; done

Q30 Stats
#########

The perl script calculates the Total_Reads, Total_Bases, Total_Q20_Bases, Total_Q30_Bases, Mean Read Length from the fastq file

time for d in $(ls */*fq); do echo $d; perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d & done 

mv Q30_Q20_readstats.txt PostTrimmed_Q30_Q20_readstats.txt

Let us plot the numbers generated from the above perl script

mkdir PostTrimmedDataStats_RPlots
mv PostTrimmed_Q30_Q20_readstats.txt PostTrimmedDataStats_RPlots
cd PostTrimmedDataStats_RPlots

# Change variables in the Rscript beginning portion and run the R script from terminal: ---->had some issues in making R run on server so running Rscript steps on my desktop

# The Rscript generates the figures and the results of the above quality stats can be found in file:
 /data02/Analysis/for_Colloborators/for_Wenjian/Feb2021/Wenjian_12_Samples_DataAnalysis/PostTrimmedDataStats_RPlots/Wenjian_samples_Sequencing Company_Post Trimming with Q30 score_aggregated_stats.txt


fastqc after cleaning adapters
#################################

time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 48 *.fq & cd ..; done

#combine all fastqc reports using multiqc
multiqc .

Observations: adapters are removed and the bases are quality trimmed in sequences. Plot looks much better now.

2. Spades
###########

time for d in $(ls -d */ | grep -v 'multiqc'); do echo $d; subdir=`echo -n $d | tr -d "/"`; cd $subdir; R1=`ls *_1_bbmap_adaptertrimmed.fq `; R2=`ls *_2_bbmap_adaptertrimmed.fq`; echo "$R1 $R2"; spades.py --pe1-1 "$R1" --pe1-2 "$R2" -o "$subdir"_spades --careful -t 48 >>log_spades; cd ..; done

[datta@ttshbio 2_AdapterTrimmed_bbduk]$ mkdir ../3_SPAdes
[datta@ttshbio 2_AdapterTrimmed_bbduk]$ mv */*_spades/ ../3_SPAdes/

#Copy and rename the assemblies
cd ../3_SPAdes/
$ for d in $(ls *_spades/contigs.fasta); do prefix=`echo "$d" | cut -f1 -d "/"`; cp "$d" "$prefix"_contigs.fasta; done

$ mkdir ../4_SPAdes_Assemblies
$ mv *contigs.fasta ../4_SPAdes_Assemblies/
$ cd ../4_SPAdes_Assemblies/

Filter for contigs >1kb
------------------------
used bioawk installed in server

$ for d in $(ls *.fasta | sed 's/_contigs.fasta//g' ); do echo $d; /storage/apps/bioawk/bioawk -c fastx 'length($seq) >=1000 {print "\>"$name"\n"$seq}' "$d"_contigs.fasta >$d.gte1kb.contigs.fasta; done
$ mkdir gteq1kb
$ mv *gte1kb* gteq1kb/
$ mkdir full_assembly
$ mv *.fasta full_assembly/

3. MLST for species calling (from here the analysis is done only on the fasta with contigsize gteq 1Kb)
####################################

$ cd gteq1kb
mlst *contigs.fasta | grep ".fasta" >>log_mlst

[datta@ttshbio gteq1kb]$ mkdir ../../5_MLST
[datta@ttshbio gteq1kb]$ mv log_mlst ../../5_MLST/


5. CGE - Bacterial analysis pipeline
####################################

$ mkdir ../6_CGE

Download the metadataform here: https://cge.cbs.dtu.dk/tools/client/testUploader/batch/php/metadataform.xlsx
(a template can be taken from /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_NovogeneAIT_2020/6_CGE)

Fill the form with the sample details 

Go to webpage: https://cge.cbs.dtu.dk/services/cge/

Login: ramadatta.88 - nb9GxPxm

Upload the metadataform and fasta files (gteq1kb size) - Click Submit! Provide email address to get intimated by CGE server when the data processing is done.

Awaiting results!


6. Kraken2 for species calling on adapter trimmed reads (prior to this project - we used only kraken1):
########################################################

cd 2_AdapterTrimmed_bbduk/

$ time for d in $(ls -d */ | tr -d '/'); do echo $d; cd $d; R1=`ls *_1_*.fq | sed 's/.fq//g'`; R2=`ls *_2_*.fq | sed 's/.fq//g'`; echo "$R1 $R2"; kraken2 --db /storage/apps/Kraken2/Kraken2_db/minikraken_8GB_20200312 --threads 48 --report "$d"_kraken2_report --paired "$R1".fq "$R2".fq --output "$d"_kraken2_result; cd ../; done

--> Top species

$ for d in $(ls */*report); do echo $d; egrep -m1 "\sS\s" $d; done >>ALL_kraken2_report_top1species 
$ for d in $(ls */*report); do echo $d; egrep -m2 "\sS\s" $d; done >>ALL_kraken2_report_top2species

mkdir ../7_KRAKEN2
mv *kraken2* ../7_KRAKEN2
mv */*kraken2* ../7_KRAKEN2

7. Kleborate (Only if all species are Klesiella Species)
#########################################################

#Not performed for this project since these are PAE samples

 mkdir 8_Kleborate

8. Core Genome alignment and variant calling
#############################################

# From MLST tool assignment, 179 samples/183 with P.aeruginosa ST308 (NUH samples from Jeanette's paper are also ST308)

# On March 23rd, Dr. Ng wanted to add 31 samples of PAE from Jeanette's study from NUH to our data and run the SNIPPY using PA01 as reference genomes
# and generate a BEAST evolutionary tree. So, downloaded the 31 samples data from NCBI link : https://www.ncbi.nlm.nih.gov/bioproject/PRJNA507901

#The details and raw-read processing of these 31 samples are listed in the bottom (please refer Ad-Hoc Analysis 1). The adapter-trimmed-quality-filtered data is used for the core-genome aligned

## Running SNIPPY on 179 ST308 samples internal and 31 Jeanette's samples to generate recombinant free polymorphic SNP sites

$ cd /storage/data/DATA4/analysis/27_Project/

$ mkdir 15_SNIPPY_internal179_Jeanette_31_ST308samples

#Generate a list of ID R1 R2 for running multiple samples for SNIPPY input (based on SNIPPY's documentation)
############################################################################################################

# 179 PAE samples
$ for d in $(cat list_PAE_ST308_179samples); do echo "$d"; find /storage/data/DATA4/analysis/27_Project/ -name "$d*.fq*" | grep 'bbmap' | sort; done | paste - - - >input.tab

# 31 NUH samples
2_AdapterTrimmed_bbduk_Q30]$ for d in $(ls -d */ | tr -d "/"); do echo "$d"; find /storage/data/DATA4/analysis/27_Project/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30 -name "$d*bbmap*.fastq" | sort; done | paste - - - >>/storage/data/DATA4/analysis/27_Project/15_SNIPPY_internal179_Jeanette_31_ST308samples/input.tab

# Create a runMe file of SNIPPY commands using the following command (based on SNIPPY's documentation)
$ /storage/apps/snippy/bin/snippy-multi input.tab --ref PAO1_reference.fasta --cpus 32 > paeruginosa_snippy_runme.sh

# Run SNIPPY pipeline
$ ./paeruginosa_snippy_runme.sh # This script generates all the core genome alignments

# Remove all the "weird" characters and replace them with N. This is required if we are running gubbins (based on SNIPPY's documentation)
$ /storage/apps/snippy/bin/snippy-clean_full_aln core.full.aln > clean.core.full.aln

# "clean.core.full.aln" file contains full genome alignment with SNPs flaged by SNIPPY. 
# Run gubbins on clean.core.full.aln to fiter the recombinant regions. Gubbins generates recombinant free polymorphic SNP sites.
$ run_gubbins.py -p gubbins clean.core.full.aln --threads 48 --verbose

#############--------------Divergence tree estimation - For BEAST analysis------------------#################

# The recombinant free polymorphic SNP sites fasta file can be used an input by 
#   1) by adding invariant sites (beast constant sites)  into XML (OR)
#   2) by adding invariant sites to sequences and sequences loaded into BEAST through XML.

# From experience, method 1 is magnitude times faster than method 2 simply because loading the whole genome into memory for computation is intensive than loading only SNP sites with constant sites.

# I am listing steps using both methods

#----------------> Divergence tree estimation - method 1

# Append date to the fasta headers

time perl /storage/apps/SNP_Validation_Scripts/tools/BEAST_DatesHeader.pl <inputfasta> <outputfile_Prefix>
time perl /storage/apps/SNP_Validation_Scripts/tools/BEAST_DatesHeader.pl postGubbins.filtered_polymorphic_sites.fasta postGubbins.filtered_polymorphic_sites # This command generates postGubbins.filtered_polymorphic_sites_withDates_forBEAST.fasta file

#----------------> Divergence tree estimation - method 2

# Steps: Mask recombination in the full alignment using maskrc-svg, 

# Generate a cleancore.full.recomb.masked.fasta alignment
$ maskrc-svg.py --aln clean.core.full.aln --out cleancore.full.recomb.masked.fasta --gubbins gubbins

# Pad Invariant Sites after Gubbins - this script takes genome calculates invariant sites, Runs trimal and removes recombinant regions but retains "-" by converting to Ns. Beast takes the input with gaps (-)/Ns
time /storage/apps/SNP_Validation_Scripts/tools/padInvariantSites_afterGubbins.sh PAO1_reference.fasta

# remove extra inforamtion other than sample names
sed -i 's/ 6091903 bp//g'  cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend.fasta

# Append date to the fasta headers
time perl /storage/apps/SNP_Validation_Scripts/tools/BEAST_DatesHeader.pl cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend.fasta cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta

# Exclude the fasta sequences without dates
grep -A1 "|" cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates.fasta >cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates.fasta

sed -i 's/--//g' cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates.fasta
sed -i '/^$/d' cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates.fasta

5 samples did not have dates so excluded from the analysis. So BEAST run on 205 samples from ST308


# SRST2 check for ICETn43716385 

SRST2 has some version issue. So, created a new SRST2 environment and installed the SRST2

conda activate srst2_environment

Extracted region :

trbl        3401482 - 3405853
ndm        3431913
radc        3463390
integrase     3473759 - 3474958



for d in $(ls /storage/data/DATA4/analysis/27_Project/2_AdapterTrimmed_bbduk_Q30/*/*_1_bbmap_adaptertrimmed.fq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Project/17_srst2/"$base"_S1_L001_R1_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Project/2_AdapterTrimmed_bbduk_Q30/*/*_2_bbmap_adaptertrimmed.fq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Project/17_srst2/"$base"_S1_L001_R2_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Project/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30/*/*_1_bbmap_adaptertrimmed.fastq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Project/17_srst2/"$base"_S1_L001_R1_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Project/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30/*/*_2_bbmap_adaptertrimmed.fastq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Project/17_srst2/"$base"_S1_L001_R2_001.fq; done

$ for d in $(ls *.fq | awk -F_ '{print $1}'|sort -u); do echo $d; mkdir $d; mv "$d"_* $d; done

$ time for d in $(ls -d */ | tr -d "/"); do echo $d; cd $d;  R1=`ls *_S1_L001_R1_001.fq`; R2=`ls *_S1_L001_R2_001.fq`; echo "R1: $R1----R2: $R2"; time srst2 --input_pe $R1 $R2 --output test_srst2 --log --gene_db /storage/data/DATA4/analysis/27_Project/16_Check_presence_of_ICETn43716385/ICETn43716385.fasta --threads 48; cd ..; done

$ cat */test_srst2__fullgenes__ICETn43716385__results.txt | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $7}' | grep -v 'coverage' | column -t


# To stitch contigs based on the mapped regions to ICETn43716385

time bash DenovoAssemble_fromBAM_AND_orderGenomeContigs.sh &

# This C01992 samples was giving a problem so removing to run BEAST

sed  '/C01992|2020.39071038251/,+1d' cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates.fasta >cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates_Exlcuded_C01992_sample.fasta

###########################################------------------------Ad-Hoc Analysis 1----------------------################################################

# On March 23rd, Dr. Ng wanted to add 31 samples of PAE from Jeanette's study from NUH to our data and run the SNIPPY using PA01 as reference genomes
# and generate a BEAST evolutionary tree. So, downloaded the 31 samples data from NCBI link : https://www.ncbi.nlm.nih.gov/bioproject/PRJNA507901

# Download the sequences using following weblink by providing SRX number concatenated with comma
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name or use fasterq-dump

$ for d in $(cat SraAccList.txt); do echo $d; fasterq-dump $d -p; done

# Rename SRR Names to Sample Names in paper
$ while read line; do echo $line; srr=`echo $line | awk '{print $1}'`; sample=`echo $line | awk '{print $2}'`; echo "srr: $srr, sample: $sample"; mv "$srr"_1.fastq "$sample"_1.fastq; mv "$srr"_2.fastq "$sample"_2.fastq; done <SRRIds_SamplesNames.list 

# Make sure the file extension is as needed (fastq,fq,fastq etc)
$ for d in $(ls *.fastq | cut -f1 -d "_"  | sort -u); do mkdir $d; mv $d*.fastq $d; done

# Run fastqc on the raw reads
$ time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 3 *.fastq & cd ..; done

#combine all fastqc reports using multiqc
$ multiqc .

#Observations: adapters are present in many sequences. So, trimmed them using BBMAP. But base quality of all the samples looking high (above Q30).

# Q30 Stats
############

# The perl script calculates the Total_Reads, Total_Bases, Total_Q20_Bases, Total_Q30_Bases, Mean Read Length from the fastq file passed as an input

$ time for d in $(ls */*fq); do echo $d; perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d & done 
mv Q30_Q20_readstats.txt rawData_Q30_Q20_readstats.txt

Let us plot the numbers generated from the above perl script

mkdir rawDataStats_RPlots
mv rawData_Q30_Q20_readstats.txt rawDataStats_RPlots
cd rawDataStats_RPlots

# Change variables in the Rscript beginning portion and run the R script from terminal: ---->had some issues in making R run on server so running Rscript steps on my desktop

# The Rscript generates the figures and the results of the above quality stats can be found in file:
 /data02/Analysis/for_Colloborators/for_Wenjian/Feb2021/Wenjian_12_Samples_DataAnalysis/rawDataStats_RPlots/Wenjian_samples_Sequencing Company_Raw data_aggregated_stats.txt


2. Adapter trimming BBDuk - Done!
####################################

#Samples_with_Topup
$ time for d in $(ls -d */); do echo $d; subdir=`echo $d`; cd $subdir; R1=`ls *_1.fastq | sed 's/.fastq//g'`; R2=`ls *_2.fastq | sed 's/.fastq//g'`; echo "$R1 $R2"; /storage/apps/bbmap/bbduk.sh -Xmx6g in1=$R1.fastq in2=$R2.fastq out1=$R1\_bbmap_adaptertrimmed.fastq out2=$R2\_bbmap_adaptertrimmed.fastq ref=/storage/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=30 minavgquality=30; cd ..; done

manually copied the bbmap command log into adapter_trimming.log

mkdir ../2_AdapterTrimmed_bbduk_Q30

# move the adapter trimmed files into another directory
mv */*bbmap_adaptertrimmed.fq /storage/data/DATA4/analysis/27_Project/2_AdapterTrimmed_bbduk_Q30/2_AdapterTrimmed_bbduk_Q30

# move to specific folder
$ for d in $(ls *.fastq| awk -F_ '{print $1}' | sort -u); do echo $d; mkdir $d; mv "$d"_* $d; done

Q30 Stats
#########

The perl script calculates the Total_Reads, Total_Bases, Total_Q20_Bases, Total_Q30_Bases, Mean Read Length from the fastq file

time for d in $(ls */*fq); do echo $d; perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d & done 

mv Q30_Q20_readstats.txt PostTrimmed_Q30_Q20_readstats.txt

Let us plot the numbers generated from the above perl script

mkdir PostTrimmedDataStats_RPlots
mv PostTrimmed_Q30_Q20_readstats.txt PostTrimmedDataStats_RPlots
cd PostTrimmedDataStats_RPlots

# Change variables in the Rscript beginning portion and run the R script from terminal: ---->had some issues in making R run on server so running Rscript steps on my desktop

# The Rscript generates the figures and the results of the above quality stats can be found in file:
 /data02/Analysis/for_Colloborators/for_Wenjian/Feb2021/Wenjian_12_Samples_DataAnalysis/PostTrimmedDataStats_RPlots/Wenjian_samples_Sequencing Company_Post Trimming with Q30 score_aggregated_stats.txt


fastqc after cleaning adapters
#################################

time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 3 *.fastq & cd ..; done

#combine all fastqc reports using multiqc
multiqc .

Observations: adapters are removed and the bases are quality trimmed in sequences. Plot looks much better now.


###########################################------------------------Ad-Hoc Analysis 1----------------------################################################



