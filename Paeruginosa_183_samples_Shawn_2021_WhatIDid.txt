Location -
/storage/data/DATA4/analysis/27_Paeurginosa_183_samples_Shawn

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

cd /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/1_Rawdata/Samples_without_Topup

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
mv */*bbmap_adaptertrimmed.fq /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/2_AdapterTrimmed_bbduk_Q30/2_AdapterTrimmed_bbduk_Q30

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

Login: ramadatta.88 - Krisarc8892

Upload the metadataform and fasta files (gteq1kb size) - Click Submit! Provide email address to get intimated by CGE server when the data processing is done.

Awaiting results!


6. Kraken2 for species calling on adapter trimmed reads (prior to this project - we used only kraken1):
########################################################

cd 2_AdapterTrimmed_bbduk/

$ time for d in $(ls -d */ | tr -d '/'); do echo $d; cd $d; R1=`ls *_1_*.fq | sed 's/.fq//g'`; R2=`ls *_2_*.fq | sed 's/.fq//g'`; echo "$R1 $R2"; kraken2 --db /storage/apps/Kraken2/Kraken2_db/minikraken_8GB_20200312 --threads 48 --report "$d"_kraken2_report --paired "$R1".fq "$R2".fq --output "$d"_kraken2_result; cd ../; done

real	59m52.226s
user	381m15.221s
sys	55m7.728s


--> Top species

$ for d in $(ls */*report); do echo $d; egrep -m1 "\sS\s" $d; done >>ALL_kraken2_report_top1species 
$ for d in $(ls */*report); do echo $d; egrep -m2 "\sS\s" $d; done >>ALL_kraken2_report_top2species

mkdir ../7_KRAKEN2
mv *kraken2* ../7_KRAKEN2
mv */*kraken2* ../7_KRAKEN2

7. Kleborate (Only if all species are Klesiella Species)
#########################################################

 Backburner...

 mkdir 8_Kleborate

8. Prokka
#########

MLST tool assignment: 

- 179 samples/183 with P.aeruginosa ST308 (same as jeanette's paper)

Out of the rest 4 samples, 

- 3 samples are assigned with P.aeruginosa, No ST assigned, 
- 1 sample with No Species, No ST

So using 179 samples for prokka because of same species and ST

$ mkdir 9_Prokka
$ cd 9_Prokka

$ ln -s /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/4_SPAdes_Assemblies/gteq1kb/*.fasta .
$ for d in $(cat log_mlst | sort -nrk3,3 | fgrep '-' | cut -f1); do echo $d; rm ../9_Prokka/$d ; done # remove No Species, No ST samples

$ for d in $(ls *.fasta | sed 's/_spades.gte1kb.contigs.fasta//g'); do echo "$d"; prokka --force --cpus 32 --outdir "$d"_prokka_out --prefix "$d" "$d"_spades.gte1kb.contigs.fasta >>prokka_log 2>>prokka_error; done

mkdir ../10_Roary

cp */*.gff ../10_Roary/

9. Roary
########

$ time roary -e --mafft -g 60000 -p 48 -cd 95 *.gff

real	41m2.056s
user	591m33.247s
sys	10m39.796s


$ mkdir ../11_SNP_Validation

--> Create a list of all the samples according to ST type to be SNP validated

# This oneliner generates list of files with the sample IDs according to the ST assigned by MLST tool

# Make sure there is a tab delimted space when using the fgrep

for d in $(cat ../5_MLST/log_mlst | awk '{print $3}' | sort -u);do 
num=`grep -P "\t$d\t" ../5_MLST/log_mlst | wc -l`; 
grep -P "\t$d\t" ../5_MLST/log_mlst | awk '{print $1 "\t" $3}' | sed 's/_spades.*//g' >list_PAE_ST"$d"_"$num"samples; 
done

$ cp list_PAE_ST*samples../11_SNP_Validation/
$ cp core_gene_alignment.aln ../11_SNP_Validation

10. SNP Validation
###################

cd ../11_SNP_Validation  

cp /storage/apps/SNP_Validation_Scripts/tools/snpvalidation.sh .

In the snpvalidation.sh script, change the path of the directory for R1 and R2 to run snippy. For example, for wenjian project here I have put:

/storage/data/DATA4/analysis/24_Wenjian_data_analysis/2_AdapterTrimmed_bbduk_Q30/$d/"$d"*_1_bbmap_adaptertrimmed.fastq #---> for R1 
/storage/data/DATA4/analysis/24_Wenjian_data_analysis/2_AdapterTrimmed_bbduk_Q30/$d/"$d"*_2_bbmap_adaptertrimmed.fastq #---> for R2 


$ time ./snpvalidation.sh list_PAE_ST-_4samples




time ./snpvalidation.sh list_PAE_ST308_179samples
real	265m40.723s
user	2240m11.863s
sys	143m41.217s



11. Gubbins (We no more do padding before gubbins) &&  Calculate SNP difference
######################################################################

Once the snpvalidation is done, change directory into the PAE_ST308_179samples, and run gubbins

cd PAE_ST308_179samples # Later repeat the following steps to the WJ_Klebsiella_ST23_4samples also

$ mkdir Gubbins

$ cp core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta Gubbins/

$ cd Gubbins/

# Convert the consensus alignment into one liner fasta format

$ perl /storage/apps/SNP_Validation_Scripts/tools/convertFastatoOneLine.pl core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta

# Assuming consensus is on top remove consensus sequence from the above fasta file
sed '1,2d' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta_OneLine.fasta >>core_gene_alignment_woConsensus_editedByGATK_renamed.fasta

$ time run_gubbins.py --prefix postGubbins --filter_percentage 100 --threads 48 core_gene_alignment_woConsensus_editedByGATK_renamed.fasta --verbose 
real	9m17.957s
user	69m1.245s
sys	0m23.285s
$ cp /storage/apps/SNP_Validation_Scripts/tools/dnaDist.r .

$ Rscript dnaDist.r 

(Repeat the above steps in 11.Gubbins for all the ST groups - for this project the second ST group is : WJ_Klebsiella_ST23_4samples)


12. Clustering (Did not do this for Wenjian - but documenting for reference)
############################################################################ 

# We want to cluster all the isolates which have a SNP difference  with certain threshold (for this case : 0 SNPs, Therefore my awk has $NF==0)

# The following one-liner takes the SNP difference file | takes pair1,pair2,SNP Diff | removes double quotes if any | extract pairs satisfyingSNP Diff threshold (Here 0 SNPs) | print rhe pair1, pair2 to a file

$ cat postGubbins.filtered_polymorphic_sites_melt_sorted.csv | awk -F',' '{print $2"\t"$3"\t"$4}' |  tr -d "\"" | awk '$NF==0' | awk '{print $1"\t"$2}' >pairsFile.list

# The Clustering Script will recruit the isolates if any of the isoaltes has link with another to the group

$ cp /storage/apps/PlasmidSeeker/Databases_aftr_04122019/ECCMID_2020/Deduplicated_Fasta_for_PlasmidSeeker_DB/ClustrPairs_v3.pl .

$ perl ClustrPairs_v3.pl pairsFile.list 


13. Visualization of tree
#########################

Using Figtree: 

$ java -jar /storage/apps/FigTree_v1.4.3/lib/figtree.jar

File -> Open -> postGubbins.final_tree.tre



To generate plots
#################

1) Pre and Post Trimmed Reads

qualStats.R (on datta's desktop)
qualStats_aggregated.R (on datta's desktop)


2) Assembly Size

Full assemblies & gteq1kb:

$ statswrapper.sh *.fasta >full_assemblies_summary_stats.txt
$ statswrapper.sh *.fasta >gteq1kb_assemblies_summary_stats.txt

Extract columns (n_contigs_full_assembly, contig_bp_full_assembly, ctg_N50_full_assembly,ctg_max_full_assembly) from the above files and create one file (full_gteq1kb_combined.csv)

Run in R console

plotAssembly_stats.R

3) MLST

Run in R console

plotMLST.R

4) CGE - (https://cge.cbs.dtu.dk/services/cge/index.php)

Run in R console

plotCGE.R

5) Kraken

6) Table with SNP results

plotSNPDiff_matrix.R --> Need full triangle to plot the heatmap 

7) Phylogenetic tree


On March 23rd, Dr. Ng wanted to add 31 samples of PAE from Jeanette's study from NUH to our data and run the SNIPPY using PA01 as reference genomes and form a BEAST evolutionary tree.
So, downloaded the 31 samples data from NCBI link : https://www.ncbi.nlm.nih.gov/bioproject/PRJNA507901

# Download the sequences using following weblink by providing SRX number concatenated with comma
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name

or use fasterq-dump

for d in $(cat SraAccList.txt; do echo $d; fasterq-dump $d -p; done

# Rename SRR Names to Sample Names in paper
$ while read line; do echo $line; srr=`echo $line | awk '{print $1}'`; sample=`echo $line | awk '{print $2}'`; echo "srr: $srr, sample: $sample"; mv "$srr"_1.fastq "$sample"_1.fastq; mv "$srr"_2.fastq "$sample"_2.fastq; done <SRRIds_SamplesNames.list 


# Make sure the file extension is as needed (fastq,fq,fastq etc)

for d in $(ls *.fastq | cut -f1 -d "_"  | sort -u); do mkdir $d; mv $d*.fastq $d; done

# Run fastqc on the raw reads
time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 3 *.fastq & cd ..; done

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
$ time for d in $(ls -d */); do echo $d; subdir=`echo $d`; cd $subdir; R1=`ls *_1.fastq | sed 's/.fastq//g'`; R2=`ls *_2.fastq | sed 's/.fastq//g'`; echo "$R1 $R2"; /storage/apps/bbmap/bbduk.sh -Xmx6g in1=$R1.fastq in2=$R2.fastq out1=$R1\_bbmap_adaptertrimmed.fastq out2=$R2\_bbmap_adaptertrimmed.fastq ref=/storage/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=30 minavgquality=30; cd ..; done

manually copied the bbmap command log into adapter_trimming.log

mkdir ../2_AdapterTrimmed_bbduk_Q30

# move the adapter trimmed files into another directory
mv */*bbmap_adaptertrimmed.fq /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/2_AdapterTrimmed_bbduk_Q30/2_AdapterTrimmed_bbduk_Q30

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

## Running SNIPPY on 179 ST308 samples internal and 31 Jeanette's samples

$ cd /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/

$ mkdir 15_SNIPPY_internal179_Jeanette_31_ST308samples

#Generate a list of ID R1 R2 for running multiple samples
##########################################################

$ for d in $(cat list_PAE_ST308_179samples); do echo "$d"; find /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/ -name "$d*.fq*" | grep 'bbmap' | sort; done | paste - - - >input.tab

2_AdapterTrimmed_bbduk_Q30]$ for d in $(ls -d */ | tr -d "/"); do echo "$d"; find /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30 -name "$d*bbmap*.fastq" | sort; done | paste - - - >>/storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/15_SNIPPY_internal179_Jeanette_31_ST308samples/input.tab

$ /storage/apps/snippy/bin/snippy-multi input.tab --ref PAO1_reference.fasta --cpus 32 > paeruginosa_snippy_runme.sh

$ ./paeruginosa_snippy_runme.sh # This script generates all the core genome alignments

# Remove all the "weird" characters and replace them with N. This is required if we are running gubbins

$ /storage/apps/snippy/bin/snippy-clean_full_aln core.full.aln > clean.core.full.aln

# Run gubbins on clean.core.full.aln

$ run_gubbins.py -p gubbins clean.core.full.aln --threads 48 --verbose

# Generate a cleancore.full.recomb.masked.fasta alignment

$ maskrc-svg.py --aln clean.core.full.aln --out cleancore.full.recomb.masked.fasta --gubbins gubbins

# Pad Invariant Sites after Gubbins

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


##SRST2 check for ICETn43716385 

SRST2 has some version issue. So, created a new SRST2 environment and installed the SRST2

conda activate srst2_environment

Extracted region :

trbl        3401482 - 3405853
ndm        3431913
radc        3463390
integrase     3473759 - 3474958



for d in $(ls /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/2_AdapterTrimmed_bbduk_Q30/*/*_1_bbmap_adaptertrimmed.fq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/17_srst2/"$base"_S1_L001_R1_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/2_AdapterTrimmed_bbduk_Q30/*/*_2_bbmap_adaptertrimmed.fq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/17_srst2/"$base"_S1_L001_R2_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30/*/*_1_bbmap_adaptertrimmed.fastq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/17_srst2/"$base"_S1_L001_R1_001.fq; done

for d in $(ls /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/14_Jeanette_31samples_NUH/2_AdapterTrimmed_bbduk_Q30/*/*_2_bbmap_adaptertrimmed.fastq); do echo "$d"; base=`echo "$d" | awk -F/ '{print $(NF-1)}'`;ln -s $d /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/17_srst2/"$base"_S1_L001_R2_001.fq; done

$ for d in $(ls *.fq | awk -F_ '{print $1}'|sort -u); do echo $d; mkdir $d; mv "$d"_* $d; done

$ time for d in $(ls -d */ | tr -d "/"); do echo $d; cd $d;  R1=`ls *_S1_L001_R1_001.fq`; R2=`ls *_S1_L001_R2_001.fq`; echo "R1: $R1----R2: $R2"; time srst2 --input_pe $R1 $R2 --output test_srst2 --log --gene_db /storage/data/DATA4/analysis/27_Paeruginosa_183_samples_Shawn/16_Check_presence_of_ICETn43716385/ICETn43716385.fasta --threads 48; cd ..; done

$ cat */test_srst2__fullgenes__ICETn43716385__results.txt | awk '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $7}' | grep -v 'coverage' | column -t


# To stitch contigs based on the mapped regions to ICETn43716385

time bash DenovoAssemble_fromBAM_AND_orderGenomeContigs.sh &

# This C01992 samples was giving a problem so removing to run BEAST

sed  '/C01992|2020.39071038251/,+1d' cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates.fasta >cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend_withDates.fasta_BEAST_withDates_ExcludedFastawithoutDates_Exlcuded_C01992_sample.fasta
