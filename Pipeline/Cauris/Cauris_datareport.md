---
title: "Candida auris Data Analysis "
author: "Prakki Sai Rama Sridatta"
date: "2021-10-12"
output:
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---



### Step 1: Raw read analysis fastqc and multiqc

Let us run `fastqc` and `multiqc` commands for analysing the quality of raw data.


```bash
# Location -
# cd /storage/data/DATA4/analysis/18_Cauris/2021/

###################################################################

# Step 1.fastqc on Raw Reads 
##############################

mkdir 1_Rawdata

cd 1_Rawdata

# Move the files to respective folders

ln -s /storage/data/DATA2/Sequencing_Company_Metadata/AITBiotech/Sep2021/September_batch2/IMH_C_auris_3_samples/X401SC21060194-Z01-F011/raw_data/*/*.fq.gz .
ln -s /storage/data/DATA2/Sequencing_Company_Metadata/AITBiotech/Sep2021/September_batch3/X401SC21060194-Z01-F016/raw_data/ADH0074/*.gz .

# Moving into fastq files into respective folders
$ time for d in $(ls *.gz | cut -f1 -d "_"  | sort -u); do mkdir $d; mv "$d"_*.gz $d; done

# A total of 139 samples have topup sequencing done
$ ls */*.fq.gz | awk -F/ '{print $1}' | sort | uniq -c | sort -n | grep "^\s*4" | wc -l
1

# Let's combine the topup sequencing files

$ time  for d in $(ls -d ADH0074/ | tr -d "/"); do cd $d; echo $d; base=`ls -1 *_1.fq.gz | awk '{print $1}' | cut -f2 -d "_" | sort -u`; echo $base; cat *_L*_1.fq.gz >"$d"_"$base"_combined_topup_L8_1.fq.gz; cd ..; done 

$ time  for d in $(ls -d ADH0074/ | tr -d "/"); do cd $d; echo $d; base=`ls -1 *_2.fq.gz | awk '{print $1}' | cut -f2 -d "_" | sort -u`; echo $base; cat *_L*_2.fq.gz >"$d"_"$base"_combined_topup_L8_2.fq.gz; cd ..; done 

# Unzip the gz files 
$ time for d in $(ls -d ADH0074/ | tr -d "/"); do echo $d; cd $d; gunzip *combined_topup*.gz; cd ..; done

# NOTE: because these are softlinks, we cannot just use gunzip. Rather, we write to an other file

$ time for d in $(ls -d ADH0075/ | tr -d "/"); do echo $d; cd $d; for e in $(ls *.fq.gz); do base=`echo $e | sed 's/.fq.gz//'`; gunzip -c $e >"$base".fq; done; cd ..; done

$ time for d in $(ls -d ADH0076/ | tr -d "/"); do echo $d; cd $d; for e in $(ls *.fq.gz); do base=`echo $e | sed 's/.fq.gz//'`; gunzip -c $e >"$base".fq; done; cd ..; done

# Run fastqc on the raw reads 

time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 48 *.fq; cd ..; done

# For generating read stats - Run this on both "Samples_without_Topup" and "Samples_with_Topup" folders
time for d in $(ls */*fq); do echo "perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d"; done >parallel_script_Q20_Q30_stats.sh
time parallel < parallel_script_Q20_Q30_stats.sh &

#combine all fastqc reports using multiqc
cd 1_Rawdata

multiqc .

# Observations from raw fastq: Adapters are present in many sequences. So, trimmed them using BBMAP.
# But base quality of all the samples looking high (above Q30).

$ cat */Q30_Q20_readstats.txt >rawData_Q30_Q20_readstats.txt

# Run qualStats.R script with the rawData_Q30_Q20_readstats.txt as an input. This should generate "_Pre-Trimmed_aggregated_stats.txt" file

$ pwd
# /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_Cleandata/illumina_bbmap_trimmed

# Run qualStats.R script with the rawData_Q30_Q20_readstats.txt as an input. This should generate "_Pre-Trimmed_aggregated_stats.txt" file

# Run qualStats.R script with the filtData_Q30_Q20_readstats.txt as an input. This should generate "_Post-Trimmed_aggregated_stats.txt" file

# Now combine both Raw and Filt Data using qualStats_aggregated.R script
```

### Step 2: Pre-process Illumina reads 

Now let us trim the raw data from illumina using `BBduk`.


```bash

$ cd /storage/data/DATA4/analysis/18_Cauris/2021/1_Rawdata

# Step 2. Adapter trimming using BBDuk 
######################################

time for d in $(ls -d */); 
do 
  echo $d; 
  subdir=`echo $d`; 
  cd $subdir; 
  R1=`ls *_1.fq | sed 's/.fq//g'`; R2=`ls *_2.fq | sed 's/.fq//g'`; 
  echo "$R1 $R2"; 
  /storage/apps/bbmap/bbduk.sh -Xmx6g in1=$R1.fq in2=$R2.fq out1=$R1\_bbmap_adaptertrimmed.fq out2=$R2\_bbmap_adaptertrimmed.fq ref=/storage/apps/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=30 minavgquality=30; 
  cd ..; 
done

mkdir ../2_CleanData

# move the adapter trimmed files into another directory
mv */*bbmap_adaptertrimmed.fq /storage/data/DATA4/analysis/18_Cauris/2021/2_CleanData

cd ../2_CleanData
 
 # move to specific folder
for d in $(ls *.fq| awk -F_ '{print $1}' | sort -u); do echo $d; mkdir $d; mv "$d"_*.fq $d; done


# For generating read stats - Run this on both "Samples_without_Topup" and "Samples_with_Topup" folders
time for d in $(ls */*fq); do echo "perl /storage/apps/SNP_Validation_Scripts/tools/Q20_Q30_Stats_wo_functions.pl $d"; done >parallel_script_Q20_Q30_stats.sh
time parallel < parallel_script_Q20_Q30_stats.sh &

# Run fastqc on the raw reads (Env + 11F)
time for d in $(ls -d */| sed 's/\///g'); do echo $d; cd $d; /storage/apps/FastQC/fastqc -t 48 *.fq; cd ..; done

#combine all fastqc reports using multiqc
cd 1_Rawdata

multiqc .

# Combine the read stats from Raw and Clean reads into one file to later use to create plot in R
$ grep "" */Q30_Q20_readstats.txt | sed 's/\/Q30_Q20_readstats.txt:/  /g' | sed 's/\//        /g' | sort -nk1,1 -nk2,2 >1_Raw_Clean_Q30_Q20_readstats.txt
```


### Step 3: Variant calling using SNIPPY


```bash

$ conda activate snippy_4.4.3

$ for d in $(ls -d */ | tr -d "/"); do echo "$d"; find /storage/data/DATA4/analysis/18_Cauris/2021/2_CleanData -name "$d*_bbmap*q" | sort; done | paste - - - >>input.tab

$ snippy-multi input.tab --ref /storage/data/DATA4/analysis/18_Cauris/2021/3_SNIPPY/B8441_PacBio_Assembly_PEKT02.1.fasta --cpus 48 > runme.sh

$ time bash runme.sh &
snippy-multi input.tab --ref /storage/data/DATA4/analysis/18_Cauris/2021/3_SNIPPY/B8441_PacBio_Assembly_PEKT02.1.fasta --cpus 48 > runme.sh

$ /storage/apps/snippy/bin/snippy-clean_full_aln core.full.aln > clean.full.aln
$ time run_gubbins.py --prefix postGubbins --filter_percentage 100 --threads 48 clean.full.aln --verbose

# fasttree generates a ML tree with bootstrap 
time FastTree -gtr -nt postGubbins.filtered_polymorphic_sites.fasta > postGubbins.filtered_polymorphic_sites.fasttree.tree

# Just running the below commands to compare
#$ snp-sites -c postGubbins.filtered_polymorphic_sites.fasta > postGubbins.filtered_polymorphic_sites.snp-sites.fasta ##removes the columns with N's
$ FastTree -gtr -nt clean.core.aln > clean.core.tree

```



```r
library("ggplot2")
library("ggtree")

setwd("/data02/Analysis/Projects/11_Cauris_2021/woRef")
tree <- read.tree(file = "postGubbins.filtered_polymorphic_sites.fasttree.tree")

p <- ggtree(tree) + 
  #xlim(0, 1) + # to allow more space for labels
  geom_treescale() # adds the scale

#p

tipcategories = read.csv("tree_meta.csv", 
                         sep = ",",
                         col.names = c("seq", "cat"), 
                         header = FALSE, 
                         stringsAsFactors = FALSE)
dd = as.data.frame(tipcategories)
#dd

p <- p %<+% dd + geom_tippoint() +
  geom_tiplab(align=F,aes(color=cat)) +
  labs(x = "",
       y = "",
       title = "ML Tree of Candida auris samples") +
    theme(legend.position = c(0.9,0.2), 
        legend.title = element_blank(), # no title
       # legend.key = element_blank()
       )
p 
```

![](Cauris_datareport_files/figure-html/chr_snps_presence-1.png)<!-- -->

```r
# Without aligned tips

p3 <- ggtree(tree,branch.length="none",layout="rectangular") +
  #xlim(0, 90) + 
  geom_tiplab(align=T,size=2, color="black") +
  geom_label2(aes(subset=!isTip, label=node), size=4, color="darkred", alpha=0.5) +
 # theme_tree2("white") +
  scale_color_continuous(low='white', high='hotpink', name="Branch length (my)") +
  theme(legend.position="bottom")
plot(p3)
```

![](Cauris_datareport_files/figure-html/chr_snps_presence-2.png)<!-- -->