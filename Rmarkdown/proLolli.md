---
title: "proLolli"
author: "Prakki Sai Rama Sridatta"
date: "2021-08-20"
output:
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---



### Step 1: Reference genome annotation

For this step, I have annotated my reference genome using Prokka.


```bash
prokka2 --force --cpus 32 --outdir  SHG10_pKPC2-km_prokka_out --prefix  SHG10_pKPC2-km  SHG10_pKPC2-km_unicycler.gte1kb.contigs.fasta >>prokka_log 2>>prokka_error; 
```

### Step 2: Generate SNPs from using multiple samples

Using the above reference, I have run SNIPPY for multiple samples


```bash
snippy commands
```

### Step 3: Extract SNPs

SNIPPY generates a VCF file for each sample and we want the SNPs.
As we can see below 3 WT files have around 10-16 SNPs


```bash
$ grep "" *WT*/snps.filt.vcf  | grep -v '#' | cut -f1 -d":" | sort | uniq -c
     15 SGH10_WT_G300__1/snps.filt.vcf
     16 SGH10_WT_G300__2/snps.filt.vcf
     10 SGH10_WT_G300__3/snps.filt.vcf
```

### Step 4: Extract SNPs from Single VCF file

Let us use SNPs from `SGH10_WT_G300__3` sample.


```bash
cat snps.filt.vcf | grep -v "#" | awk '{print $1 "\t" $2}' >SGH10_WT_G300__2_snps.list
```

### Step 5: Generate Lollipop Plot R script

I wrote a shell script which generate a R script to generate lollipop plot. 
Basically using SNPs from `SGH10_WT_G300__2_snps.list` file, let us bedtools intersect and extract the regions where the SNPs are falling under.

Once the regions and coordinates are found, we pass this to `trackViewer` which generates lollipop plot


```bash
bash lollipop_rcode_generator.sh >proLolli.R
```

### Step 6: Generate Lollipop Plot using R script


```r
setwd("/home/prakki/Documents/LeaRn/proLolli")

source(file = "proLolli.R")
```

![](proLolli_files/figure-html/snp1-1.png)<!-- -->

