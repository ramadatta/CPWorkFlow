---
title: "Data Visualization of Pre and Post Trimmmed Nanopore reads: Statistics from bbtools statswrapper"
author: "Prakki Sai Rama Sridatta"
date: "2021-06-16"
output:
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---



### Step 1: Filter Long reads 

All the raw fastq reads for the samples are processed using the tool "filtlong" using the following command.


```bash
time for d in $(ls *.fastq | sed 's/.fastq//g'); 
do 
filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 "$d".fastq > "$d".filtlong.fastq & 
done
```

### Step 2: Generate read statistics
Now that we have both raw and filtered fastq files run bbtools statswrapper.sh


```bash
time statswrapper.sh *.fastq >bbstats_wrapper.txt 
```

### Step 3: Generate R plots
Let us load the **bbstats_wrapper.txt** file into R and generate some figures from the data using R ggplot function.

```r
setwd("/data02/Analysis/Projects/8_Aqueos_samples/11F/bbstats/")

library(dplyr)
library(ggplot2)
library(ggunchained)
library(kableExtra)

bbstats <- read.table("bbstats_wrapper.txt", sep = "\t", header = TRUE)

# Since the table is very big outputting the contents into table with scroll bar
head(bbstats) %>%
  kbl() %>%
   kable_paper() %>%
  scroll_box(width = "100%", height = "200px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:200px; overflow-x: scroll; width:100%; "><table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> n_scaffolds </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> n_contigs </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_bp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> contig_bp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> gap_pct </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_L50 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_N50 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ctg_L50 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ctg_N50 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_L90 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_N90 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ctg_L90 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ctg_N90 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_max </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ctg_max </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_n_gt50K </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> scaf_pct_gt50K </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> gc_avg </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> gc_std </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> filename </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Sample </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Path </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 7518 </td>
   <td style="text-align:right;"> 7518 </td>
   <td style="text-align:right;"> 82157895 </td>
   <td style="text-align:right;"> 82157895 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1198 </td>
   <td style="text-align:right;"> 20696 </td>
   <td style="text-align:right;"> 1198 </td>
   <td style="text-align:right;"> 20696 </td>
   <td style="text-align:right;"> 4189 </td>
   <td style="text-align:right;"> 5129 </td>
   <td style="text-align:right;"> 4189 </td>
   <td style="text-align:right;"> 5129 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:right;"> 135 </td>
   <td style="text-align:right;"> 10.168 </td>
   <td style="text-align:right;"> 0.54364 </td>
   <td style="text-align:right;"> 0.04361 </td>
   <td style="text-align:left;"> A1089_NP </td>
   <td style="text-align:left;"> A1089 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A1089_Nanopore.fastq </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4475 </td>
   <td style="text-align:right;"> 4475 </td>
   <td style="text-align:right;"> 73944083 </td>
   <td style="text-align:right;"> 73944083 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1014 </td>
   <td style="text-align:right;"> 22956 </td>
   <td style="text-align:right;"> 1014 </td>
   <td style="text-align:right;"> 22956 </td>
   <td style="text-align:right;"> 3176 </td>
   <td style="text-align:right;"> 8100 </td>
   <td style="text-align:right;"> 3176 </td>
   <td style="text-align:right;"> 8100 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:right;"> 133 </td>
   <td style="text-align:right;"> 11.121 </td>
   <td style="text-align:right;"> 0.54459 </td>
   <td style="text-align:right;"> 0.03822 </td>
   <td style="text-align:left;"> A1089_NP.filtlong </td>
   <td style="text-align:left;"> A1089 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A1089_Nanopore.filtlong.fastq </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 135106 </td>
   <td style="text-align:right;"> 135106 </td>
   <td style="text-align:right;"> 1442374512 </td>
   <td style="text-align:right;"> 1442374512 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 21266 </td>
   <td style="text-align:right;"> 19686 </td>
   <td style="text-align:right;"> 21266 </td>
   <td style="text-align:right;"> 19686 </td>
   <td style="text-align:right;"> 76962 </td>
   <td style="text-align:right;"> 5405 </td>
   <td style="text-align:right;"> 76962 </td>
   <td style="text-align:right;"> 5405 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:right;"> 2473 </td>
   <td style="text-align:right;"> 11.611 </td>
   <td style="text-align:right;"> 0.54790 </td>
   <td style="text-align:right;"> 0.05361 </td>
   <td style="text-align:left;"> A1416_NP </td>
   <td style="text-align:left;"> A1416 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A1416_Nanopore.fastq </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 18128 </td>
   <td style="text-align:right;"> 18128 </td>
   <td style="text-align:right;"> 500097148 </td>
   <td style="text-align:right;"> 500097148 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 5447 </td>
   <td style="text-align:right;"> 30520 </td>
   <td style="text-align:right;"> 5447 </td>
   <td style="text-align:right;"> 30520 </td>
   <td style="text-align:right;"> 14425 </td>
   <td style="text-align:right;"> 15937 </td>
   <td style="text-align:right;"> 14425 </td>
   <td style="text-align:right;"> 15937 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:right;"> 1461 </td>
   <td style="text-align:right;"> 19.792 </td>
   <td style="text-align:right;"> 0.55250 </td>
   <td style="text-align:right;"> 0.03752 </td>
   <td style="text-align:left;"> A1416_NP.filtlong </td>
   <td style="text-align:left;"> A1416 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A1416_Nanopore.filtlong.fastq </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 25032 </td>
   <td style="text-align:right;"> 25032 </td>
   <td style="text-align:right;"> 349105749 </td>
   <td style="text-align:right;"> 349105749 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4072 </td>
   <td style="text-align:right;"> 25524 </td>
   <td style="text-align:right;"> 4072 </td>
   <td style="text-align:right;"> 25524 </td>
   <td style="text-align:right;"> 14526 </td>
   <td style="text-align:right;"> 7422 </td>
   <td style="text-align:right;"> 14526 </td>
   <td style="text-align:right;"> 7422 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:right;"> 1039 </td>
   <td style="text-align:right;"> 19.769 </td>
   <td style="text-align:right;"> 0.55321 </td>
   <td style="text-align:right;"> 0.04948 </td>
   <td style="text-align:left;"> A472_NP </td>
   <td style="text-align:left;"> A472 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A472_Nanopore.fastq </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16698 </td>
   <td style="text-align:right;"> 16698 </td>
   <td style="text-align:right;"> 314218586 </td>
   <td style="text-align:right;"> 314218586 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 3502 </td>
   <td style="text-align:right;"> 27360 </td>
   <td style="text-align:right;"> 3502 </td>
   <td style="text-align:right;"> 27360 </td>
   <td style="text-align:right;"> 11711 </td>
   <td style="text-align:right;"> 8186 </td>
   <td style="text-align:right;"> 11711 </td>
   <td style="text-align:right;"> 8186 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:right;"> 1003 </td>
   <td style="text-align:right;"> 21.218 </td>
   <td style="text-align:right;"> 0.55385 </td>
   <td style="text-align:right;"> 0.04438 </td>
   <td style="text-align:left;"> A472_NP.filtlong </td>
   <td style="text-align:left;"> A472 </td>
   <td style="text-align:left;"> /storage/data/DATA4/analysis/30_Aqueous_Env_11F_data_analysis/2_AdapterTrimmed_bbduk_Q30/nanopore_qcat_adaptertrimmed_filtlong/11F/A472_Nanopore.filtlong.fastq </td>
  </tr>
</tbody>
</table></div>

```r
# Subsetting columns to a smaller dataframe for ease of analysis
bbstats_filtered <- bbstats %>% select(n_contigs, contig_bp, ctg_N50,ctg_max,filename,Sample)

#renaming column names
colnames(bbstats_filtered) <- c("ReadCount","Totalbases", "N50","LongestRead","Filename", "Sample")

# Since the table is small outputting the contents into kable table
 head(bbstats_filtered) %>%
  kbl() %>%
  kable_classic_2(full_width = F)
```

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;"> ReadCount </th>
   <th style="text-align:right;"> Totalbases </th>
   <th style="text-align:right;"> N50 </th>
   <th style="text-align:right;"> LongestRead </th>
   <th style="text-align:left;"> Filename </th>
   <th style="text-align:left;"> Sample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 7518 </td>
   <td style="text-align:right;"> 82157895 </td>
   <td style="text-align:right;"> 20696 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:left;"> A1089_NP </td>
   <td style="text-align:left;"> A1089 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4475 </td>
   <td style="text-align:right;"> 73944083 </td>
   <td style="text-align:right;"> 22956 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:left;"> A1089_NP.filtlong </td>
   <td style="text-align:left;"> A1089 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 135106 </td>
   <td style="text-align:right;"> 1442374512 </td>
   <td style="text-align:right;"> 19686 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:left;"> A1416_NP </td>
   <td style="text-align:left;"> A1416 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 18128 </td>
   <td style="text-align:right;"> 500097148 </td>
   <td style="text-align:right;"> 30520 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:left;"> A1416_NP.filtlong </td>
   <td style="text-align:left;"> A1416 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 25032 </td>
   <td style="text-align:right;"> 349105749 </td>
   <td style="text-align:right;"> 25524 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:left;"> A472_NP </td>
   <td style="text-align:left;"> A472 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16698 </td>
   <td style="text-align:right;"> 314218586 </td>
   <td style="text-align:right;"> 27360 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:left;"> A472_NP.filtlong </td>
   <td style="text-align:left;"> A472 </td>
  </tr>
</tbody>
</table>

```r
# Parsing data  
bbstats_filtered <- bbstats_filtered %>% 
      mutate(data = case_when(grepl("_NP.filtlong", Filename) ~ "Filtered",
                              grepl("_NP", Filename, ignore.case = TRUE) ~"Raw")) %>% 
          select(-Filename) 

# Since the table is small outputting the contents into kable table
 head(bbstats_filtered) %>%
  kbl() %>%
  kable_classic_2(full_width = F)
```

<table class=" lightable-classic-2" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;"> ReadCount </th>
   <th style="text-align:right;"> Totalbases </th>
   <th style="text-align:right;"> N50 </th>
   <th style="text-align:right;"> LongestRead </th>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:left;"> data </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 7518 </td>
   <td style="text-align:right;"> 82157895 </td>
   <td style="text-align:right;"> 20696 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:left;"> A1089 </td>
   <td style="text-align:left;"> Raw </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4475 </td>
   <td style="text-align:right;"> 73944083 </td>
   <td style="text-align:right;"> 22956 </td>
   <td style="text-align:right;"> 119488 </td>
   <td style="text-align:left;"> A1089 </td>
   <td style="text-align:left;"> Filtered </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 135106 </td>
   <td style="text-align:right;"> 1442374512 </td>
   <td style="text-align:right;"> 19686 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:left;"> A1416 </td>
   <td style="text-align:left;"> Raw </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 18128 </td>
   <td style="text-align:right;"> 500097148 </td>
   <td style="text-align:right;"> 30520 </td>
   <td style="text-align:right;"> 167349 </td>
   <td style="text-align:left;"> A1416 </td>
   <td style="text-align:left;"> Filtered </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 25032 </td>
   <td style="text-align:right;"> 349105749 </td>
   <td style="text-align:right;"> 25524 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:left;"> A472 </td>
   <td style="text-align:left;"> Raw </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16698 </td>
   <td style="text-align:right;"> 314218586 </td>
   <td style="text-align:right;"> 27360 </td>
   <td style="text-align:right;"> 153451 </td>
   <td style="text-align:left;"> A472 </td>
   <td style="text-align:left;"> Filtered </td>
  </tr>
</tbody>
</table>

```r
#N50 Dumbbell plot 
ggplot(bbstats_filtered, aes(x=N50, y=reorder(Sample,N50))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="N50",
       y="Sample",
       title = "N50 Pre and Post trimming (Nanopore)")
```

![](bbstats_files/figure-html/bb-1.png)<!-- -->

```r
#TotalBases Dumbbell plot
ggplot(bbstats_filtered, aes(x=Totalbases, y=reorder(Sample,Totalbases))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="Totalbases",
       y="Sample",
       title = "Totalbases Pre and Post trimming (Nanopore)")
```

![](bbstats_files/figure-html/bb-2.png)<!-- -->

```r
#ReadCounts Dumbbell plot
ggplot(bbstats_filtered, aes(x=ReadCount, y=reorder(Sample,ReadCount))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="ReadCount",
       y="Sample",
       title = "ReadCount Pre and Post trimming (Nanopore)")
```

![](bbstats_files/figure-html/bb-3.png)<!-- -->

```r
#LongestRead Dumbbell plot
ggplot(bbstats_filtered, aes(x=LongestRead, y=reorder(Sample,LongestRead))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="Read Length",
       y="Sample",
       title = "LongestRead Pre and Post trimming (Nanopore)")
```

![](bbstats_files/figure-html/bb-4.png)<!-- -->

## Easy Way to combine multiple figures into single figure

We will write the all the ggplot functions into respective ggplot variable and then provides all the plots to the function ggarrange() [in ggpubr] as below. We can create a common unique legend for multiple plots.



```r
#N50 Dumbbell plot 
N50_dumb <- ggplot(bbstats_filtered, aes(x=N50, y=reorder(Sample,N50))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="N50",
       y="Sample",
       title = "N50 Pre and Post trimming (Nanopore)")
#class(N50_dumb)

#TotalBases
TB_dumb <- ggplot(bbstats_filtered, aes(x=Totalbases, y=reorder(Sample,Totalbases))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="Totalbases",
       y="Sample",
       title = "Totalbases Pre and Post trimming (Nanopore)")

#ReadCounts
RC_dumb <- ggplot(bbstats_filtered, aes(x=ReadCount, y=reorder(Sample,ReadCount))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="ReadCount",
       y="Sample",
       title = "ReadCount Pre and Post trimming (Nanopore)")

#LongestRead
LongCong_dumb <- ggplot(bbstats_filtered, aes(x=LongestRead, y=reorder(Sample,LongestRead))) + 
  geom_line(aes(group = Sample))+
  geom_point(aes(color=data), size=1) +
  theme_light()+
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size = 6), axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  scale_color_brewer(palette="Dark2", direction=-1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  labs(x="Read Length",
       y="Sample",
       title = "LongestRead Pre and Post trimming (Nanopore)")

library(ggpubr)

ggarrange(N50_dumb, TB_dumb, RC_dumb, LongCong_dumb, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,common.legend = TRUE)
```

![](bbstats_files/figure-html/ggarrange-1.png)<!-- -->

```r
# To save in a high quality figure can use the following command
ggsave("bbstats_nanopore.png",dpi = 300, width = 20, height = 10, units = "in")
```
