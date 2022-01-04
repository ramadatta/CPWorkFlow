# This script will take the prokka gff files and can generate the +N and -N flanking genes in the gggenes format


setwd("/home/prakki/Desktop/renamed/Prokka_Annotation/Local_Prokka")
library(dplyr)
library(stringr)
library(tidyverse)
local_df <- read.table("Local_PAE_genomes.gff", sep = "\t", header = FALSE, quote = "")
head(local_df)

colnames(local_df) <- c("Sample_Ctg","Prodigal_Annotation","CDS","START","END","Dot","Direction","UnknownCol","GeneProduct")
head(local_df)

local_df_subset <- local_df %>% select(Sample_Ctg,START,END,Direction,GeneProduct)
head(local_df_subset)

local_df_subset_mod <- local_df_subset %>%
  mutate(Strand = case_when(
    grepl('\\+', Direction) ~ 'forward',
    grepl('\\-', Direction) ~ 'reverse')) %>%  
  mutate(Direction2 = case_when(
    grepl('\\+', Direction) ~ '1',
    grepl('\\-', Direction) ~ '-1')) %>% 
  select(-Direction) %>% 
  rename(Direction = Direction2) %>% 
  mutate(gene = str_extract(GeneProduct, "product=.+$")) %>% 
  select(-GeneProduct) %>% 
  mutate(gene = str_replace(gene, "product=", "")) %>% 
  mutate(gene = str_replace_all(gene, " ", "_")) %>%
  mutate(Sample_Ctg = str_replace(Sample_Ctg, "_depth=.*", "")) %>% 
  filter(!grepl('hypothetical_protein', gene)) %>% 
  filter(!grepl('C02340_unicycler.fasta', Sample_Ctg)) %>% 
  select(Sample_Ctg,START,END,Direction,gene,Strand) 
  

colnames(local_df_subset_mod) <- c("molecule","start","end","direction","gene","strand")
head(local_df_subset_mod)
#view(local_df_subset_mod)        
#NDM_rows <- local_df_subset_mod %>% filter(grepl('NDM', gene)) %>% filter(!grepl('C02340_unicycler.fasta', molecule))
#NDM_rows

###### This generates dataframe required for the gggenes package ##########
# Gene list:
#quinolone_resistance_pentapeptide_repeat_protein_QnrVC1
#subclass_B1_metallo-beta-lactamase_NDM-1
target=c("quinolone_resistance_pentapeptide_repeat_protein_QnrVC1")
match = which(local_df_subset_mod$gene %in% target)
getThese = unique(as.vector(mapply(seq,match-10,match+10))) # number of rows before and after
getThese

getThese = getThese[getThese > 0 & getThese <= nrow(local_df_subset_mod)]
gene_flanking_coords <- local_df_subset_mod[getThese,]
gene_flanking_coords 
write.csv(NDM_flanking_coords,file = "PAE_local_gggenes_input_80.txt")

###### This generates list ##########
# library(SOfun)
# df <- getMyRows(local_df_subset_mod, which(local_df_subset_mod$gene == "subclass B1 metallo-beta-lactamase NDM-1"), -10:10, TRUE)
# head(df)
# 
# sink("mylist.txt")
# print(df)
# sink()




