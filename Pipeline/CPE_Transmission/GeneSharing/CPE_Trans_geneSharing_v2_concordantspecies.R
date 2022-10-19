library(dplyr)
library(UpSetR)
library(tidyverse)
library(splitstackshape)
library(reshape2)
library(randomcoloR)
library(ComplexHeatmap)
library(data.table)
library(kableExtra)



setwd("/data02/Analysis/Projects/2_CPE_Transmission/VennDiagram_of_genes_between_species/CPE_Transmission_Data_Analysis")


hybridAssemb_1198list <- readLines("/data02/Analysis/Projects/2_CPE_Transmission/VennDiagram_of_genes_between_species/CPE_Transmission_Data_Analysis/CPE_Hyrbid_Assemblies_1198.list")
head(hybridAssemb_1198list)

###@@@@@@@@@@@@@@ Dataframe1 - PlasClass Predictions

plasClass_Prob_df <- data.table::fread(file = "1_plasClass_probs/CPE_Trans_PlasClass_Predictions_comb.tab", sep = "\t", header = FALSE)
head(plasClass_Prob_df)

colnames(plasClass_Prob_df) <- c("Fasta", "Contig", "Probability")

plasClass_Prob_df <- 
  plasClass_Prob_df %>%  
  mutate(Fasta = stringr::str_replace(Fasta, ".plasclass.probs.out", "")) %>% 
  mutate(Fasta_Contig=paste(Fasta, Contig, sep = '#')) %>% 
  mutate(Classification = if_else(Probability >= 0.6, "Plasmid", "Chromosome")) 

head(plasClass_Prob_df)

###@@@@@@@@@@@@@@ Dataframe2 - MLST

mlst_df <- data.table::fread(file = "5_MLST/mlst_log2.19", sep = "\t", header = FALSE, fill = TRUE)
head(mlst_df) # mlst_df now has 1289 assembly information
colnames(mlst_df) <- c("Fasta", "Species", "ST", "Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7")
mlst_df <- mlst_df %>% filter(Fasta %in% hybridAssemb_1198list) # mlst_df now has only 1198 assembly information

Species_Counts <- as.data.frame(table(mlst_df$Species))
colnames(Species_Counts) <- c("Species","SampleCount")
head(Species_Counts)

#----------------Step1: Combining MLST and plasClass dataframes

plasClass_Prob_df <- left_join(plasClass_Prob_df, mlst_df, by = c("Fasta")) %>% 
    select(Fasta,Contig,Fasta_Contig,Classification,Species, ST) %>% 
    mutate(Species_ContigClass=paste(Species,Classification, sep = '#'))

head(plasClass_Prob_df)

###@@@@@@@@@@@@@@ Dataframe2 - Bacant Predictions

bacAnt_df <- data.table::fread(file = "2_bacant_annot/CPE_Trans_BacAnt_Transposons.tab", sep = "\t", header = TRUE)
head(bacAnt_df)

colnames(bacAnt_df) <- c("Fasta", "Contig", "Fasta_Contig", "Transposon")

#----------------Step2: Combining plasClass dataframe with bacant dataframe for plotting Transposons

bacant_plasClass_df <- left_join(bacAnt_df, plasClass_Prob_df, by = c("Fasta_Contig")) %>% 
  select(Fasta.x,Contig.x,Fasta_Contig,Transposon,Classification,Species_ContigClass)

head(bacant_plasClass_df)

bacant_plasClass_df <- bacant_plasClass_df %>% filter(Fasta.x %in% hybridAssemb_1198list) # bacant_plasClass_df now has only 1198 assembly information

bacant_plasClass_df %>%
  select(Transposon,Species_ContigClass) %>%
  table() %>% # frequency of unique transposons
  kbl(caption = "Transposons and their counts in Chromsomes and Plasmids across Species") %>%
  kable_classic(full_width = F, html_font = "Cambria")

# bacant_plasClass_df %>%
#  # filter(Transposon=="Tn5403|EU287476" & Species_ContigClass=="cfreundii#Plasmid") %>% 
#   select(Fasta.x,Transposon,Species_ContigClass) %>%
#   unique() %>% 
#   select(Transposon,Species_ContigClass) %>% 
#   table() %>% 
#   kbl(caption = "Sample counts for the Transposons in Chromsomes and Plasmids across Species") %>% 
#   kable_classic(full_width = F, html_font = "Cambria")

# percentages of transposons
bacant_plasClass_df %>% 
  select(Fasta.x,Transposon,Species_ContigClass) %>%
  unique() %>% 
  select(Transposon,Species_ContigClass) %>% 
  count(Transposon, Species_ContigClass) %>% 
  rename(Freq = n) %>%  
  drop_na(Species_ContigClass) %>% 
  mutate(Species_ContigClass2=Species_ContigClass) %>% 
  separate(Species_ContigClass2, sep = "#", into = c("Species", "ContigClass")) %>% 
  left_join(Species_Counts, by = c("Species")) %>% 
  mutate(SamplePrcnt_wGeneElement=Freq*100/SampleCount) %>% 
  select(Transposon,Species_ContigClass,SamplePrcnt_wGeneElement) %>% 
  data.table::dcast(Transposon~Species_ContigClass, value.var="SamplePrcnt_wGeneElement") %>% 
  replace(is.na(.), 0) %>% 
  mutate_if(is.numeric, round, digits=1) %>% 
  kbl(caption = "Sample count percentages for the Transposons in Chromsomes and Plasmids across Species") %>% 
  kable_classic(full_width = F, html_font = "Cambria")
  
head(bacant_plasClass_df)

bacant_plasClass_df %>% 
  select(Transposon,Species_ContigClass) %>% 
  unique() %>% 
  table() %>% # frequency of unique transposons
  kbl(caption = "Transposons and their counts in Chromsomes and Plasmids across Species") %>% 
  kable_classic(full_width = F, html_font = "Cambria")

#----------------Step3: UpsetR Plot Transposons
 
tn_upsetR_df <- bacant_plasClass_df %>% select(Transposon,Species_ContigClass)
#View(tn_upsetR_df)
tn_upsetR_lt <- split(tn_upsetR_df$Transposon, tn_upsetR_df$Species_ContigClass) 
tn_comb_mat = make_comb_mat(tn_upsetR_lt)

tn_comb_mat
#str(tn_comb_mat)

as.data.frame(tn_comb_mat)
# comb_mat
 UpSet(tn_comb_mat)
# UpSet(t(comb_mat))
# UpSet(comb_mat, top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
# right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE))

 
###---------------Setting a color for upset plot rows-----------###  
# set.seed(123)
# 
# upset_row_col <- as.data.frame(stack(set_size(tn_comb_mat))) %>% 
#   arrange(ind) %>% 
#   mutate(Species = stringr::str_replace(ind, "#.*", "")) %>% 
#    arrange(values) %>% 
#   group_by(Species) %>% 
#   dplyr::mutate(Color = randomcoloR::randomColor(length(unique(Species,luminosity="light")))) %>%
#  # dplyr::mutate(Color = ) %>%
#   dplyr::ungroup() %>%  as.data.frame()
# #class(upset_row_col)
# upset_row_col$Color
# UpSet(tn_comb_mat, 
#       row_title = "Species#SeqClassification", column_title = "Transposon Intersection across species",
#       pt_size = unit(6, "pt"),
#       lwd = unit(2, "pt"),
#     # comb_col = "red",
#      bg_col = upset_row_col$Color,
#     #  bg_pt_col = "#CCCCFF",
#       top_annotation = upset_top_annotation(tn_comb_mat, gp = gpar(fill = "#009797"), add_numbers = TRUE,bar_width = 0.5, annotation_name_rot = 90),
#       right_annotation = upset_right_annotation(tn_comb_mat, add_numbers = TRUE, gp = gpar(fill = "#009797")),
#       show_row_names = TRUE,
#       row_names_gp = gpar(fontsize = 9), #changes font size of "set size" labels
#       width = unit(700, units = "pt"), height = unit(10, "cm"))

col_size = comb_size(tn_comb_mat)
row_size = set_size(tn_comb_mat)

ht = UpSet(tn_comb_mat, 
      row_title = "Species#SeqClassification", column_title = "Transposon Intersection across species",
      pt_size = unit(6, "pt"),
      lwd = unit(2, "pt"),
      #comb_col = "red",
       bg_col = c("#d6b23b","#81dce2","#acf49a","#a784d8","#acf49a","#81dce2","#1a51dd","#1d7a01","#683cd8","#d62c4c","#683cd8","#1d7a01" ,
      "#dd71c0","#1a51dd","#dd71c0","#d62c4c"),
      #bg_col = upset_row_col$Color,
      top_annotation = upset_top_annotation(tn_comb_mat, gp = gpar(fill = "#009797"), add_numbers = FALSE, 
                                            bar_width = 0.5, annotation_name_rot = 90),
      right_annotation = upset_right_annotation(tn_comb_mat, add_numbers = TRUE, gp = gpar(fill = "#009797")),
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 9), #changes font size of "set size" labels
      width = unit(700, units = "pt"), height = unit(10, "cm"))

ht = draw(ht)
ht
col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], 
            seq_len(length(col_size)), 
            unit(col_size[col_od], "native") + unit(2, "mm"), 
            default.units = "native", just = "bottom",
            gp = gpar(fontsize = 8))
})

###@@@@@@@@@@@@@@ Dataframe3 - AMRfinderplus - AMR Genes

amrFinderplus_df <- data.table::fread(file = "3_amrfinderplus_results/CPE_Trans_AMRFinderplus_AMR_Genes.tab", sep = "\t", header = TRUE)
head(amrFinderplus_df)

colnames(amrFinderplus_df) <- c("Fasta", "Contig", "AMR_GENE", "AMR_GENE_CLASS","Coverage","Identity")

amrFinderplus_df <- amrFinderplus_df %>% filter(Fasta %in% hybridAssemb_1198list) # bacant_plasClass_df now has only 1198 assembly information

amrFinderplus_df <- 
  amrFinderplus_df %>%  
  mutate(Fasta_Contig=paste(Fasta, Contig, sep = '#'))

#----------------Step2: Combining plasClass dataframe with AMRFinderplus dataframe for plotting AMR genes

amrFinderplus_df <- left_join(amrFinderplus_df, plasClass_Prob_df, by = c("Fasta_Contig")) %>% 
  select(Fasta.x,Contig.x,AMR_GENE,AMR_GENE_CLASS,Classification,Species_ContigClass)


# percentages of AMR_GENEs
amrFinderplus_df %>% 
  select(Fasta.x,AMR_GENE,Species_ContigClass) %>%
  unique() %>% 
  select(AMR_GENE,Species_ContigClass) %>% 
  count(AMR_GENE, Species_ContigClass) %>% 
  rename(Freq = n) %>%  
  drop_na(Species_ContigClass) %>% 
  mutate(Species_ContigClass2=Species_ContigClass) %>% 
  separate(Species_ContigClass2, sep = "#", into = c("Species", "ContigClass")) %>% 
  left_join(Species_Counts, by = c("Species")) %>% 
  mutate(SamplePrcnt_wGeneElement=Freq*100/SampleCount) %>% 
  select(AMR_GENE,Species_ContigClass,SamplePrcnt_wGeneElement) %>% 
  data.table::dcast(AMR_GENE~Species_ContigClass, value.var="SamplePrcnt_wGeneElement") %>% 
  replace(is.na(.), 0) %>% 
  mutate_if(is.numeric, round, digits=1) %>% 
  kbl(caption = "Sample count percentages for the AMR GENEs in Chromsomes and Plasmids across Species") %>% 
  kable_classic(full_width = F, html_font = "Cambria")


#----------------Step3: UpsetR Plot AMR Genes

amr_upsetR_df <- amrFinderplus_df %>% select(AMR_GENE,Species_ContigClass)
View(amr_upsetR_df)
amr_upsetR_lt <- split(amr_upsetR_df$AMR_GENE, amr_upsetR_df$Species_ContigClass) 
amr_comb_mat = make_comb_mat(amr_upsetR_lt)

col_size = comb_size(amr_comb_mat)
row_size = set_size(amr_comb_mat)

ht = UpSet(amr_comb_mat, 
           row_title = "Species#SeqClassification", column_title = "AMR Genes Intersection across species",
           pt_size = unit(6, "pt"),
           lwd = unit(2, "pt"),
           #comb_col = "red",
           # bg_col = c("#d6b23b","#81dce2","#acf49a","#a784d8","#acf49a","#81dce2","#1a51dd","#1d7a01","#683cd8","#d62c4c","#683cd8","#1d7a01" ,
           #            "#dd71c0","#1a51dd","#dd71c0","#d62c4c"),
           #bg_col = upset_row_col$Color,
           top_annotation = upset_top_annotation(amr_comb_mat, gp = gpar(fill = "#009797"), add_numbers = FALSE, 
                                                 bar_width = 0.5, annotation_name_rot = 90),
           right_annotation = upset_right_annotation(amr_comb_mat, add_numbers = TRUE, gp = gpar(fill = "#009797")),
           show_row_names = TRUE,
           row_names_gp = gpar(fontsize = 9), #changes font size of "set size" labels
           width = unit(900, units = "pt"), height = unit(12, "cm"))

ht = draw(ht)
ht
col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], 
            seq_len(length(col_size)), 
            unit(col_size[col_od], "native") + unit(2, "mm"), 
            default.units = "native", just = "bottom",
            gp = gpar(fontsize = 8))
})

###@@@@@@@@@@@@@@ Dataframe4 - Abricate - Virulence Factor Genes

Abricate_vf_df <- data.table::fread(file = "4_abricate_vf_results/CPE_Trans_Abricate_VF_Genes.tab", sep = "\t", header = TRUE)
head(Abricate_vf_df)

##########-------------- Modifying dataframe START -------------------################

Abricate_vf_df <- Abricate_vf_df %>% select(`#FILE`,SEQUENCE,START,END,STRAND,GENE)  
colnames(Abricate_vf_df ) <- c("Fasta", "Contig","START","END","STRAND","VF_GENE")

Abricate_vf_df <- 
  Abricate_vf_df %>%  
  mutate(Fasta_Contig=paste(Fasta, Contig, sep = '#'))

head(Abricate_vf_df)

#----------------Step2: Combining plasClass dataframe with abricate_vf_dataframe for plotting virulence genes

Abricate_vf_df <- left_join(Abricate_vf_df, plasClass_Prob_df, by = c("Fasta_Contig")) %>% 
  select(Fasta.x,Contig.x,VF_GENE,Classification,Species_ContigClass)

# percentages of AMR_GENEs
Abricate_vf_df %>% 
  select(Fasta.x,VF_GENE,Species_ContigClass) %>%
  unique() %>% 
  select(VF_GENE,Species_ContigClass) %>% 
  count(VF_GENE, Species_ContigClass) %>% 
  rename(Freq = n) %>%  
  drop_na(Species_ContigClass) %>% 
  mutate(Species_ContigClass2=Species_ContigClass) %>% 
  separate(Species_ContigClass2, sep = "#", into = c("Species", "ContigClass")) %>% 
  left_join(Species_Counts, by = c("Species")) %>% 
  mutate(SamplePrcnt_wGeneElement=Freq*100/SampleCount) %>% 
  select(VF_GENE,Species_ContigClass,SamplePrcnt_wGeneElement) %>% 
  data.table::dcast(VF_GENE~Species_ContigClass, value.var="SamplePrcnt_wGeneElement") %>% 
  replace(is.na(.), 0) %>% 
  mutate_if(is.numeric, round, digits=1) %>% 
  kbl(caption = "Sample count percentages for the AMR GENEs in Chromsomes and Plasmids across Species") %>% 
  kable_classic(full_width = F, html_font = "Cambria")

#----------------Step3: UpsetR Plot VFGenes

abr_upsetR_df <- Abricate_vf_df %>% select(VF_GENE,Species_ContigClass)
#View(abr_upsetR_df)

abr_upsetR_lt <- split(abr_upsetR_df$VF_GENE, abr_upsetR_df$Species_ContigClass) 
abr_comb_mat = make_comb_mat(abr_upsetR_lt)

col_size = comb_size(abr_comb_mat)
row_size = set_size(abr_comb_mat)

ht = UpSet(abr_comb_mat, 
           row_title = "Species#SeqClassification", column_title = "AMR Genes Intersection across species",
           pt_size = unit(6, "pt"),
           lwd = unit(2, "pt"),
           #comb_col = "red",
           # bg_col = c("#d6b23b","#81dce2","#acf49a","#a784d8","#acf49a","#81dce2","#1a51dd","#1d7a01","#683cd8","#d62c4c","#683cd8","#1d7a01" ,
           #            "#dd71c0","#1a51dd","#dd71c0","#d62c4c"),
           #bg_col = upset_row_col$Color,
           top_annotation = upset_top_annotation(abr_comb_mat, gp = gpar(fill = "#009797"), add_numbers = FALSE, 
                                                 bar_width = 0.5, annotation_name_rot = 90),
           right_annotation = upset_right_annotation(abr_comb_mat, add_numbers = TRUE, gp = gpar(fill = "#009797")),
           show_row_names = TRUE,
           row_names_gp = gpar(fontsize = 9), #changes font size of "set size" labels
           width = unit(900, units = "pt"), height = unit(12, "cm"))

ht = draw(ht)
ht
col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("intersection_size", {
  grid.text(col_size[col_od], 
            seq_len(length(col_size)), 
            unit(col_size[col_od], "native") + unit(2, "mm"), 
            default.units = "native", just = "bottom",
            gp = gpar(fontsize = 8))
})
