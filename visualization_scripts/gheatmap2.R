library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(tidytree)
library(ggstar)
library(RColorBrewer)
library(reshape2)

Plot_PGT_SNPMat <- function(newickFile, meltfile, organism, offset_var, SNP_limits, legend_limits, wid=15, hei=10){
#organism <- "Serratia marcescens" # Serratia marcescens, Klebsiella Pneumoniae
#organism <- "Klebsiella Pneumoniae" 
#offset_var=3000 # Whereshould the heatmap start - if Maximum SNPs is 18000 I have set 3000 offset, if MaxSNps is 11, offset I have set is 1
#SNP_limits = c(0,100)
b <- legend_limits
#b <- c(0,10,20,30,40,50,60,70,100,1000,15000,20000)
#wid=15
#hei=7


tr <- read.tree(newickFile)
#tr <- read.tree("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/KP_postGubbins.final_tree.tre")
#tr <- read.tree("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/SM_postGubbins.final_tree.tre")
#tr <- read.tree("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/EC_postGubbins.final_tree.tre")
tr <- read.tree("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/ECLOC_postGubbins.final_tree.tre")
thetree <- ggtree(tr=tr, open.angle=15, size=0.2) #+ #removed this in ggtree command: , aes(colour=Clade)
thetree

cge_log <-  read.table("/data02/Analysis/Projects/8_Aqueos_samples/4_CGE_Resfinder4.1/cge_resistome_species_ST.txt",sep = "\t", header = TRUE, quote="")
ARG_Cat <- read.csv("/data02/Analysis/Projects/8_Aqueos_samples/4_CGE_Resfinder4.1/ARG_Categories_withESBL.csv",sep = ",", header = TRUE, check.names = FALSE)

sample_ARG_long <- as.data.frame(separate_rows(cge_log, resistance, sep=','))
head(sample_ARG_long)

sample_ARG_Cat_long <- as.data.frame(sample_ARG_long %>% right_join(ARG_Cat, by = 'resistance')%>% drop_na()) 
head(sample_ARG_Cat_long)

sample_ARG_Cat_long$SampleType[sample_ARG_Cat_long$SampleType == '11F'] <- 'Patient'
sample_ARG_Cat_long$SampleType[sample_ARG_Cat_long$SampleType == 'Env'] <- 'Environment'

#snp_table <- read.table(meltfile,sep = ",", header = TRUE)
#snp_table <- read.table("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix/postGubbins.filtered_polymorphic_sites_KP_full_triangle.csv", sep = ",", header = TRUE)
#snp_table <- read.table("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix/postGubbins.filtered_polymorphic_sites_SM_full_triangle.csv", sep = ",", header = TRUE)
#snp_table <- read.table("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix/postGubbins.filtered_polymorphic_sites_EC_full_triangle.csv", sep = ",", header = TRUE)
snp_table <- read.table("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix/postGubbins.filtered_polymorphic_sites_EntSp_full_triangle.csv", sep = ",", header = TRUE)

head(snp_table)

colnames(snp_table) <- c("No", "Sample1", "Sample2", "SNPDiff")
complete_snp_table <- na.omit(snp_table)
tail(complete_snp_table)
complete_snp_table[1] <- NULL
class(complete_snp_table)
complete_snp_table

## To generate a heatmap with SNP Difference values and sorted by the SNP values we will use pheatmap

library(reshape2)

wide_complete_snp_table <- acast(complete_snp_table, Sample1~Sample2, value.var="SNPDiff")
head(wide_complete_snp_table)

distances <- as.data.frame(as.matrix(wide_complete_snp_table))
head(distances)
class(distances)

class(thetree$data$label)
thetree_datasubset <- subset(thetree$data,thetree$data$isTip==TRUE)
thetree_datasubset
thetree_datasubset_sorted <- thetree_datasubset[order(thetree_datasubset$y, decreasing = TRUE),]
sorted_labels <- unlist(thetree_datasubset_sorted[,4])
head(sorted_labels)
reordered_distances <- distances[sorted_labels]
reordered_distances
class(reordered_distances)

# # lets add color for the tip points- based on the sampleType
# nb.cols <- length(unique(complete_snp_table$SNPDiff))
# nb.cols
# customPal <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols) # Set2, Accent, Paired
# #class(customPal)
# customPal



sampleType <- sample_ARG_Cat_long %>% select(sample_name,category,SampleType) %>%  unique()
head(sampleType)

sub_sampleType <- sampleType[sampleType$sample_name %in% thetree$data$label,]
class(thetree$data$label[!is.na(thetree$data$label)])
head(sub_sampleType)

# lets add color for the tip points- based on the sampleType
nbtip.cols <- length(unique(sample_ARG_Cat_long$SampleType))
nbtip.cols
#customPal <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols) # Set2, Accent, Paired
customPal.tip <- colorRampPalette(c("red", "#24868EFF")) # Set2, Accent, Paired
#class(customPal)
customPal.tip(nbtip.cols)

p1 <- thetree %<+% sub_sampleType +
  geom_tippoint(aes(colour=SampleType),
                alpha=1) +
  geom_tiplab(aes(colour=SampleType),
              align=TRUE,
              linetype=3,
              size=4,
              linesize=1, #def:0.2
              show.legend=FALSE
  ) +
  scale_colour_manual(
    name="SampleType",
    values=customPal.tip(nbtip.cols),
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=2,
                       override.aes=list(size=2,alpha=1))
  )
p1


#b <- c(0,10,20,30,40,50,60,70,100,1000,15000,20000)

p2 <- gheatmap(p1,reordered_distances, offset=offset_var, width=2,
               colnames_angle=90, colnames_offset_y = -0.2,font.size = 3) +
  # scale_fill_gradient(low="#ffd662ff", high="tomato3",name="genotype")
 # scale_fill_gradient(low="#ffd662ff", high="#24868EFF",name="genotype")
  scale_fill_gradientn(limits = SNP_limits, colours=c("#ffd662ff", "#24868EFF"),
                       breaks=b, labels=format(b), na.value="tomato3", name="SNP Difference") + 
  scale_fill_manual(values=NA) +              
  guides(colour=guide_legend("SNP Difference >100", override.aes=list(colour="tomato3"))) +
  labs(x="", y="", title = paste0(organism," Gubbins Phylogenetic Tree and SNP Difference Matrix"))

p2

#ggsave( paste0(organism,"_PGTree_SNPMat.pdf"), width = 70, height = 40, units = "cm", limitsize = FALSE)
ggsave(paste0(organism,"_PGTree_SNPMat.png"),dpi = 100, width = wid, height = hei, units = "in")
}

setwd("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis")

#KP

# Plot_PGT_SNPMat(
#   newickFile = "KP_postGubbins.final_tree.tre",
#   meltfile = "postGubbins.filtered_polymorphic_sites_KP_full_triangle.csv",
#   organism = "Klebsiella Pneumoniae",
#   offset_var = 3000,
#   SNP_limits = c(0, 100),
#   legend_limits = c(0, 10, 20, 30, 40, 50, 60, 70, 100, 1000, 15000, 20000),
#   wid = 15,
#   hei = 10
# )

#SM
Plot_PGT_SNPMat(
  newickFile = "SM_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_SM_full_triangle.csv",
  organism = "Serratia Marcescens",
  offset_var = 1,
  SNP_limits = c(0, 12),
  legend_limits = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
  wid = 18,
  hei = 13
)

#EC
Plot_PGT_SNPMat(
  newickFile = "EC_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_EC_full_triangle.csv",
  organism = "Escherichia Coli",
  offset_var = 1,
  SNP_limits = c(0, 12),
  legend_limits = c(0, 1, 2, 10, 1000, 4000),
  wid = 10,
  hei = 5
)

#ECLOC
Plot_PGT_SNPMat(
  newickFile = "ECLOC_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_EC_full_triangle.csv",
  organism = "Enterobacter Species",
  offset_var = 1,
  SNP_limits = c(0, 40),
  legend_limits = c(0, 1, 2, 10, 20, 30, 40, 100),
  wid = 12,
  hei = 5
)

