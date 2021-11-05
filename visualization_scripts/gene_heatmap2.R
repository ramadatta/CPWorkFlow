library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(tidytree)
library(ggstar)
library(RColorBrewer)
library(reshape2)
library(tidyverse)

setwd("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/PGTree_Gene_PAMatrix")

cge_log <-  read.table("/data02/Analysis/Projects/8_Aqueos_samples/4_CGE_Resfinder4.1/cge_resistome_species_ST.txt",sep = "\t", header = TRUE, quote="")
ARG_Cat <- read.csv("/data02/Analysis/Projects/8_Aqueos_samples/4_CGE_Resfinder4.1/ARG_Categories_withESBL.csv",sep = ",", header = TRUE, check.names = FALSE)
tr <- read.tree("/data02/Analysis/Projects/8_Aqueos_samples/9_SNP_matrix_redo_analysis/ECLOC_postGubbins.final_tree.tre")
head(cge_log)

sample_ARG_long <- as.data.frame(separate_rows(cge_log, resistance, sep=','))
head(sample_ARG_long)

sample_ARG_Cat_long <- as.data.frame(sample_ARG_long %>% right_join(ARG_Cat, by = 'resistance')%>% drop_na()) 
head(sample_ARG_Cat_long)

sample_ARG_Cat_long$SampleType[sample_ARG_Cat_long$SampleType == '11F'] <- 'Patient'
sample_ARG_Cat_long$SampleType[sample_ARG_Cat_long$SampleType == 'Env'] <- 'Environment'

# lets add color for the tip points- based on the sampleType
nb.cols <- length(unique(sample_ARG_Cat_long$SampleType))
#nb.cols
customPal <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols) # Set2, Accent, Paired
#class(customPal)
#customPal

sample_ARG_Cat_long$SampleTypeColor <- ifelse(sample_ARG_Cat_long$SampleType %in% "Env", 
                       customPal[1], 
                       ifelse(sample_ARG_Cat_long$SampleType %in% "11F", customPal[2], NA))

head(sample_ARG_Cat_long)
tail(sample_ARG_Cat_long)

# sample_ARG_Cat_long %>% 
#   pivot_wider(names_from=sample_name, values_from=resistance)


# dat1 <- sample_ARG_Cat_long %>% select(c("sample_name", "SampleType", "SampleTypeColor"))
# head(dat1)
# 
# # For the heatmap layer
# dat2 <- sample_ARG_Cat_long %>% select(c("sample_name", "resistance", "category"))
# head(dat2)

# Presence/Absence Matrix
mat <- acast(sample_ARG_Cat_long, sample_name~resistance, length)
head(mat)

thetree <- ggtree(tr=tr, open.angle=15, size=0.2) #+ #removed this in ggtree command: , aes(colour=Clade)
#thetree <- ggtree(tr=tr, open.angle=15, size=0.2) + ggplot2::ylim(NA, NA)
thetree

sampleType <- sample_ARG_Cat_long %>% select(sample_name,category,SampleType) %>%  unique()
head(sampleType)

sub_sampleType <- sampleType[sampleType$sample_name %in% thetree$data$label,]
class(thetree$data$label[!is.na(thetree$data$label)])
head(sub_sampleType)

# lets add color for the tip points- based on the sampleType
nb.cols <- length(unique(sample_ARG_Cat_long$SampleType))
nb.cols
#customPal <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols) # Set2, Accent, Paired
customPal <- colorRampPalette(c("red", "#24868EFF")) # Set2, Accent, Paired
#class(customPal)
customPal(nb.cols)

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
    values=customPal(nb.cols),
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=2,
                       override.aes=list(size=2,alpha=1))
  )
p1

distances <- as.data.frame(as.matrix(mat))
head(distances)
class(distances)
class(thetree$data$label)

#class(distances[na.omit(thetree$data$label),])
# The below line takes the tip labels first and subsets the distance matrix
# So, only samples in particular tree will be utilised 
# Then we omit the NA rows
# Subset the matrix again by removing the (genes) columns with all 0 and finally write into same matrix distances
distances <- subset(distances[na.omit(thetree$data$label),], select=colSums(distances[na.omit(thetree$data$label),]) > 0) 

# thetree_datasubset <- subset(thetree$data,thetree$data$isTip==TRUE)
# thetree_datasubset 
# 
# thetree_datasubset_sorted <- thetree_datasubset[order(thetree_datasubset$y, decreasing = TRUE),] 
# sorted_labels <- unlist(thetree_datasubset_sorted[,4])
# head(sorted_labels)
# class(sorted_labels)

# Sorting the presence/absence matrix for genes
distances.transpose <- t(distances)
distances.transpose
reordered_distances.transpose <- distances.transpose[do.call(order,as.data.frame(distances.transpose)),]
reordered_distances <- t(reordered_distances.transpose)

reordered_distances

EntSp_p2 <- gheatmap(p1,reordered_distances, offset=70000, width=2,
                  colnames_angle=90, colnames_offset_y = 0.1,font.size = 3) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "Enterobacter Species Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")

EntSp_p2

ggsave("Enterobacter Species_PGTree_Gene_PAMat.png",dpi = 600, width = 10, height = 8, units = "in")



KP_p2 <- gheatmap(p1,reordered_distances, offset=3000, width=2,
                  colnames_angle=90, colnames_offset_y = -0.2,font.size = 3) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "Klebsiella Pneumoniae Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")

KP_p2

ggsave("Klebsiella Pneumoniae_PGTree_Gene_PAMat.png",dpi = 600, width = 15, height = 12, units = "in")


SM_p2 <- gheatmap(p1,reordered_distances, offset=1, width=2,
                  colnames_angle=0, colnames_offset_y = 0.1,font.size = 3) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "Serratia Marcescens Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")

SM_p2

ggsave("Serratia Marcescens_PGTree_Gene_PAMat.png",dpi = 600, width = 15, height = 12, units = "in")


EC_p2 <- gheatmap(p1,reordered_distances, offset=45000, width=2,
               colnames_angle=90, colnames_offset_y = 0.1,font.size = 4) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "Escherichia Coli Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")

EC_p2

ggsave("Escherichia Coli_PGTree_Gene_PAMat.png",dpi = 600, width = 12, height = 8, units = "in")

ggsave("Serratia Marcescens_PGTree_Gene_PAMat.png",dpi = 600, width = 15, height = 12, units = "in")

#KP
Plot_PGT_SNPMat(
  newickFile = "SM_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_SM_full_triangle.csv",
  organism = "Serratia Marcescens",
  offset=3000,
  SNP_limits = c(0,100),
  legend_limits = c(0,100),
  wid = 15,
  hei = 10,
  colnames_offset_y = -0.2,
  font.size = 2
)


#SM
Plot_PGT_SNPMat(
  newickFile = "SM_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_SM_full_triangle.csv",
  organism = "Serratia Marcescens",
  offset=1,
  SNP_limits = c(0,1),
  legend_limits = c(0,1),
  wid = 15,
  hei = 12,
  colnames_offset_y = -0.1,
  font.size = 2.5
)

#SM
Plot_PGT_SNPMat(
  newickFile = "SM_postGubbins.final_tree.tre",
  meltfile = "postGubbins.filtered_polymorphic_sites_SM_full_triangle.csv",
  organism = "Serratia Marcescens",
  offset=1,
  SNP_limits = c(0,1),
  legend_limits = c(0,1),
  wid = 15,
  hei = 12,
  colnames_offset_y = -0.1,
  font.size = 2.5
)
