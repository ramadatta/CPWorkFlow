
library(ggtree)
library(tidytree)
library(treeio)
library(lubridate)
library(ggplot2)
library(reshape2)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)
library(cowplot)
library(tidyverse)
library(data.table)
library(aplot)


########################################-----DIFFERENT ANALYSIS - PHYLOGENETIC TREE + TTSH SAMPLES -------#########################################
# Load BEAST Tree
setwd("/data02/Analysis/Projects/12_Ralstonia_Insidiosa_2022")

beast_tree <- read.tree(file = "postGubbins.filtered_polymorphic_sites.fasta.treefile")

#beast_tree <- read.tree(file = "postGubbins.final_tree.tre")
plot(beast_tree)

ggtree(beast_tree, ladderize=FALSE) + geom_treescale()
ggtree(beast_tree) + geom_treescale() + theme_tree2()
ggtree(beast_tree) + geom_treescale() + geom_tiplab()
ggtree(beast_tree) + geom_treescale() + xlim(0, 0.1) + 
  geom_tiplab(size=2,align=F,linesize=0.5) + 
  geom_rootpoint()

# Let us re-root the tree

bootTree_rooted <- ape::root(beast_tree,outgroup = c('T00892'),resolve.root=TRUE) 
ggtree(bootTree_rooted) + geom_treescale() + geom_tiplab() 
ggtree(bootTree_rooted) + geom_treescale() + geom_tiplab() + geom_rootpoint()
ggtree(bootTree_rooted, right = FALSE) + geom_treescale() + geom_tiplab() + geom_rootpoint()

thetree <- ggtree(bootTree_rooted, ladderize = TRUE,right = FALSE) + geom_treescale() + 
  geom_tiplab(size=4,align=F,linesize=0.5,color = "purple" ) + geom_rootpoint() # without xlim the labels are cutoff

thetree <- ggtree(bootTree_rooted, ladderize = TRUE,right = FALSE) + xlim(0,0.1) + geom_treescale() + 
  geom_tiplab(size=4,align=F,linesize=0.5,color = "purple" ) + geom_rootpoint() # without xlim the labels are cutoff

thetree

#thetree <- ggtree(beast_tree, open.angle=15, colour="#FF6F61FF", ladderize = TRUE, right = TRUE) + 
#  geom_treescale() + geom_tippoint(size=1) + geom_rootpoint() +  geom_tiplab(size=3,align=F,linesize=0.5,color = "purple" ) + xlim(0, 0.01)

# Let us add a Matrix in the right side

# Plot SNP Matrix along with Phylogenetic tree

snp_table <- read.table("/data02/Analysis/Projects/12_Ralstonia_Insidiosa_2022/postGubbins.filtered_polymorphic_sites_full_triangle.csv", sep = ",", header = TRUE)
head(snp_table)

colnames(snp_table) <- c("No", "Sample1", "Sample2", "SNPDiff")
complete_snp_table <- na.omit(snp_table)
complete_snp_table <- snp_table %>% filter(Sample1!="ADH0090_Duplicate") %>% filter(Sample2!="ADH0090_Duplicate")

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
reordered_distances <- distances[sorted_labels,]
head(reordered_distances)
class(reordered_distances)


# Crucial part
#SNPBySample_df_long <- melt(as.matrix(reordered_distances)
SNPBySample_df_long <- melt(as.table(as.matrix(reordered_distances)))
head(SNPBySample_df_long)
colnames(SNPBySample_df_long) <- c("Pair1","Pair2","Value")
head(SNPBySample_df_long)

ggplot(SNPBySample_df_long, aes (x = Pair1, y = Pair2, fill = Value, label = Value)) +
  geom_tile(aes(fill= Value)) +
  geom_text() +
 # scale_fill_manual(values=(brewer.pal(3, "RdYlBu")), na.value="white")
  scale_fill_gradient2(low = "#ffd662ff", mid= "#24868EFF", midpoint = 20, high = "tomato3", name='Gene Presence/Absence') 


snpmatrix_p5 <- ggplot(SNPBySample_df_long, aes(Pair1,Pair2,fill = Value, label = Value)) + 
  geom_tile(aes(fill=Value)) +
  geom_text() +
  #geom_text(aes(label=SNPDiff)) +
  scale_fill_gradient2(low = "#ffd662ff", mid= "#24868EFF", midpoint = 20, high = "tomato3", name='SNP Difference') +
  #scale_fill_gradient2(low = "#ffd662ff", mid = "#24868EFF", high="tomato3", name='Gene Presence/Absence') +
  # scale_fill_gradientn(limits = c(0, 40), colours=c("#ffd662ff", "#24868EFF"),
  #                      breaks=c(0, 1, 2, 10, 20, 30, 40), na.value="tomato3", name="SNP Difference") + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 8), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) 

snpmatrix_p5 

# But, I just want tree on left and heatmap on right

PT_SNPMatrix_Comb <- snpmatrix_p5 %>% insert_left(thetree,width=1)
PT_SNPMatrix_Comb 

pdf("Ralstonia_insidiosaPGTree_SNPMatrix.pdf", width=6, height=3)
PT_SNPMatrix_Comb

dev.off()

# Let's try gheatmap






