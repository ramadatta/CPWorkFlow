setwd("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST")

library(ggtree)
library(tidytree)
library(treeio)
library(lubridate)
library(ggplot2)
library(ggtree)
library(reshape2)

# Load BEAST Tree
beast_tree <- read.beast(file = "postGubbins.filtered_polymorphic_sites_BEAST_withDates_treecombined_samplestate_10000_treeannotator.tree")
beast_tree

# Add Meta-data into dataframe
tipcategories = read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST/Combined_TTSH_NUH_SGH_SpecimenSite.csv", 
                           sep = ",",
                           header = TRUE)

head(tipcategories)

dd = as.data.frame(tipcategories)

head(dd)

# Basic Tree
thetree <- ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07") +  
#thetree <- ggtree(beast_tree, open.angle=180, layout="circular",mrsd = "2020-11-07")  + # If needed circular
  geom_range(range='height_0.95_HPD', color='red', alpha=.2, size=2) + # Bars
  geom_treescale(x=2005, y=250, offset=2, fontsize = 3) +
  geom_rootpoint()

thetree

# Add Details to tree
thetree %<+% dd + geom_tiplab(align=TRUE, linesize=.1, size =2, aes(col=Institution)) + 
 # geom_tippoint( aes(shape = Institution), color = "black", size = 1, show.legend = FALSE) + # To add tippoints
   ggplot2::xlim(1999, 2021) + theme_tree2() +
  # geom_text(aes(x=branch, label=height_0.95_HPD), size=2, vjust=-.3, color="firebrick") + # This will help find the values 2020.85245901639−11.178206=2009.6, 2020.85245901639−18.850128= 2002.006
   # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.1, x=branch), vjust=1,nudge_x = 0.2,size=3) + #Bootstrap values
  scale_color_manual(values=c("#E7B800", "#FC4E07", "darkgreen", "darkblue", "purple")) +
  scale_x_continuous(breaks=c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020, 2021))

# Storing the above tree in an object for adding other layers later
p1 <- thetree %<+% dd + geom_tiplab(align=TRUE, linesize=.1, size =2, aes(col=Institution)) + 
  # geom_tippoint( aes(shape = Institution), color = "black", size = 1, show.legend = FALSE) +
  ggplot2::xlim(1999, 2021) + theme_tree2() +
  # geom_text(aes(x=branch, label=height_0.95_HPD), size=2, vjust=-.3, color="firebrick") + # This will help find the values 2020.85245901639−11.178206=2009.6, 2020.85245901639−18.850128= 2002.006
  # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.1, x=branch), vjust=1,nudge_x = 0.2,size=3) + #Bootstrap values
  scale_color_manual(values=c("#E7B800", "#FC4E07", "darkgreen", "darkblue", "purple")) +
  scale_x_continuous(breaks=c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020, 2021)) +
  labs(x="",
       y="",
       title = "Pseudomonas aeruginosa Divergence Tree and Resistant Gene Presence/Absence Matrix")

# Save tree
setwd("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST")
#ggsave("pae_samples_multiple_hospitals.pdf",dpi = 300, width = 20, height = 25, units = "in")

# Adding resistant genes matrix to the right side of PGTree
cge_log <-  read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/Abricate/Abricate_270_samplesresults.csv",sep = ",", header = TRUE, quote="")

#view(cge_log)
#head(cge_log)
tail(cge_log)

# Presence/Absence Matrix
#mat <- acast(cge_log, BEAST_Header~GENE,length) #  this gives count matrix. But i want only pres/abs 
mat <- as.data.frame(with(cge_log, table(BEAST_Header, GENE)) > 0L) +0L
head(mat)
class(mat)
distances <- as.data.frame(as.matrix(mat))
head(distances)
class(distances)
class(thetree$data$label)

# The below line takes the tip labels first and subsets the distance matrix
# So, only samples in particular tree will be utilised 
# Then we omit the NA rows
# Subset the matrix again by removing the (genes) columns with all 0 and finally write into same matrix distances
distances <- subset(distances[na.omit(thetree$data$label),], select=colSums(distances[na.omit(thetree$data$label),]) > 0) 

# Sorting the presence/absence matrix for genes and restoring to BEAST_header matrix
distances.transpose <- t(distances)
distances.transpose
reordered_distances.transpose <- distances.transpose[do.call(order,as.data.frame(distances.transpose)),]
reordered_distances <- t(reordered_distances.transpose)
reordered_distances

# gheatmap(p1,reordered_distances, offset=15, width=2,
#                      colnames_angle=90, colnames_offset_y = 265,font.size = 2) +
#   scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
#                        breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
#   labs(x="",
#        y="",
#        title = "Enterobacter Species Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")
# 
# ggsave("Enterobacter Species_PGTree_Gene_PAMat.pdf",dpi = 600, width = 10, height = 8, units = "in")

library(aplot)
# p1
#geneBySample_df_long <- reshape2:::melt(as.matrix(mat))
geneBySample_df_long <- melt(reordered_distances)
colnames(geneBySample_df_long) <- c("BEAST_Header","Gene","Value")
head(geneBySample_df_long)

#geneBySample_df_long %>% filter(Value>=2)

arg_freq_df1 <- geneBySample_df_long %>% filter(Value == 1) 

arg_freq_df2 <- as.data.frame(table(arg_freq_df1$Gene))

colnames(arg_freq_df2) <- c("ARG","Frequency")

arg_freq <- ggplot(arg_freq_df2,aes(x=ARG,y=Frequency,fill=ARG)) +
  #geom_col(stat = count)
  geom_bar(stat = "identity") +
  geom_text(aes(label=Frequency), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), axis.text.x=element_blank(),axis.title.x = element_blank())

arg_freq 

p2 <- ggplot(geneBySample_df_long,aes(x=Gene, y=BEAST_Header)) + 
  geom_tile(aes(fill=Value)) + 
  #scale_fill_viridis_c() + 
  scale_fill_gradient(low = "#ffd662ff", high = "#24868EFF", name='Gene Presence/Absence', ) +
  theme_minimal() + xlab(NULL) + ylab(NULL) + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

p2 # This will contain heatmap of all the samples in the tree

# First lets put bar chart on top
# arg_freq_top_heatmap_down <- p2 %>% insert_top(arg_freq)
# p3 <- arg_freq_top_heatmap_down %>% insert_left(p1)
# p3

# But, I just want tree on left and heatmap on right

p4 <- p2 %>% insert_left(p1)
p4

ggsave("test.pdf",dpi = 600, width = 20, height = 18, units = "in")

#ggsave("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST/test.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

# Plot SNP Matrix along with Phylogenetic tree

snp_table <- read.table("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/Gubbins/postGubbins.filtered_polymorphic_sites_full_triangle_BEAST_headers.csv", sep = ",", header = TRUE)
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
reordered_distances <- distances[sorted_labels,]
head(reordered_distances)
class(reordered_distances)
p1

# Crucial part
#SNPBySample_df_long <- melt(as.matrix(reordered_distances)
SNPBySample_df_long <- melt(as.table(as.matrix(reordered_distances)))
head(SNPBySample_df_long)
colnames(SNPBySample_df_long) <- c("Pair1","Pair2","Value")
head(SNPBySample_df_long)

snpmatrix_p5 <- ggplot(SNPBySample_df_long, aes(Pair1,Pair2)) + 
  geom_tile(aes(fill=Value)) +
  #geom_text(aes(label=SNPDiff)) +
  scale_fill_gradient2(low = "#ffd662ff", mid= "#24868EFF", midpoint = 20, high = "tomato3", name='Gene Presence/Absence') +
  #scale_fill_gradient2(low = "#ffd662ff", mid = "#24868EFF", high="tomato3", name='Gene Presence/Absence') +
  # scale_fill_gradientn(limits = c(0, 40), colours=c("#ffd662ff", "#24868EFF"),
  #                      breaks=c(0, 1, 2, 10, 20, 30, 40), na.value="tomato3", name="SNP Difference") + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 2), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

snpmatrix_p5
# But, I just want tree on left and heatmap on right

PT_SNPMatrix_Comb <- snpmatrix_p5 %>% insert_left(p1)
PT_SNPMatrix_Comb
