setwd("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST")

library(ggtree)
library(tidytree)
library(treeio)
library(lubridate)
library(ggplot2)
library(ggtree)
library(reshape2)

#-----> Go to the bottom of the script

# beast_tree <- read.beast(file = "3IndependentRuns_1Billion_logcombined_downsampled10000_treeannotater.tree")
# ggtree(beast_tree, mrsd="2020-11-07") + theme_tree2()
# 
# ## separate the tree by host species
# tip <- get.tree(beast_tree)$tip.label
# tip
# beast_tree <- groupOTU(beast_tree, tip[grep("C", tip)], 
#                        group_name = "host")
# # 1
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#   theme_classic() + theme(legend.position='right') +
#   scale_color_manual(values=c("blue", "red"), 
#                      labels=c("Institution A", "Institution B")) +
#   ylab("Bla")
# 
# #2 
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#   geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,offset = 2) +
#   theme_classic() + theme(legend.position='right') +
#   scale_color_manual(values=c("steelblue", "tomato3"), 
#                      labels=c("Institution A", "Institution B")) +
#   ylab("Bla")
# 
# #3
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#   geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,offset = 2) +
#   geom_range("length_0.95_HPD", color='orange', size=1, alpha=.5) +
#   theme_classic() + theme(legend.position='right') +
#   scale_color_manual(values=c("steelblue", "tomato3"), 
#                      labels=c("Institution A", "Institution B")) +
#   ylab("Samples")
# 
# #4
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#   geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,offset = 2) +
#   geom_range("length_0.95_HPD", color='orange', size=1, alpha=.5) +
#   theme_classic() + theme(legend.position='right') +
#   scale_color_manual(values=c("steelblue", "tomato3"), 
#                      labels=c("Institution A", "Institution B")) +
# geom_text2(aes(label=round(as.numeric(posterior), 2), 
#                subset=as.numeric(posterior)> 0.9, 
#                x=branch), vjust=0) 
# 
# #5
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#   geom_tippoint(aes(color=host), size=1) +
#   geom_tiplab(aes(color=host), align=TRUE, linetype='dashed', size=1.5,linesize=.3, offset = 2) +
#  # geom_range("length_0.95_HPD", color='orange', size=1, alpha=.5) +
#   theme_tree2()+ theme(legend.position='right') +
#   scale_color_manual(values=c("steelblue", "tomato3"), 
#                      labels=c("Institution B", "Institution A")) +
#  # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.9, x=branch), vjust=0) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
# 
# #6 
# #get the node number
# ggtree(beast_tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
# 
# MRCA(beast_tree, tip=c('345|2015.33424657534', 'C02434|2020.39617486339'))
# 
# ggtree(beast_tree, aes(color=host), mrsd="2020-11-07") + 
#  #geom_treescale(x=0, y=45, width=1, color='red') +
#   geom_rootpoint() + 
#   geom_rootedge(rootedge = 5) +
#   geom_tippoint(aes(color=host), size=1) +
#   geom_tiplab(aes(color=host), align=TRUE, linetype='dashed', size=3,linesize=.1, offset = 1) +
#    geom_range("length_0.95_HPD", color='orange', size=1, alpha=.5) +
#   theme_tree2()+ theme(legend.position='right') +
#   #scale_x_continuous( limit  = c(2005, 2021)) +
#   # scale_x_continuous(breaks = c(2005:2021), 
#   #                    labels = factor(2005:2021), 
#   #                    limits = c(2005,2021))
#   scale_x_continuous(expand=expansion(0.2)) + # make more room for the labels
#   #scale_y_tree() +
#   scale_color_manual(values=c("steelblue", "tomato3"), 
#                      labels=c("Institution B", "Institution A")) +
#   # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.9, x=branch), vjust=0) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
# 
# beast_tree
# #7
# # from https://github.com/yaotli/h5nx_mrca/blob/master/figure.md

#8 #===> Worked

beast_tree <- read.beast(file = "postGubbins.filtered_polymorphic_sites_BEAST_withDates_treecombined_samplestate_10000_treeannotator.tree")
beast_tree

## separate the tree by host species
# tip <- get.tree(beast_tree)$tip.label
# tip
# tip[grep("C", tip)]
# beast_tree <- groupOTU(beast_tree, tip[grep("C", tip)], 
#                        group_name = "Institute")
# 
# ggtree(beast_tree, mrsd = "2020-11-07")  +
#   geom_treescale(x=2000, y=5, offset=2, fontsize = 3) +
#   geom_rootpoint() + 
#   geom_tiplab(aes(color=Institute), align=TRUE, linetype='dashed', size=1,linesize=.3, offset = 1) + ggplot2::xlim(2000, 2025) +
#  # geom_range(range='height_0.95_HPD', color='orange', alpha=.6, size=1) +
#   theme_tree2() + scale_color_manual(values=c("#1B9E77", "#D95F02"), 
#                                      labels=c("Institution B", "Institution A")) +
#   geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.9, x=branch), vjust=1,nudge_x = 0.2,size=3) + #https://www.researchgate.net/post/How_do_I_display_bootstrap_values_using_BEAST_and_Figtree_on_my_phylogenetic_tree
#   geom_tippoint( aes(shape = Institute), color = "black", size = 1, show.legend = FALSE) 

ggsave("3IndependentRuns_1Billion_logcombined_downsampled10000_treeannotater.tree.png",dpi = 300, width = 8, height = 8, units = "in")


tipcategories = read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST/Combined_TTSH_NUH_SGH_SpecimenSite.csv", 
                           sep = ",",
                           header = TRUE)

head(tipcategories)

dd = as.data.frame(tipcategories)

head(dd)


thetree <- ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07") +  
#thetree <- ggtree(beast_tree, open.angle=180, layout="circular",mrsd = "2020-11-07")  +
  geom_range(range='height_0.95_HPD', color='red', alpha=.2, size=2) + # Bars

  geom_treescale(x=2005, y=250, offset=2, fontsize = 3) +
  geom_rootpoint()

#thetree

# Add Details to tree
thetree %<+% dd + geom_tiplab(align=TRUE, linesize=.1, size =2, aes(col=Institution)) + 
 # geom_tippoint( aes(shape = Institution), color = "black", size = 1, show.legend = FALSE) +
   ggplot2::xlim(1999, 2021) + theme_tree2() +
  # geom_text(aes(x=branch, label=height_0.95_HPD), size=2, vjust=-.3, color="firebrick") + # This will help find the values 2020.85245901639−11.178206=2009.6, 2020.85245901639−18.850128= 2002.006
   # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.1, x=branch), vjust=1,nudge_x = 0.2,size=3) + #Bootstrap values
  scale_color_manual(values=c("#E7B800", "#FC4E07", "darkgreen", "darkblue", "purple")) +
  scale_x_continuous(breaks=c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020, 2021))

# For adding other layers later
p1 <- thetree %<+% dd + geom_tiplab(align=TRUE, linesize=.1, size =2, aes(col=Institution)) + 
  # geom_tippoint( aes(shape = Institution), color = "black", size = 1, show.legend = FALSE) +
  ggplot2::xlim(1999, 2021) + theme_tree2() +
  # geom_text(aes(x=branch, label=height_0.95_HPD), size=2, vjust=-.3, color="firebrick") + # This will help find the values 2020.85245901639−11.178206=2009.6, 2020.85245901639−18.850128= 2002.006
  # geom_text2(aes(label=round(as.numeric(posterior), 2), subset=as.numeric(posterior)> 0.1, x=branch), vjust=1,nudge_x = 0.2,size=3) + #Bootstrap values
  scale_color_manual(values=c("#E7B800", "#FC4E07", "darkgreen", "darkblue", "purple")) +
  scale_x_continuous(breaks=c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020, 2021))


# Save tree
setwd("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST")
#ggsave("pae_samples_multiple_hospitals.pdf",dpi = 300, width = 20, height = 25, units = "in")

# Adding resistant genes matrix to the right side of PGTree

cge_log <-  read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/Abricate/Abricate_270_samplesresults.csv",sep = ",", header = TRUE, quote="")

#view(cge_log)
head(cge_log)
tail(cge_log)

# Presence/Absence Matrix
#mat <- acast(cge_log, BEAST_Header~GENE,length) #  this gives count matrix. But i want only pres/abs 
mat <- as.data.frame(with(cge_log, table(BEAST_Header, GENE)) > 0L) +0L
head(mat)

class(mat)
sub_sampleType <- cge_log[cge_log$BEAST_Header %in% thetree$data$label,]
sub_sampleType
class(thetree$data$label[!is.na(thetree$data$label)])
head(sub_sampleType)

distances <- as.data.frame(as.matrix(mat))
head(distances)
class(distances)
class(thetree$data$label)

# The below line takes the tip labels first and subsets the distance matrix
# So, only samples in particular tree will be utilised 
# Then we omit the NA rows
# Subset the matrix again by removing the (genes) columns with all 0 and finally write into same matrix distances
distances <- subset(distances[na.omit(thetree$data$label),], select=colSums(distances[na.omit(thetree$data$label),]) > 0) 

# Sorting the presence/absence matrix for genes
distances.transpose <- t(distances)
distances.transpose
reordered_distances.transpose <- distances.transpose[do.call(order,as.data.frame(distances.transpose)),]
reordered_distances <- t(reordered_distances.transpose)

reordered_distances

gheatmap(p1,reordered_distances, offset=15, width=2,
                     colnames_angle=90, colnames_offset_y = 265,font.size = 2) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "Enterobacter Species Gubbins Phylogenetic Tree and Resistant Gene Presence/Absence Matrix")

ggsave("Enterobacter Species_PGTree_Gene_PAMat.pdf",dpi = 600, width = 10, height = 8, units = "in")

library(aplot)
p1
class(mat)
head(mat)
#geneBySample_df_long <- reshape2:::melt(as.matrix(mat))
geneBySample_df_long <- melt(reordered_distances)
colnames(geneBySample_df_long) <- c("BEAST_Header","Gene","Value")
head(geneBySample_df_long)
#geneBySample_df_long %>% filter(Value>=2)

p2 <- ggplot(geneBySample_df_long,aes(x=Gene, y=BEAST_Header)) + 
  geom_tile(aes(fill=Value)) + scale_fill_viridis_c() + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

p2

p3 <- p2 %>% insert_left(p1)
p3

#ggsave("test.pdf",dpi = 600, width = 20, height = 18, units = "in")
#ggsave("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST/test.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

