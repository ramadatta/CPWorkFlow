setwd("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST")

library(ggtree)
library(tidytree)
library(treeio)
library(lubridate)
library(ggplot2)
library(ggtree)
library(reshape2)
library(ggtreeExtra)
library(ggnewscale)
library(phytools)

# Load BEAST Tree
beast_tree <- read.beast(file = "postGubbins.filtered_polymorphic_sites_BEAST_withDates_treecombined_samplestate_10000_treeannotator.tree")
beast_tree

# Add Meta-data into dataframe
tipcategories = read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/BEAST/Combined_TTSH_NUH_SGH_SpecimenSite.csv", 
                           sep = ",",
                           header = TRUE)

head(tipcategories)

# Basic Tree
ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07") 

# Let's add root point
ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07") + 
  geom_treescale(x=2005, y=250, offset=2, fontsize = 3) + geom_rootpoint()


thetree <- ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07") + 
    geom_treescale(x=2005, y=250, offset=2, fontsize = 3) +  geom_rootpoint()

thetree  <- thetree + 
            theme_tree2() +  
            #ggplot2::xlim(1999, 2021) +
            geom_range(range='height_0.95_HPD', color='pink1', alpha=.20, size=2) + # Bars
            scale_x_continuous(breaks=c(2000,2001,2002,2003,2004,2005,2006,2007,
                                        2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020, 2021)) + 
  labs(x="", y="", title = "Pseudomonas aeruginosa Divergence Tree and Resistant Gene Presence/Absence Matrix")
thetree
# Add heatmap

meta_gheatmap <- tipcategories %>% select(BEAST_header,Institution, Specimen_Site2) %>% filter(BEAST_header!="#N/A") # let's not disturb original dataframe, create a new one

meta_gheatmap <- data.frame(column_to_rownames(meta_gheatmap, var = "BEAST_header"),check.names = FALSE) #assign ID to rownames for compatibility for gheatmap
head(meta_gheatmap)

meta_gheatmap <- meta_gheatmap %>% mutate(SampleType = case_when(
  grepl("ENV", Institution) ~ 'Environment',
  !grepl("ENV", Institution) ~ 'Patient'))

meta_gheatmap_sampType <- meta_gheatmap %>% select(SampleType)

p1 <- gheatmap(thetree, meta_gheatmap_sampType, offset=0.0001, width=0.01, low="darkorange1", high="blue", colnames = FALSE) +
  scale_fill_manual(values=c("Environment"="darkorange1","Patient"="blue"), name="Sample Type")

p1

meta_gheatmap_Inst <- meta_gheatmap %>% select(Institution)
head(meta_gheatmap_Inst)
p2 <- p1 + new_scale_fill()
p2

p3 <- gheatmap(p2, meta_gheatmap_Inst, offset=0.5, width=0.01, colnames = FALSE) + scale_fill_manual(values=c("NUH"="#E7B800", "SGH"="#FC4E07", "TTSH"="darkgreen", "TTSH_ENV"="darkblue"), name="Institution")
p3


# Adding resistant genes matrix to the right side of PGTree
cge_log <-  read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/Abricate/Abricate_270_samplesresults.csv",sep = ",", header = TRUE, quote="")

#view(cge_log)
#head(cge_log)
tail(cge_log)

cge_log %>% select(FILE,GENE) %>% filter(GENE=="blaNDM-1")

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

p4 <- p3 + new_scale_fill()
p4


####METHOD1 : to have heatmap (Label names of genes are not properly shown in linear chart - circular have to adjust a bit)

gheatmap(p4,reordered_distances, offset=0.5, width=1,
         colnames_angle=90, colnames_offset_y = -18,font.size = 3.5) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "BEAST Phylogenetic Tree and Resistant Gene Presence/Absence Matrix") + vexpand(.1, -1)
head(reordered_distances)
class(reordered_distances)
data.frame(reordered_distances)

gheatmap(p4,data.frame(factor(reordered_distances)), offset=0.5, width=1,
         colnames_angle=90, colnames_offset_y = -18,font.size = 3.5) +
   scale_color_manual(values=c("1"="#ffd662ff", "2"="#24868EFF"), name="Gene Presence/Absence") +
  # scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
  #                      breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
  labs(x="",
       y="",
       title = "BEAST Phylogenetic Tree and Resistant Gene Presence/Absence Matrix") + vexpand(.1, -1)


