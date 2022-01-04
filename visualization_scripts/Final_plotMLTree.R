
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


########################################-----DIFFERENT ANALYSIS - PHYLOGENETIC TREE + TTSH SAMPLES -------#########################################
# Load BEAST Tree
setwd("/data02/Analysis/Projects/6_Paeruginosa_ST308_Singapore_vs_Global/Roary/CoreGenome")

beast_tree <- read.tree(file = "postGubbins.filtered_polymorphic_sites_wo2_P_PA.oneline.fasta.treefile")
plot(beast_tree)

ggtree(beast_tree) + geom_treescale()
ggtree(beast_tree) + geom_treescale() + theme_tree2()
ggtree(beast_tree) + geom_treescale() + geom_tiplab()
ggtree(beast_tree) + geom_treescale() + xlim(0, 10) + geom_tiplab(size=2,align=F,linesize=0.5)

# Let us re-root the tree

bootTree_rooted <- ape::root(beast_tree, outgroup = c('PAO1_reference'), edgelabel = TRUE) 
ggtree(bootTree_rooted) + geom_treescale() + geom_tiplab() 
ggtree(bootTree_rooted) + geom_treescale() + geom_tiplab() + geom_rootpoint()
ggtree(bootTree_rooted, right = FALSE) + geom_treescale() + geom_tiplab() + geom_rootpoint()
ggtree(bootTree_rooted, right = TRUE) + xlim(0, 0.025) + geom_treescale() + geom_tiplab(size=2,align=F,linesize=0.5) + geom_rootpoint()

# Since the PAO1 reference is too far from rest of samples, for time being it is pruned from the tree object

to_drop <- c("PAO1_reference")
thetree_wo_ref <- as.phylo(bootTree_rooted)
nhx_reduced <- drop.tip(thetree_wo_ref, to_drop)

# Default layout
thetree <- ggtree(nhx_reduced) + geom_treescale() + geom_tiplab(size=2,align=F,linesize=0.5) + geom_rootpoint()

# layout = roundrect
thetree <- ggtree(nhx_reduced,layout='roundrect', ladderize = TRUE, right = TRUE) + 
            geom_treescale() + geom_tiplab(size=2,align=F,linesize=0.5) + geom_rootpoint()

# layout = circular
thetree <- ggtree(nhx_reduced,layout='circular', ladderize = TRUE, right = TRUE) + 
  xlim(0, 0.008) + geom_treescale() + geom_tiplab(size=2,align=F,linesize=0.5) + geom_rootpoint()

# layout = roundrect
thetree <- ggtree(nhx_reduced,layout='roundrect', colour="purple", ladderize = TRUE, right = TRUE) + 
   geom_treescale() + geom_tiplab(size=2,align=F,linesize=0.5) + geom_rootpoint()

thetree

# Lets add metadata - institution

# Add Meta-data into dataframe
meta = read.table("/data02/Analysis/Projects/6_Paeruginosa_ST308_Singapore_vs_Global/Roary/CoreGenome/metadata.txt", 
                         sep = "\t", header = TRUE)

head(meta)

thetree_byCountry <- thetree %<+% meta + geom_tiplab(align=F, linesize=.1, size =2, aes(col=Country)) 
thetree_byLocation <- thetree %<+% meta + geom_tiplab(align=F, linesize=.1, size =2, aes(col=Location)) 

# Let's combine the plots using cowplot package function
plot_grid(thetree_byCountry, thetree_byLocation, labels = c('By Country', 'Local vs International'), label_size = 12)

# Add heatmap


meta_gheatmap <- meta # let's not disturb original dataframe, create a new one

meta_gheatmap <- data.frame(column_to_rownames(meta, var = "Sample"),check.names = FALSE) #assign ID to rownames for compatibility for gheatmap
head(meta_gheatmap)

meta_gheatmap_loc <- meta_gheatmap %>% select(Location)

#-----------> Adding the Location - Is the isolate local or international

gheatmap(thetree_byLocation, meta_gheatmap_loc, offset = 0, width = 1,
         colnames_position = "top",
         colnames_offset_y = 1) +
  scale_fill_viridis_d(option = "D", name = "Clade") + coord_cartesian(clip = "off")

# Let's beautify the plot a bit more

# With tip labels
p10 <- ggtree(nhx_reduced,layout='roundrect', colour="purple", ladderize = TRUE, right = TRUE) + 
  geom_treescale() + geom_tiplab(size=2,align=T,linesize=0.005,offset=0.0002) + geom_rootpoint()

# add heatmap
p10 <- gheatmap(p10, meta_gheatmap_loc, offset=0.0001, width=0.01, low="darkorange1", high="blue", colnames = FALSE, font.size=2, color="black") +
  scale_fill_manual(values=c("darkorange1","blue"),name="Location")

#p10 

# Without tip labels
p11 <- ggtree(nhx_reduced,layout='roundrect', colour="purple", ladderize = TRUE, right = TRUE) + 
  geom_treescale() + geom_rootpoint()

# add heatmap
p11 <- gheatmap(p11, meta_gheatmap_loc, offset=0.0001, width=0.01, low="darkorange1", high="blue", colnames = FALSE) +
  scale_fill_manual(values=c("Global"="darkorange1","Local"="blue"), name="Location")

p11

# Let's combine the plots using cowplot package function
plot_grid(p10, p11, labels = c('With Tip Labels', 'Without Tip Labels'), label_size = 12)

#-----------> Adding the Country 

# I want to add on another layer called Country, so I will take previous pretty p11 plot
meta_gheatmap_ctry <- meta_gheatmap %>% select(Country)

p12 <- p11 + new_scale_fill()


# add heatmap
p13 <- gheatmap(p12, meta_gheatmap_ctry, offset=0.0003, width=0.01, colnames = FALSE) +
  scale_fill_manual(values=c("Germany"="red","Spain"="skyblue","France"="Green", 
                             "Singapore_NUH"="magenta", "Singapore_TTSH"="orange", "Singapore_TTSH_ENV"="tomato3",
                             "Singapore_SGH"="yellow"), 
                    name="Country") 
p13

p14 <- gheatmap(p12, meta_gheatmap_ctry, offset=0.0003, width=0.01, colnames = FALSE) +
  scale_fill_manual(values=c("Germany"="red","Spain"="skyblue","France"="Green", 
                             "Singapore_NUH"="magenta", "Singapore_TTSH"="orange", "Singapore_TTSH_ENV"="tomato3",
                             "Singapore_SGH"="yellow"), 
                    breaks=c("Germany","Singapore_NUH","Spain","Singapore_TTSH","Singapore_TTSH_ENV","France","Singapore_SGH"),name="Country") 

p14

# Let's combine the plots using cowplot package function
plot_grid(p13, p14, labels = c('With Breaks', 'Without breaks'), label_size = 12)

#-------> Spread the countries across multiple layers

meta_gh_color <- meta
head(meta_gh_color)
# Let us assign a color value for each country
meta_gh_color <- meta_gh_color %>%
  mutate(Color = case_when(
    grepl("Germany", Country) ~ 'red',
    grepl("Spain", Country) ~ 'skyblue',
    grepl("France", Country) ~ 'Green',
    grepl("Singapore_NUH", Country) ~ 'magenta',
    grepl("^Singapore_TTSH$", Country) ~ 'orange',
    grepl("Singapore_TTSH_ENV", Country) ~ 'red',
    grepl("Germany", Country) ~ 'tomato3',
    grepl("Singapore_SGH", Country) ~ 'yellow'))

meta_gh_color <- meta_gh_color %>%
  mutate(Color = case_when(
    grepl("Global", Location) ~ 'red',
    grepl("Local", Location) ~ 'blue',
    ))

view(meta_gh_color)
head(meta_gh_color)
meta_gh_color %>% filter(Country=="Singapore_TTSH_ENV")

p15 <- p12 + geom_fruit(data=meta_gh_color,geom=geom_tile,
                 mapping = aes(y=Sample, x=Color, fill=Country),
                 offset=0.024,pwidth=0.07) #+
   # scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
   #                           "#800000", "#006400","#800080","#696969"))

p15

# Adding axis.params - the labels are not fitting in properly
p16 <- p12 + geom_fruit(data=meta_gh_color,geom=geom_tile,
                        mapping = aes(y=Sample, x=Country, fill=Country),
                        offset=0.024,pwidth=0.07,
                      axis.params = list(axis="x", text.size=2, vjust=1, line.size=0, hjust=1, text.angle=-90)) 
  
  
p16

p16.1 <- p12 + geom_fruit(data=meta_gh_color,geom=geom_tile,
                 mapping = aes(y=Sample, x=Country, fill=Country),
                 offset=0.024,pwidth=0.07,
                 axis.params = list(axis="x", text.size=2, vjust=0, line.size=0, hjust=0, text.angle=-90)) + vexpand(.1, -1)

p16.1

# Let's combine the plots using cowplot package function
plot_grid(p16, p16.1, labels = c('Without vexpand', 'With vexpand'), label_size = 12)

# Let's try circularizing the plot
# layout = circular
ggtree(nhx_reduced,layout="fan", open.angle=0)
ggtree(nhx_reduced,layout="fan", open.angle=10)
ggtree(nhx_reduced,layout="fan", open.angle=30)
ggtree(nhx_reduced,layout="fan", open.angle=60)
ggtree(nhx_reduced,layout="fan", open.angle=90)

thecirc_tree <- ggtree(nhx_reduced,layout="fan", open.angle=40)  

p17 <- gheatmap(thecirc_tree, meta_gheatmap_loc, offset=0.0001, width=0.1, low="darkorange1", high="blue", colnames = FALSE) +
  scale_fill_manual(values=c("Global"="darkorange1","Local"="blue"), name="Location")

p17
p18 <- p17 + new_scale_fill()

p19 <- p18 + geom_fruit(data=meta_gh_color,geom=geom_tile,
                 mapping = aes(y=Sample, x=Country, fill=Country),
                 offset=0.24,pwidth=0.7,
                 axis.params = list(axis="x", text.size=3, vjust=1, line.size=0, hjust=1, text.angle=-90)) 
p19

# Labels does not nice with circular tree 
p20 <- p18 + geom_fruit(data=meta_gh_color,geom=geom_tile,
                 mapping = aes(y=Sample, x=Country, fill=Country),
                 offset=0.16,pwidth=0.7) 


p20.rotated <- rotate_tree(p20, -125)
p20.rotated
# Let's combine the plots using cowplot package function
plot_grid(p20, p20.rotated, labels = c('Raw Gheatmap', 'Rotated Tree'), label_size = 12)

# For Manuscript

# Figure 1: Non-dated ML tree 
head(meta_gh_color)
view(meta_gh_color)

ggtree(nhx_reduced,colour="#FF6600") %<+% meta_gh_color + geom_tippoint(aes(color=I(Color)))

ggtree(nhx_reduced,colour="#FF6600") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("red","blue"), name="Location")

ggtree(nhx_reduced,colour="black") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#949398FF","#F4DF4EFF"), name="Location")

ggtree(nhx_reduced,colour="black") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#FC766AFF","#5B84B1FF"), name="Location")


ggtree(nhx_reduced,colour="black") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#ADEFD1FF"), name="Location")

ggtree(nhx_reduced,colour="#5B84B1FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#ADEFD1FF"), name="Location")

ggtree(nhx_reduced,colour="#FC766AFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#ADEFD1FF"), name="Location")

ggtree(nhx_reduced,colour="#FC766AFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#ADEFD1FF"), name="Location")

ggtree(nhx_reduced,colour="#FC766AFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#2BAE66FF"), name="Location")

ggtree(nhx_reduced,colour="#FC766AFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#ADEFD1FF"), name="Location")

ggtree(nhx_reduced,colour="#5B84B1FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#FC766AFF"), name="Location")

ggtree(nhx_reduced,colour="#5B84B1FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#990011FF"), name="Location")

ggtree(nhx_reduced,colour="#FC766AFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00203FFF","#990011FF"), name="Location")

ggtree(nhx_reduced,colour="#FFD662FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00539CFF","#990011FF"), name="Location")

ggtree(nhx_reduced,colour="#990011FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00539CFF","#F95700FF"), name="Location")

ggtree(nhx_reduced,colour="#00539CFF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#5B84B1FF","#F95700FF"), name="Location")

ggtree(nhx_reduced,colour="#F95700FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#00539CFF","#990011FF"), name="Location")

ggtree(nhx_reduced,colour="#011936FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#FEAE51FF","#ED254EFF"), name="Location")


ggtree(nhx_reduced,colour="#FF6F61FF") %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  scale_color_manual(values=c("#011936FF","#ED254EFF"), name="Location")
# 
# 
# 
# ggtree(nhx_reduced,colour="#3E282BFF", ladderize = TRUE, right = TRUE) %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
#   geom_treescale() +
#   geom_rootpoint() +
#   scale_color_manual(values=c("#008C76FF","#D34F73FF"), name="Location")

ggtree(nhx_reduced,colour="#FF6F61FF", ladderize = TRUE, right = TRUE) %<+% meta_gh_color + geom_tippoint(aes(color=Location)) +
  geom_treescale() +
  geom_rootpoint() +
  scale_color_manual(values=c("#011936FF","#ED254EFF"), name="Location")


# Figure 2: Dated BEAST tree 

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


thetree <- ggtree(beast_tree, open.angle=15, mrsd = "2020-11-07",colour="#FF6F61FF", ladderize = TRUE, right = TRUE) + 
  geom_treescale(x=2007, y=50, offset=2, fontsize = 3) + geom_tippoint(size=1) + geom_rootpoint()

thetree 

thetree  <- thetree + 
  theme_tree2() +  
  #ggplot2::xlim(1999, 2021) +
  geom_range(range='height_0.95_HPD', color='skyblue', alpha=.30, size=2) + # Bars
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

# p1 <- gheatmap(thetree, meta_gheatmap_sampType, offset=0.0001, width=0.01, low="darkorange1", high="blue", colnames = FALSE) +
#   scale_fill_manual(values=c("Environment"="darkorange1","Patient"="blue"), name="Sample Type")
# 
# 
# p1 <- gheatmap(thetree, meta_gheatmap_sampType, offset=0.0001, width=0.01, low="#F9DC5CFF", high="#ED254EFF", colnames = FALSE) +
#   scale_fill_manual(values=c("Environment"="#F9DC5CFF","Patient"="#ED254EFF"), name="Sample Type")

p1 <- gheatmap(thetree, meta_gheatmap_sampType, offset=0.0001, width=0.01, low="#F9DC5CFF", high="#ED254EFF", colnames = TRUE,colnames_angle=90, colnames_offset_y = -15,font.size = 3.5) +
  scale_fill_manual(values=c("Environment"="#ED254EFF","Patient"="#6DAC4FFF"), name="Sample Type")
p1

p2 <- p1 + new_scale_fill()
p2

head(meta_gheatmap)

meta_gheatmap <- meta_gheatmap %>% mutate(Hospital = case_when(
  grepl("NUH", Institution) ~ 'A',
  grepl("SGH", Institution) ~ 'B',
  grepl("TTSH", Institution) ~ 'C'))

# meta_gheatmap_Inst <- meta_gheatmap %>% select(Institution)
# p3 <- gheatmap(p2, meta_gheatmap_Inst, offset=0.5, width=0.01, colnames = FALSE) + 
#   scale_fill_manual(values=c("NUH"="#E7B800", "SGH"="#FC4E07", "TTSH"="darkgreen", "TTSH_ENV"="darkblue"), name="Institution")
# 
# p3 <- gheatmap(p2, meta_gheatmap_Inst, offset=0.5, width=0.01, colnames = FALSE) + 
#   scale_fill_manual(values=c("NUH"="#F93822FF", "SGH"="#006B38FF", "TTSH"="#FFD653FF", "TTSH_ENV"="black"), name="Institution")
# p3


# For Anonymized Institution A - NUH, B- SGH, C-TTSH
#head(meta_gheatmap_Inst)

meta_gheatmap_AnonInst <- meta_gheatmap %>% select(Hospital)

# variation
# p3 <- gheatmap(p2, meta_gheatmap_AnonInst, offset=0.5, width=0.01, colnames = FALSE) + 
#   scale_fill_manual(values=c("NUH"="#E7B800", "SGH"="#FC4E07", "TTSH"="darkgreen", "TTSH_ENV"="darkblue"), name="Institution")

p3 <- gheatmap(p2, meta_gheatmap_AnonInst, offset=0.5, width=0.01, colnames = TRUE, colnames_angle=90, colnames_offset_y = -10,font.size = 3.5) + 
  scale_fill_manual(values=c("A"="#F93822FF", "B"="#006B38FF", "C"="#FFD653FF"), name="Hospital")
p3


# Adding resistant genes matrix to the right side of PGTree
cge_log <-  read.csv("/data02/Analysis/Projects/6_Paeruginosa_199_TTSH_31_NUH_40_SGH/Abricate/Abricate_270_samplesresults.csv",sep = ",", header = TRUE, quote="")

#view(cge_log)
#head(cge_log)
tail(cge_log)

#cge_log %>% select(FILE,GENE) %>% filter(GENE=="blaNDM-1")

# rename fosA in cge_log 
library(stringr)
cge_log$GENE = str_replace(cge_log$GENE,"fosA-354827590","fosA")


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
view(head(reordered_distances))
class(reordered_distances)

names(reordered_distances)[names(reordered_distances) == 'fosA-354827590'] <- 'fosA'
head(reordered_distances)
#nrow(reordered_distances)
p4 <- p3 + new_scale_fill()
p4


####METHOD1 : to have heatmap (Label names of genes are not properly shown in linear chart - circular have to adjust a bit)

# gheatmap(p4,reordered_distances, offset=1, width=0.5,
#          colnames_angle=90, colnames_offset_y = -18,font.size = 3.5) +
#   scale_fill_gradientn(limits = c(0,1), colours=c("#ffd662ff", "#24868EFF"),
#                        breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Presence/Absence") +
#   labs(x="",
#        y="",
#        title = "BEAST Phylogenetic Tree and Resistant Gene Presence/Absence Matrix") + vexpand(.1, -1)
gheatmap(p4,reordered_distances, offset=1, width=0.4,
         colnames_angle=90, colnames_offset_y = -18,font.size = 3.5,color=NA) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#FBEEE6", "#68B2A0"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Resistant Gene Pres/Abs") +
  labs(x="",
       y="",
       title = "Phylogenetic tree with divergence date estimates of the 261 ST308 Pseudomonas aeruginosa samples and Resistant Gene Presence/Absence Matrix") + vexpand(.1, -1)

localtree_ggplot_heatmap <- gheatmap(p4,reordered_distances, offset=1, width=0.4,
         colnames_angle=90, colnames_offset_y = -18,font.size = 3.5,color=NA) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#FBEEE6", "#68B2A0"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Resistant Gene Pres/Abs") +
  labs(x="",
       y="",
       title = "Phylogenetic tree with divergence date estimates of the 261 ST308 Pseudomonas aeruginosa samples and Resistant Gene Presence/Absence Matrix") + vexpand(.1, -1)


pdf("Three_Hospitals_SNIPPY_ST308_PAE_PGTree_ARG_PAMatrix_28122021.pdf", width=15, height=10)
localtree_ggplot_heatmap
dev.off()

# library(Cairo)
# cairo_pdf("Three_Hospitals_SNIPPY_ST308_PAE_PGTree_ARG_PAMatrix_28122021_cairo.pdf", width=15, height=10, family = "Times")
# localtree_ggplot_heatmap
# dev.off()


# # Figure 3: Compare closed genomes between the local and global  .. in progress

# # Figure 4: Compare closed genomes between the local and global  .. in progress

thetree

# Add PAN genome heatmap


# Adding PAN Genome from roary output to the right side of PGTree
gene_presence <- read.table("/data02/Analysis/Projects/6_Paeruginosa_ST308_Singapore_vs_Global/Roary/PanGenome/gene_presence_absence.header_renamed.Rtab", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
row.names(gene_presence) <- gene_presence[,1]
gene_presence <- gene_presence[,-1]
gene_presence = t(gene_presence)
view(head(gene_presence))


gheatmap(thetree,gene_presence, offset=1, width=1, colnames = FALSE,
         colnames_angle=90, colnames_offset_y = -18,font.size = 3.5,color=NA) +
  scale_fill_gradientn(limits = c(0,1), colours=c("#FBEEE6", "#68B2A0"),
                       breaks=c(0,1), labels=format(c(0,1)), na.value="tomato3", name="Gene Pres/Abs") +
  labs(x="",
       y="",
       title = "Phylogenetic tree with divergence date estimates of the 261 ST308 Pseudomonas aeruginosa samples and Pan genome Presence/Absence Matrix") + vexpand(.1, -1)



