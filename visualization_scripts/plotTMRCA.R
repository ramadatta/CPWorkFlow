require(ggplot2)
require(ggtree)
library(treeio)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/14_Beast_205samples/Final_Batch9_1Billion_184samples")

# Read the data
tree <- read.beast("Run1_2_1Billion_treecombined_downsampled10000_treeannotater.tree")

# supply a most recent sampling date so you get the dates
# and add a scale bar
ggtree(tree, mrsd="2020-11-07") + theme_tree2() 

# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2020-11-07") + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5, size =3) + xlim(1990, 2021) 

# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2020-11-07") + 
  theme_tree2() + 
  geom_tiplab(align=FALSE, linesize=.5) + xlim(1990, 2021) 


##----PLotting Newick Tree From Gubbins Polymorphic sites for all 205 samples----#####
nwk <- system.file("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/6SNPDifftable_SNIPPY/", "gubbins.final_tree.tre", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree) 
p + geom_nodepoint(color="#b5e521", alpha=1/4, size=10)
p + geom_text(aes(label=label), size=3, color="purple", hjust=-0.3)
ggtree(tree) + geom_text(aes(label=label), size=3, color="purple", hjust=-0.3)


setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/14_Beast_205samples/Final_Batch9_1Billion_184samples/")
x <- read.beast("Run1_2_1Billion_treecombined_downsampled10000_treeannotater.tree")
x
cols <- scale_color(x, by="height")
pdf("myfile.pdf",width = 20, height = 24)
ggtree(x, right=TRUE, mrsd="2020-11-07") + theme_tree2() +
  geom_text(aes(x=max(x), label=label), size=3, hjust=-.3) +
  scale_x_continuous(breaks=c(2006,2007,2010,2015,2016,2017,2018,2019,2020)) +
  geom_segment(aes(xend=max(x)+.20, yend=y), linetype="dotted", size=.1, color=cols) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) 
dev.off()

viewClade(p, MRCA(p, "C02364|2020.82513661202", "345|2015.33424657534"))

##Original

# ggtree(x, right=FALSE, mrsd="2020-11-07", color=cols) + theme_tree2() +
#   geom_text(aes(x=max(x), label=label), size=3, color=cols, hjust=-.3) +
#   scale_x_continuous(breaks=c(2019,2020,2021), minor_breaks=seq(1)) +
#   geom_segment(aes(xend=max(x)+.20, yend=y), linetype="dotted", size=.1, color=cols) +
#   theme(panel.grid.major   = element_line(color="black", size=.2),
#         panel.grid.minor   = element_line(color="grey", size=.2),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank()) 

beast_file <- system.file("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/14_Beast_205samples/Run1_2_Combined_Annotated.trees", 
                          package="ggtree")
beast_tree <- read.beast(beast_file)
ggtree(beast_tree, mrsd="2020-11-07") + theme_tree2()



