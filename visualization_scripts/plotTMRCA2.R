require(ggplot2)
require(ggtree)
library(treeio)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/14_Beast_205samples/Batch1_300Million_205samples/")

# Read the data
beast_tree <- read.beast("Run1_2_Combined_Annotated.trees")
file.exists("core.aln_withDates.ExcludedSamples_withoutDates_oneline.fasta")
fasta <- "core.aln_withDates.ExcludedSamples_withoutDates_oneline.fasta"
fasta
msaplot(ggtree(beast_tree), fasta) 
msaplot(ggtree(beast_tree), fasta,window=c(150, 200)) + coord_polar(theta='y')

ggtree(beast_tree, ladderize=TRUE)

p <- ggtree(beast_tree, mrsd="2020-11-07") + geom_treescale(x=2008, y=1)
p <- p + geom_tiplab(size=3)
p

q <- ggtree(beast_tree, mrsd="2020-01-01") + geom_tiplab(size=3, align=TRUE) + theme_tree2()
qq <- (q + scale_y_continuous(expand=c(0, 0.3))) %>%
    scale_x_ggtree()
qq + theme(legend.position="right")
qq
# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2020-11-07") + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5, size =3) + xlim(1990, 2021) 


