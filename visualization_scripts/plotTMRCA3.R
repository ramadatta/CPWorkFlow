setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/14_Beast_205samples/Final_Batch9_1Billion_184samples/")
x <- read.beast("Run1_2_1Billion_treecombined_downsampled10000_treeannotater.tree")
x
cols <- scale_color(x, by="CAheight_mean", low="#0072B2", high="#D55E00",interval=seq(0, 1.5, length.out=100))
cols
ggtree(x, right=TRUE, mrsd="2005-04-02", color=cols) + theme_tree2() +
  geom_text(aes(x=max(x), label=label), size=1, color=cols, hjust=-.3) +
  scale_x_continuous(breaks=c(1992, 1995, 1997, 2000, 2002, 2005), minor_breaks=seq(1992, 2005, 1)) +
  geom_segment(aes(xend=max(x)+.20, yend=y), linetype="dotted", size=.1, color=cols) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) 