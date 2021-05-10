library(waffle)
library(extrafont)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/13_Prediction_Chr_Plasmids/rfplasmids")
parts <- c("Plasmid" = 101, "Chromosome" = 80, "No CPgene" = 2)

waffle_graph <- waffle(parts,
       rows=5, legend_pos = "bottom", xlab = ""
       # colors=c("tomato3", "yellow", "steelblue")
       #colors=c("#CC0000", "#006600", "#669999", "#00CCCC", "#FF9999", "#FF9900", "blue")
       )
ggsave("wafflePlot.png",dpi = 150, width = 5, height = 2, units = "in")
