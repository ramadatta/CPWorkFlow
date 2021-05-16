setwd("~/R/VisABRIcate/")

library(reshape)
library(ggplot2)
library(dplyr)
library(esquisse)

abricate_full_tbl <- read.csv("PAE_183_abricate_ncbi_resistome.csv")
head(abricate_full_tbl)

abricate_subset_tbl <- abricate_full_tbl %>% filter(X.COVERAGE == "100") %>% select(X.FILE,GENE,RESISTANCE)
head(abricate_subset_tbl)

ggplot(abricate_subset_tbl,aes(GENE,X.FILE)) + 
  geom_tile(aes(fill = RESISTANCE)) +
  coord_equal(ratio = 1) +
  scale_fill_gradient(low="white", high="#009C95") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 9), 
        legend.position = "none") +
  ggtitle("Paeruginosa data abricate-Resistome Results") + # for the main title
  xlab("Samples") + # for the x axis label
  ylab("Resistant Genes") #+ coord_flip()  
#esquisse

ggplot(abricate_subset_tbl, aes(x = X.FILE, y = interaction(RESISTANCE,GENE), fill = "red")) +
         geom_tile(aes(width=0.5,height=0.5)) + scale_fill_hue(direction = 1) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 5)) +
  theme(legend.position = "none")

ggplot(abricate_subset_tbl) +
  aes(x = X.FILE, y = interaction(RESISTANCE,GENE), fill = RESISTANCE) +
  geom_tile(aes(width=1,height=2)) + coord_fixed(ratio = 1) +
  scale_fill_hue(direction = 1) +
  theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) +
   theme(legend.position = "none")

ggsave("abricate_Resistome_hq.png",dpi = 600, width = 20, height = 20, units = "in")
