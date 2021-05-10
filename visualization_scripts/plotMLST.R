setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/3mlst")

mlst_log <- read.table("log_mlst",sep = "\t", header = FALSE)
head(mlst_log)

names(mlst_log) <- c("filename", "organism", "ST", "gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7")
head(mlst_log)

library(ggplot2)
# Don't map a variable to y
ggplot(mlst_log, aes(x=factor(ST)))+
  geom_bar(stat="count", width=0.1, fill="tomato3") +
  geom_text(stat='count',aes(label=..count..),vjust=-1,size=3) +
  theme_minimal() + ggtitle("Pseudomonas aeruginosa data MLST Results") + # for the main title
xlab("Pseudomonas aeruginosa ST") + # for the x axis label
ylab("Count")

ggsave("MLST.png",dpi = 600, width = 8, height = 6, units = "in")

# Don't map a variable to y
ggplot(mlst_log, aes(x=factor(organism)))+
  geom_bar(stat="count", width=0.1, fill="tomato3") +
  geom_text(stat='count',aes(label=..count..),vjust=-1,size=3) +
  theme_minimal() + ggtitle("Species Counts") + # for the main title
  xlab("Species") + # for the x axis label
  ylab("Count")

ggsave("organism.png",dpi = 600, width = 8, height = 6, units = "in")
