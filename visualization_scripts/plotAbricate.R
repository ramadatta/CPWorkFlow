setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/4Abricate/")

library(reshape)
library(ggplot2)

abricate_log <- read.table("PAE_183_abricate_vfdb.summary_dotsReplacedwithZero.csv",sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(abricate_log)

abricate_Modlog <- abricate_log %>% select(-one_of("NUM_FOUND"))

abricate_Modlog[-1] <- as.integer(abricate_Modlog[-1] != 0)
head(abricate_Modlog)
abricate_Modlog.melt <- melt(abricate_Modlog)
head(abricate_Modlog.melt)
tail(abricate_Modlog.melt)
dim(abricate_Modlog.melt)

ggplot(abricate_Modlog.melt,aes(variable,FILE)) + 
  geom_tile(aes(fill = value)) +
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

ggsave("abricate_Resistome_hq.png",dpi = 300, width = 20, height = 12, units = "in")




# ggplot(melt(res),aes(sample_name, variable, fill=value)) +
#   # tile with black contour
#   geom_tile() + 
#   # B&W theme, no grey background
#   theme_bw() + 
#   # square tiles
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
#         axis.text.y = element_text(angle = 0, hjust = 1, size = 5), 
#         legend.position = "none") +
#   coord_equal() + 
#   # Green color theme for `fill`
#   scale_fill_distiller(palette="Greens", direction=1) + 
#   # printing values in black
#  # geom_text(aes(label=value), color="black") +
#   # removing legend for `fill` since we're already printing values
#   guides(fill=F) +
#   # since there is no legend, adding a title
#   labs(title = "Count of fruits per person") #+ coord_flip()   
