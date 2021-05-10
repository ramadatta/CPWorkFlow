setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/4CGE")

library(splitstackshape)
library(dplyr)
library(stringr)
library(qdapTools)
library(ggplot2)
library(reshape2)

cge_log <- read.table("user-result-summary.txt",sep = "\t", header = TRUE)
head(cge_log)

cge_rest_split <- cSplit(cge_log, "resistance", ",")
head(cge_rest_split)

#df1 <- cge_rest_split %>% mutate(ST = as.numeric(str_extract(mlst, "[0-9]+"))) %>% arrange(ST) %>% select("sample_name","ST",contains("resistance_"))
df1 <- cge_rest_split %>% select("sample_name",contains("resistance_"))
head(df1)
df1 <- as.data.frame(df1)
class(df1)
df1
#res <- cbind(df1[1:2], mtabulate(as.data.frame(t(df1[-1:-2])))) # Print ST too
res <- cbind(df1[1], mtabulate(as.data.frame(t(df1[-1]))))
res
row.names(res) <- NULL
res

#If you prefer a different order, you can order them by hand:
res$sample_name <- factor(res$sample_name, levels=c("C02481","C02434","C02433","C02432","C02431","C02430","C02429","C02397","C02385","C02382","C02374","C02364","C02361","C02360","C02358","C02357","C02349","C02345","C02340","C02329","C02322","C02318","C02288","C02263","C02252","C02247","C02224","C02223","C02191","C02170","C02164","C02123","C02120","C02099","C02092","C02091","C02071","C02021","C01996","C01992","C01978","C01965","C01958","C01922","C01921","C01914","C01907","C01890","C01889","C01888","C01883","C01877","C01873","C01872","C01871","C01863","C01845","C01844","C01843","C01835","C01813","C01812","C01809","C01776","C01758","C01756","C01747","C01727","C01716","C01695","C01683","C01681","C01657","C01623","C01603","C01591","C01578","C01572","C01571","C01528","C01505","C01495","C01494","C01448","C01447","C01446","C01429","C01412","C01399","C01395","C01391","C01390","C01340","C01339","C01338","C01337","C01336","C01305","C01295","C01267","C01248","C01245","C01244","C01243","C01238","C01237","C01206","C01201","C01196","C01187","C01157","C01148","C01121","C01107","C01094","C01076","C01074","C00954","C00950","C00916","C00881","C00851","C00819","C00818","C00814","C00811","C00787","C00765","C00746","C00718","C00716","C00707","C00703","C00696","C00669","C00651","C00625","C00584","C00569","C00507","C00499","C00496","C00488","C00487","C00476","C00472","C00461","C00432","C00401","C00400","C00391","C00389","C00387","C00386","C00369","C00337","C00336","C00318","C00309","C00281","C00280","C00270","C00220","C00219","C00194","C00186","C00183","C00177","C00136","C00132","C00131","C00130","C00118","C00073","C00068","C00056","C00042","C00027","C00024","C01020","C00255","C02244","C00189"))


# ggplot(melt(res), aes(variable,sample_name, fill = value, alpha = value)) + 
#   geom_tile(colour = "gray50") +
#   scale_alpha_identity(guide = "none") +
#   coord_equal(expand = 0) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 5), legend.position = "none") +
#   ggtitle("Paeruginosa data CGE-Resistome Results") + # for the main title
#   xlab("Resistant Genes") + # for the x axis label
#   ylab("Samples") + coord_flip()
# 
# ggsave("CGE_Resistome.png",dpi = 600, width = 15, height = 6, units = "in")

ggplot(melt(res),aes(sample_name,variable)) + 
  geom_tile(aes(fill = value)) +
  coord_equal(ratio = 1) +
  scale_fill_gradient(low="white", high="#009C95") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 9), 
        legend.position = "none") +
  ggtitle("Paeruginosa data CGE-Resistome Results") + # for the main title
  xlab("Samples") + # for the x axis label
  ylab("Resistant Genes") #+ coord_flip()  

ggsave("CGE_Resistome_hq.png",dpi = 300, width = 20, height = 4, units = "in")




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
