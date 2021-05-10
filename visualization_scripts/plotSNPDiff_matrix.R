setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/6SNPDifftable_SNIPPY")

library(splitstackshape)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gplots)

half_snp_table <- read.table("postGubbins.filtered_polymorphic_sites_melt_sorted.csv",sep = ",", header = TRUE)
head(half_snp_table)
colnames(half_snp_table) <- c("No","Sample1","Sample2","SNPDiff")
head(half_snp_table)
tail(half_snp_table)

half_snp_table <- na.omit(half_snp_table)
tail(half_snp_table)
#half_snp_table[1] <- NULL  #remove 1st column
class(half_snp_table)
head(half_snp_table)
tail(half_snp_table)

####Cumulative Frequency####################
# table(half_snp_table$SNPDiff) -> one way

# Let's try another way to generate frequency
half_snp_table %>% 
  group_by(SNPDiff) %>% 
  summarise(Freq = n()) %>%
#  ungroup %>% 
  mutate(CumFreq = cumsum(Freq))

Freq.df <-  half_snp_table %>%
  group_by(SNPDiff) %>%
  summarise(Freq = n())

tail(Freq.df)

# Plot full frequency

ggplot(data = Freq.df, mapping = aes(x = SNPDiff, y = Freq)) +
  geom_line(color="tomato3", size=0.7) + scale_x_continuous(breaks = scales::pretty_breaks(n = 30), limits = c(0, NA)) +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(0, NA)) +
 ggtitle("P.aeruginosa ST308 210 samples SNP Difference Frequency") + # for the main title +
 xlab("SNP Difference") + # for the x axis label
 ylab("Frequency")

ggsave("ST308 SNP Frequency.png",dpi = 600, width = 15, height = 6, units = "in")

# Plot partial frequency

ggplot(data = Freq.df, mapping = aes(x = SNPDiff, y = Freq)) +
  geom_line(color="tomato3", size=0.7) + scale_x_continuous(breaks = scales::pretty_breaks(n = 30), limits = c(0, 100)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(0, NA)) +
  ggtitle("P.aeruginosa ST308 210 samples SNP Difference Frequency") + # for the main title +
  xlab("SNP Difference") + # for the x axis label
  ylab("Frequency") +  geom_vline(aes(xintercept = 50), 
                                  linetype = "dashed", color="black",size = 0.5)
ggsave("ST308 SNP Frequency X-axis cutoff.png",dpi = 600, width = 15, height = 6, units = "in")

####Cumulative Sum##########################

# half_snp_table$CSNPDiff <- cumsum(half_snp_table$SNPDiff)
# head(half_snp_table)
# tail(half_snp_table)
# 
# library(reshape)
# df_snpdifflteq10 <- half_snp_table %>% filter(SNPDiff <=50)
# head(df_snpdifflteq10)
# tail(df_snpdifflteq10)
# ggplot(df_snpdifflteq10) + geom_point(aes(x = No, y = CSNPDiff), color = 'blue', size = 0.1) +
# scale_x_continuous(breaks = scales::pretty_breaks(n = 30), limits = c(0, NA)) +
# scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(0, NA)) +
# ggtitle("P.aeruginosa ST308 SNPDiff Cumulative Sum (SNP Cutoff <=1)") + # for the main title +
# xlab("Pair") + # for the x axis label
# ylab("Cumulative SNP Diff Sum")
# 
# 
# ggplot(half_snp_table) + geom_point(aes(x = No, y = CSNPDiff), color = 'blue', size = 0.1) + scale_x_continuous(breaks = scales::pretty_breaks(n = 30), limits = c(0, NA)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 20), limits = c(0, NA))
#############################

snp_table <- read.table("postGubbins.filtered_polymorphic_sites_full_triangle.csv",sep = ",", header = TRUE)
head(snp_table)
colnames(snp_table) <- c("No","Sample1","Sample2","SNPDiff")
head(snp_table)
tail(snp_table)

complete_snp_table <- na.omit(snp_table)
tail(complete_snp_table)
complete_snp_table[1] <- NULL 
class(complete_snp_table)
complete_snp_table

# Convert long snp table to wide snp table
reshape(complete_snp_table, idvar = "Sample1", timevar = "Sample2", direction = "wide")
library(reshape)
wide_complete_snp_table <- cast(complete_snp_table, Sample1 ~ Sample2)

## To generate a heatmap with SNP Difference values
library(pheatmap)

pheatmap(wide_complete_snp_table, display_numbers = F,fontsize = 6, main = "P.aeruginosa ST308 SNPDiff matrix")

SNPDiff_hm <- pheatmap(wide_complete_snp_table, display_numbers = F,fontsize = 6, main = "P.aeruginosa ST308 SNPDiff matrix")

save_pheatmap_png <- function(x, filename, width=2500, height=1500, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(SNPDiff_hm, "ST308_SNPDiff_Matrix_pheatmap.png")




####Clustering with distance matrix ###########

# distMat <- acast(complete_snp_table, Sample1~Sample2, value.var="SNPDiff")
# eucl_dist <- dist(distMat,method="euclidean")
# 
# hc <- hclust(eucl_dist)
# hc
# 
# heatmap.2(as.matrix(eucl_dist,Rowv=as.dendrogram(hc),Colv=as.dendrogram(hc)))
# 
# #plot(hc)
# 
# pdf("Heatmap_1071_small_similarity.pdf", height=10, width=20)
# heatmap.2(as.matrix(eucl_dist), Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc), col=rev(rich.colors(100)), tracecol=NA,  key = TRUE, keysize = 1, key.xlab = "Distance", cexRow=0.3, cexCol=0.3, margins=c(10,10))
# dev.off()
# 
# wide_complete_snp_table <- acast(complete_snp_table, Sample1~Sample2, value.var="SNPDiff")
# #head(melt(res))
# head(complete_snp_table)
# 
# ggplot(complete_snp_table, aes(Sample1, Sample2, fill = SNPDiff, alpha = SNPDiff)) + 
#   geom_tile() +
#   #geom_text(aes(label=SNPDiff)) +
#   scale_alpha_identity(guide = "none") +
#   coord_equal(expand = 0) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ggtitle("P.aeruginosa ST308 SNPDiff matrix") + # for the main title
#   xlab("Samples") + # for the x axis label
#   ylab("Samples")+ scale_fill_gradient(low = "white",  high = "#009C95")
# 
# ggsave("ST308_SNPDiff_Matrix.png",dpi = 600, width = 8, height = 6, units = "in")

