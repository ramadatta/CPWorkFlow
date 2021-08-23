setwd("/home/prakki/Documents/LeaRn/proLolli/MultiVCF")
library(trackViewer)

snpCount_df <- read.table("snpCount_multiVCF.list",sep = "\t", header = FALSE)
colnames(snpCount_df) <- c("SNP","Count")
head(snpCount_df)
nrow(snpCount_df)
features <- GRanges("chr1", IRanges(c(1211336,1393791,1889465,1889903,1891393,3365108,4840935), width=c(2898,1187,401,1223,1649,2003,1340), names=c("23S_ribosomal_RNA","mdh_2","hisI_2","yciC_1","feoB_2","tktA_2","lpfC"), fill=c("black", "gold1", "brown", "#CAB2D6", "steelblue4", "dodgerblue2", "yellow4"), height = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08)))

SNP <- c(1214150,1214158,1394646,1889686,1889758,1889767,1889908,1891472,3366305,4841537,5431837)
length(SNP)

sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)),
                     color = sample.int(6, length(SNP), replace=TRUE),
                     score = snpCount_df$Count
                     )
lolliplot(sample.gr, features)


