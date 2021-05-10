library(Biostrings)
library(ggmsa)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/18_abacus_results/MUSCLE")
nucl_sequences <- "PAE_183samples_combined_abacus_samples_with_refICE_faster.afa"
aln = readDNAMultipleAlignment(nucl_sequences)
ggmsa(nucl_sequences, start = 265, end = 300) 


#Set the reference as the 1st sequence, some Rattus, you can also use the consensus with consensusString() :
  
aln = unmasked(aln)
names(aln)[2]
ref = aln[2]

#Here we iterate through the sequence and make the binary for where the sequences are the same as the reference:

bm = sapply(1:length(aln),function(i){
  as.numeric(as.matrix(aln[i])==as.matrix(ref))
})

bm = t(bm)
rownames(bm) = names(aln)

# The plot you see above has sequences reversed, so to visualize the same thing we reverse it, and also subset on 265 - 300:
  
library(pheatmap)
pdf(file = "example.pdf", height = 10, width = 10)
pheatmap(bm[nrow(bm):1,1:90819],cluster_rows=FALSE,cluster_cols=FALSE)
dev.off()
