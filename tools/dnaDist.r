library(ape)
library(dplyr)
library(reshape2)

## Read in fasta file
    data <- read.FASTA(file = "postGubbins.filtered_polymorphic_sites.fasta") # Read input multiple fasta file

## Calculate the pair-wise distance
# Route 1 #
out <-  dist.dna(data,model="N",pairwise.deletion=TRUE,as.matrix=T) ## Full matrix
D_out_melt = melt(as.matrix(out), varnames = c("row", "col"))
D_out_melt_sorted = arrange(D_out_melt, value)
write.csv(D_out_melt_sorted, "postGubbins.filtered_polymorphic_sites_full_triangle.csv")

out[lower.tri(out,diag=T)] <- NA ## take upper triangular matrix, when needed

# Route 2 #
# out <-  dist.dna(data,model="raw",pairwise.deletion=TRUE)

## Plot
options(scipen=10000)
# hist(as.matrix(out),breaks=1000,na.rm=T)
# plot(density(as.matrix(out),na.rm=T))
# qplot(value, data  = filtered, geom = "density",col = Grouping2) + labs(x= "SNP") + geom_vline(xintercept = 2, col = "blue") + facet_grid(Grouping2~.)
# qplot(value, data  = SNP, geom = "histogram", binwidth = 1000, fill = Grouping1, facets = . ~ Grouping1) + geom_vline(xintercept = 2) + labs(x= "SNP")
# ggplot(SNP, aes(value, colour = Grouping2)) + stat_ecdf() + geom_vline(xintercept = 2) + labs(x = "SNP", "") + geom_histogram(aes(y = 3*..density..), alpha = 0.2, binwidth = 3)



D_out_melt = melt(as.matrix(out), varnames = c("row", "col"))
D_out_melt_sorted = arrange(D_out_melt, value)
write.csv(D_out_melt_sorted, "postGubbins.filtered_polymorphic_sites_melt_sorted.csv")
