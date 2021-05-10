## Calculate Pairwise SNP differences and Run David Eyre's Poisson script

library(ape)
library(dplyr)
library(reshape2)
library(igraph)

  setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/6SNPDifftable_SNIPPY/")

STGroup <- "PAE_ST308_210_Samples" ###Change the group
beast_mu <- 87.902116928

## Read in fasta file
  data <- read.FASTA(file = "postGubbins.filtered_polymorphic_sites.fasta") # Read input multiple fasta file

## Calculate the pair-wise distance
# Route 1 #
out <-  dist.dna(data,model="N",pairwise.deletion=TRUE,as.matrix=T) ## Full matrix
out[lower.tri(out,diag=T)] <- NA ## take upper triangular matrix, when needed

# Route 2 #
# out <-  dist.dna(data,model="raw",pairwise.deletion=TRUE)

## Plot
options(scipen=10000)

D_out_melt = melt(as.matrix(out), varnames = c("row", "col"))
D_out_melt_sorted = arrange(D_out_melt, value)

#Ignore the NA records
fullrecords <-  D_out_melt_sorted[complete.cases(D_out_melt_sorted),] 

#Replacing unwanted strings in the filebasenames
fullrecords_NameStripped <- as.data.frame(sapply(fullrecords,function(x) {x <- gsub("_Ns_converted_to_RefBase","",x)}))

#Ignore rows with Reference Sequence. The reference header contains pattern called "length"
fullrecords_NameStripped <- fullrecords_NameStripped[grep("length", fullrecords_NameStripped$row, invert = TRUE), ]
fullrecords_NameStripped
fullrecords_NameStripped <- fullrecords_NameStripped[grep("length", fullrecords_NameStripped$col, invert = TRUE), ]
fullrecords_NameStripped

head(beast_mu)

#Adding pair.id column
fullrecords_NameStripped$PairID <- paste0(fullrecords_NameStripped$row,"#",fullrecords_NameStripped$col ) #paste and paste0 has different functions
head(fullrecords_NameStripped)

# Add a column to dataframe with BEAST SNP rate with header name "beast_mu"
fullrecords_NameStripped$beast_mu <- beast_mu
fullrecords_NameStripped$STGroup <- STGroup

head(fullrecords_NameStripped)
tail(fullrecords_NameStripped)



#write.csv(D_out_melt_sorted, "postGubbins_CCGroup3_ST14_460_Samples_subset_57_Samples_forGubbins.filtered_polymorphic_sites_melt.csv")

# David's Poisson Script Portion

#setwd("/data01/Server_copies/1CPE_Transmission/7SNP_Pipeline_2019/kpneumo/CCgroup3_460_Samples/7_Mapping_ALL_Samples_Gubbins_subset_57_samples_gubbins_bModelTest_RelaxedLorgNormal_clock/Poisson")
#df = read.csv("Kpneumo_Gubbins_57_Samples_Poisson_Input.csv", sep =',', header=T)

#Importing the date of cultures
cultureDate <- read.table("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/6SNPDifftable_SNIPPY/IsolateDOC2_including_Jeanettes_NUHsamples.txt",header = TRUE);
head(cultureDate)

# Using match function, retrieved the date of cultures for sample 1 and sample 2
fullrecords_NameStripped$Sample1_Date <-cultureDate$cultureDate[match(fullrecords_NameStripped$row,cultureDate$Isolate)]
fullrecords_NameStripped$Sample2_Date <-cultureDate$cultureDate[match(fullrecords_NameStripped$col,cultureDate$Isolate)]

head(fullrecords_NameStripped)
tail(fullrecords_NameStripped)

##Calculated Absolute time difference between samples
fullrecords_NameStripped$timeDiff <- abs(as.numeric(as.character(fullrecords_NameStripped$Sample1_Date)) - as.numeric(as.character(fullrecords_NameStripped$Sample2_Date)))
head(fullrecords_NameStripped)
tail(fullrecords_NameStripped)

##Calculate Expected SNP
fullrecords_NameStripped$expected_snp <- fullrecords_NameStripped$timeDiff*fullrecords_NameStripped$beast_mu
head(fullrecords_NameStripped)

##Calculate if there is a Tranmission with original cutoff
fullrecords_NameStripped$TranmissionByBEASTcutoff <- as.numeric(as.character(fullrecords_NameStripped$value)) < as.numeric(as.character(fullrecords_NameStripped$expected_snp))
head(fullrecords_NameStripped)

##Rarranging columns 

fullrecords_NameStripped <- fullrecords_NameStripped[,c("STGroup","row","col","PairID","value","beast_mu","Sample1_Date","Sample2_Date","timeDiff","expected_snp","TranmissionByBEASTcutoff")]
head(fullrecords_NameStripped)
#keep only first 9 columns

#df = df[,1:9]
#head(df)
#rename columns
colnames(fullrecords_NameStripped)
colnames(fullrecords_NameStripped) = c("st_cp_grouping", "id1", "id2", "pair.id", "snps","mu","Sample1_Date","Sample2_Date","time.diff", "expected.snp", "original.cutoff")

#check format
#head(df)

#add probability of number of snps or fewer
cutoff = 0.95

#df = cbind(df, ppois(df$snps,df$time.diff*df$mu, lower.tail=T), ppois(df$snps-1,df$time.diff*df$mu, lower.tail=T), ppois(df$snps-1,df$time.diff*df$mu, lower.tail=T)<=cutoff)
#colnames(df) = c("st_cp_grouping", "id1", "id2", "time.diff", "pair.id", "snps", "mu", "expected.snp", "original.cutoff", "prob.snps.or.less", "prob.snps.minus.1.or.less", "proposed.cutoff")

fullrecords_NameStripped$max.snp = qpois(0.95, fullrecords_NameStripped$time.diff*fullrecords_NameStripped$mu)
fullrecords_NameStripped$max.snp 

fullrecords_NameStripped$transmission.plausible = ifelse(as.numeric(as.character(fullrecords_NameStripped$snps)) <= as.numeric(as.character(fullrecords_NameStripped$max.snp)), 1, 0)

#check format
head(fullrecords_NameStripped)

#write output

#write.csv(fullrecords_NameStripped, "postGubbins_CCgroup3_ST11_33_Samples_PairwiseSNP.csv")
write.csv(fullrecords_NameStripped, "postGubbins_ST308_PAE_210samples_Poisson_Output_linkedPairs.csv")

# Get connected components based on "transmission.plausible==1"
TranPlausPairs <- fullrecords_NameStripped %>% filter(transmission.plausible==1) %>% select(id1,id2)

  g1 <- graph.data.frame(TranPlausPairs, directed = FALSE)
  g1
  plot(g1)
  
  #Fancy figure start
  n = 100
  p = 1.5/n
  g = erdos.renyi.game(n, p)
  coords = layout.fruchterman.reingold(g)
  plot(g1, layout=coords, vertex.size = 3, vertex.label=NA)
  plot(g1, vertex.label=NA, vertex.size=2, vertex.color="#0CCF02")
  #Fancy figure end
  
  cl1 <- clusters(g1)
  
  tbl1 <- cbind( V(g1)$name, cl1$membership )
  class(tbl1)
  tbl1 <- as.data.frame(tbl1)
  class(tbl1)
  head(tbl1)
  colnames(tbl1) = c("sample", "Cluster")
  
  write.table(tbl1,file="Samples_ClusterNumber.txt",row.names=FALSE) # drops the rownames and write to text file
  #return(tbl1)


