
library(tidyr)
library(dplyr)

library(ggplot2)
library(scales)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/1Reads/postTrimmedDataStats_RPlots")

stats <- read.table("PAE_183samples_PostTrim_Q30_Q20_readstats.txt",sep = "\t")
head(stats)

############--------CHANGE THESE VARIABLES---IMPORTANT!!------##########

samples_or_project_name <- "Paeruginosa"
seqCompany <- "AITBiotech" #Internal/AITBiotech
seqData_trimmed_or_not <- "Post-Trimmed" # use "Raw data" before trimming else use "After Trimming with Q30 score"  
genomeSize <- 6264404 #Assuming the bacterial genome size is 5 Mb


############--------ADD COLNAMES TO TABLE ------##########
names(stats) <- c("filename","TotalReads","TotalBases","Q20Bases","Q30Bases","MeanReadLength")
head(stats)

stats_wHeader <- stats %>% separate(filename, c('Isolate', 'ReadFile'), sep = '/', convert = TRUE)
head(stats_wHeader)

aggregated_stats <- as.data.frame(stats_wHeader 
                                  %>% group_by(Isolate) 
                                  %>% summarise(TotalReads = sum(TotalReads),
                                                TotalBases = sum(TotalBases),
                                                Q20Bases=sum(Q20Bases),
                                                Q30Bases=sum(Q30Bases),
                                                MeanReadLength=mean(MeanReadLength)
                                                )
                                  )

head(aggregated_stats)

# Percentage of Q20 Bases in both R1 and R2 fastq
aggregated_stats$Q20_Prcnt <- aggregated_stats$Q20Bases/aggregated_stats$TotalBases

# Percentage of Q30 Bases in both R1 and R2 fastq
aggregated_stats$Q30_Prcnt <- aggregated_stats$Q30Bases/aggregated_stats$TotalBases

#Genome Coverage
aggregated_stats$Genome_Coverage <- aggregated_stats$TotalBases/genomeSize # Pseudomonas aeruginosa PAO1, complete genome (https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2)

head(aggregated_stats)


# Plotting

###----> Genome Coverage

ggplot(aggregated_stats, aes(x = reorder(Isolate, Genome_Coverage), y = Genome_Coverage)) +
  geom_point(colour="purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  xlab("Isolates") +
  ylab("Genome Coverage ") +
  ggtitle(paste0("Genome Coverage of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("GenomeCoverage.png",dpi = 600, width = 15, height = 6, units = "in")

###----> Total Reads

ggplot(aggregated_stats, aes(x = reorder(Isolate, TotalReads), y = TotalReads)) +
  geom_point(colour="purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("Total Reads") +
  ggtitle(paste0("Total Reads of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))

ggsave("TotalReads.png",dpi = 600, width = 15, height = 6, units = "in")

###----> Mean Read length

ggplot(aggregated_stats, aes(x = reorder(Isolate, MeanReadLength), y = MeanReadLength)) +
  geom_point(colour="purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  xlab("Isolates") +
  ylab("Mean Read Length ") +
  ggtitle(paste0("Mean Read Length of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))


ggsave("MeanReadLength.png",dpi = 600, width = 15, height = 6, units = "in")

  
###----> Total Bases
ggplot(aggregated_stats, aes(x = reorder(Isolate, TotalBases), y = TotalBases)) +
  geom_point(colour="purple") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("Total Bases") +
  ggtitle(paste0("Total Bases of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))


ggsave("TotalBases.png",dpi = 600, width = 15, height = 6, units = "in")

write.csv(aggregated_stats,paste0(samples_or_project_name,"_samples_",seqCompany,"_",seqData_trimmed_or_not,"_aggregated_stats.txt"))
