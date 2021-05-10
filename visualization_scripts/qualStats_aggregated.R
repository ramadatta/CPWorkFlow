
library(tidyr)
library(dplyr)

library(ggplot2)
library(scales)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/1Reads")

stats <- read.csv("Paeruginosa_183samples_AITBiotech_Pre_AND_Post_Trimming_Stats.csv", header = TRUE, sep = ",")
head(stats)

############--------CHANGE THESE VARIABLES---IMPORTANT!!------##########

samples_or_project_name <- "Paeruginosa"
seqCompany <- "AITBiotech" #Internal/AITBiotech
seqData_trimmed_or_not <- "Pre and Post-Trimmed" # use "Raw data" before trimming else use "After Trimming with Q30 score"  
genomeSize <- 6264404 


# Plotting

###----> Genome Coverage

ggplot(stats[order(pmax(stats$Genome_Coverage_Pre_Trimming, stats$Genome_Coverage_Post_Trimming)),], aes (x=factor(Isolate, levels=Isolate), Genome_Coverage)) +  
  geom_point(aes(y=Genome_Coverage_Pre_Trimming, color = "Genome_Coverage_Pre_Trimming"), size=1) + 
  geom_point(aes(y=Genome_Coverage_Post_Trimming, color = "Genome_Coverage_Post_Trimming"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  xlab("Isolates") +
  ylab("Genome Coverage ") +
  ggtitle(paste0("Genome Coverage of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("GenomeCoverage.png",dpi = 600, width = 15, height = 6, units = "in")

###----> Total Reads

ggplot(stats[order(pmax(stats$TotalReads_Pre_Trimming, stats$TotalReads_Post_Trimming)),], aes (x=factor(Isolate, levels=Isolate), Total_Reads)) +  
  geom_point(aes (y= TotalReads_Pre_Trimming, color = "TotalReads_Pre_Trimming"), size=1) + 
  geom_point(aes(y=TotalReads_Post_Trimming, color = "TotalReads_Post_Trimming"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("Total Reads") +
  ggtitle(paste0("Total Reads of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))

ggsave("TotalReads.png",dpi = 600, width = 15, height = 6, units = "in")

###----> Mean Read length

ggplot(stats[order(pmax(stats$MeanReadLength_Pre_Trimming, stats$MeanReadLength_Post_Trimming)),], aes (x=factor(Isolate, levels=Isolate), MeanReadLength)) +  
  geom_point(aes (y= MeanReadLength_Pre_Trimming, color = "MeanReadLength_Pre_Trimming"), size=1) + 
  geom_point(aes(y=MeanReadLength_Post_Trimming, color = "MeanReadLength_Post_Trimming"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  xlab("Isolates") +
  ylab("Mean Read Length ") +
  ggtitle(paste0("Mean Read Length of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))

ggsave("MeanReadLength.png",dpi = 600, width = 15, height = 6, units = "in")

  
###----> Total Bases

ggplot(stats[order(pmax(stats$TotalBases_Pre_Trimming, stats$TotalBases_Post_Trimming)),], aes (x=factor(Isolate, levels=Isolate), TotalBases)) +  
  geom_point(aes (y= TotalBases_Pre_Trimming, color = "TotalBases_Pre_Trimming"), size=1) + 
  geom_point(aes(y=TotalBases_Post_Trimming, color = "TotalBases_Post_Trimming"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("Total Bases") +
  ggtitle(paste0("Total Bases of ",samples_or_project_name," samples from ",seqCompany," (",seqData_trimmed_or_not,")"))

ggsave("TotalBases.png",dpi = 600, width = 15, height = 6, units = "in")

