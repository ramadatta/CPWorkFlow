
library(tidyr)
library(dplyr)

library(ggplot2)
library(scales)

setwd("/data02/Analysis/Projects/6_Paeruginosa_183_samples_Shawn/2Assembly")

stats <- read.csv("full_gteq1kb_combined.csv", header = TRUE, sep = ",")
head(stats)

############--------CHANGE THESE VARIABLES---IMPORTANT!!------##########

samples_or_project_name <- "Paeruginosa"
seqCompany <- "AITBiotech" #Internal/AITBiotech
seqData_trimmed_or_not <- "Full and gteq1kb assemblies" # use "Raw data" before trimming else use "After Trimming with Q30 score"  
genomeSize <- 6264404 

# Plotting

###----> Number_of_Contigs

ggplot(stats[order(pmax(stats$n_contigs_full_assembly, stats$n_contigs_gteq1kb)),], aes (x=factor(Isolate, levels=Isolate), Contigs_per_Assembly)) +  
  geom_point(aes (y= n_contigs_full_assembly, color = "n_contigs_full_assembly"), size=1) + 
  geom_point(aes(y=n_contigs_gteq1kb, color = "n_contigs_gteq1kb"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  scale_y_continuous(breaks = scales::pretty_breaks(n=20))  +
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  xlab("Isolates") +
  ylab("Contigs_per_Assembly ") + 
  ggtitle(paste0("Contigs per Assembly from ",samples_or_project_name," data"," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("Contigs_per_Assembly.png",dpi = 600, width = 15, height = 6, units = "in")

###----> Assembly Size

ggplot(stats[order(pmax(stats$contig_bp_full_assembly, stats$contig_bp_gteq1kb)),], aes (x=factor(Isolate, levels=Isolate), AssemblySize)) +  
  geom_point(aes (y= contig_bp_full_assembly, color = "contig_bp_full_assembly"), size=1) + 
  geom_point(aes(y=contig_bp_gteq1kb, color = "contig_bp_gteq1kb"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("AssemblySize") +
  ggtitle(paste0("Assembly Size ",samples_or_project_name," data"," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("AssemblySize.png",dpi = 600, width = 15, height = 6, units = "in")

###----> N50 Size

ggplot(stats[order(pmax(stats$ctg_N50_full_assembly, stats$ctg_N50_gteq1kb)),], aes (x=factor(Isolate, levels=Isolate), N50Size)) +  
  geom_point(aes (y= ctg_N50_full_assembly, color = "ctg_N50_full_assembly"), size=1) + 
  geom_point(aes(y=ctg_N50_gteq1kb, color = "ctg_N50_gteq1kb"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("N50Size") +
  ggtitle(paste0("N50 Size ",samples_or_project_name," data"," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("N50Size.png",dpi = 600, width = 15, height = 6, units = "in")

#ctg_max

###----> Max Contig Size

ggplot(stats[order(pmax(stats$ctg_max_full_assembly, stats$ctg_max_gteq1kb)),], aes (x=factor(Isolate, levels=Isolate), LargestContigSize)) +  
  geom_point(aes (y= ctg_max_full_assembly, color = "ctg_max_full_assembly"), size=1) + 
  geom_point(aes(y=ctg_max_gteq1kb, color = "ctg_max_gteq1kb"), size=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),legend.position = "right") + #Rotating and spacing axis labels
  expand_limits(x = c(0, NA), y = c(0, NA)) +
  scale_y_continuous(labels = unit_format(unit = "MB", scale = 1e-6)) + ## need to load "scales" library - very readable
  xlab("Isolates") +
  ylab("Largest Contig Size") +
  ggtitle(paste0("Largest Contig Size ",samples_or_project_name," data"," (",seqData_trimmed_or_not,")")) # ggtitle("Genome Coverage of Steno Samples from AITBiotech (After Trimming)")

ggsave("LargestContigSize.png",dpi = 600, width = 15, height = 6, units = "in")




