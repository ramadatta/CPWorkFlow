
library(ggtree)
library(plotly)
library(ggplotlyExtra)
library(ggtreeExtra)
library(tidyverse)
library(glue)
library(tidyr)


Plasmid_SNPtree_SNPlocations_Heatmap <- function(location, treeInput, modifiedVCF,genePresAbsCSV, genePresAbsTAB, plot_title,
                                                 GenePA_offset,SpeciesST_offset,Plasmidlen_offset,SNPs_offset,Hospital_offset,
                                                 GenePA_pwidth,SpeciesST_pwidth,Plasmidlen_pwidth,SNPs_pwidth,Hospital_pwidth,var_width) {
  
  # Load the modified vcf file
  vcf <- data.table::fread(file = modifiedVCF, fill = TRUE, header = TRUE, sep = "\t", check.names = FALSE) # need input
  
  vcf_reqcols <- vcf %>% select(-"CHROM",-"ID",-"ALT",-"QUAL",-"FILTER",-"INFO",-"FORMAT")
  
  # Extracting data for SNP plot
  snp_data_ggtree <- reshape2::melt(vcf_reqcols, id=c("POS","REF"),variable.name = "Plasmid", value.name = "BASE") %>%
    filter(REF != BASE) %>% 
    #filter(BASE != "-") %>% 
    #filter(BASE != "N") %>% 
    select(Plasmid,POS) %>% 
    arrange(Plasmid,POS)
  
  # Add metadata
  #meta <- read.table("/data02/Analysis/Projects/2_CPE_Transmission/Plasmid_roary/Plasmids_1071_PCA_plot/Closed_CP_Plasmids_1071_complete_v4_Plasmid_meta_for_gubbins.csv", sep = ",", header = TRUE)
  
  meta <- data.table::fread("/data02/Analysis/Projects/2_CPE_Transmission/Reference_Excels/Closed_CP_Plasmids_metadata_20230213.csv", header= TRUE, sep = ",")
  
  # Mutate a new data column with both species and ST
  meta <- meta %>% mutate(Species_ST=paste0(Lab_Species,"_ST",ST)) 
  
  # To load and add annotation to the genes since the roary only puts group_names
  heatmapData <- read.table(genePresAbsTAB, header = TRUE, check.names = FALSE, sep = "\t")
  GeneAnnot_mat <- read.csv(genePresAbsCSV, header = TRUE, check.names = FALSE, sep = ",")
  
  heatmapData2 <- GeneAnnot_mat %>% 
    select(Gene, Annotation) %>% 
    right_join(heatmapData, by = 'Gene') %>% unite(Gene_Annotation, c("Gene", "Annotation"))
  
  write.csv(heatmapData2,file = glue('{location}geneAnnotation_presence_absence_intermediate1.csv'))
  
  #View(head(heatmapData2))
  
  tmp <- as.data.frame(t(heatmapData2[,-1]))
  colnames(tmp) <- heatmapData2$Gene
  
  write.csv(tmp,file = glue('{location}geneAnnotation_presence_absence_intermediate2.csv'))
  #View(head(tmp))
  
  tmp2 <- tmp %>% 
    tibble::rownames_to_column(var="Plasmid") %>% 
    left_join(meta %>% select(Plasmid, Plasmid_Clusters_renamed), by = 'Plasmid') %>% 
    relocate(Plasmid_Clusters_renamed, .after = Plasmid) 
  
  write.csv(tmp2,file = glue('{location}geneAnnotation_presence_absence_intermediate3.csv'))
  
  tmp3 <- tmp %>% 
    tibble::rownames_to_column(var="Plasmid") %>% 
    left_join(meta %>% select(Plasmid, Plasmid_Clusters_renamed), by = 'Plasmid') %>% 
    relocate(Plasmid_Clusters_renamed, .after = Plasmid) %>% group_by(Plasmid_Clusters_renamed) %>% summarize_if(is.numeric, sum, na.rm=TRUE)
  
  write.csv(tmp3,file = glue('{location}geneAnnotation_presence_absence_intermediate4.csv'))
  
  # For geom_fruit it is easier If I have long format, so generating tmp_long dataframe
  tmp_long <- tmp %>% 
    tibble::rownames_to_column("fasta") %>% 
    reshape2::melt(variable = "Gene", value.name = "Pres_Abs") 
  
  # Load the postgubbins tree
  #thetree <- read.tree(glue('{location}postGubbins.final_tree.tre')) # need input
  thetree <- read.tree(glue('{location}{treeInput}')) # need input
  
  q <- ggtree(thetree, colour = "orange", ladderize = TRUE, right = TRUE, size = 0.25) +
 #  q <-  ggtree(thetree,branch.length='none') +
    #geom_tiplab(size=2) +
    geom_tippoint(size = 1,shape=0, color="red") + 
    geom_rootpoint() + 
    theme_tree("white") +
    geom_treescale() +
    labs(title=glue('{plot_title}')) # need input
    

#ggplotly(temp)

#=========
  # Plotting and adding layers using geom_fruit

  # Plot1 Gene Presence Absence Plot

  qq1 <- q + geom_fruit(data=tmp_long,geom=geom_tile,
                        mapping = aes(y=fasta, x=Gene, fill=factor(Pres_Abs)),
                        # color = "white",
                        offset=GenePA_offset,pwidth=GenePA_pwidth,
                        axis.params = list(axis="x", title = "Gene Presence Absence Matrix", title.size = 3, text.size=2,
                                           vjust=0.5, line.size=0, hjust = 1, text.angle=90),
                        # grid.params=list(color="white")
  ) +
    scale_fill_manual(values=c("1"="#0000FF","0"="#FFA500"),
                      guide=guide_legend(ncol = 1, order=1,title.position = "top")) +
    labs(fill="Gene Presence/Absence") +
    theme(legend.position= 'bottom', plot.title = element_text(size=15, hjust = 0.5))

  #View(qq1)
  
  qq2 <- qq1

  ## Plot Cluster layer
  # qq2 <- qq1 + ggnewscale::new_scale_fill() +
  #   geom_fruit(data=meta,geom=geom_tile,
  #              mapping = aes(y=Plasmid,fill=factor(Plasmid_Clusters_renamed)),
  #              offset=Hospital_offset,pwidth=Hospital_pwidth, width=var_width,
  #              axis.params = list(axis="x", title = "Cluster", title.size = 3, text.size=2,
  #                                 vjust=0.5, line.size=0, hjust = 1, text.angle=90)) +
  #   # scale_fill_manual(values = c('red', 'blue', 'green','orange', "#00AFBB","#000000"),
  #   #                   guide = guide_legend(ncol=1, order = 5,title.position = "top")) +
  #   labs(fill="Plasmid_Cluster",)
  
  ## # Plot2 Species and ST Layer
  # qq2 <- qq1 + ggnewscale::new_scale_fill() +
  #    geom_fruit(data=meta,geom=geom_tile,
  #               mapping = aes(y=Plasmid,x=Species_ST,fill=Lab_Species),
  #               offset=SpeciesST_offset,pwidth=SpeciesST_pwidth,
  #               axis.params = list(axis="x", title = "Species and ST", title.size = 3, text.size=2,
  #                                  vjust=0.5, line.size=0, hjust = 1, text.angle=90),grid.params=list()) +
  #                scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  #                                           "#00AFBB", "#E7B800", "#FC4E07","#4dbfbc"),
  #                                  guide = guide_legend(ncol=2,order = 2,title.position = "top")) +
  #                labs(fill="Species")

  #qq2

  # Plot3 Plasmid Lengths layer
  qq3 <- qq2 + ggnewscale::new_scale_fill() +
    geom_fruit(data=meta,geom=geom_point,
               mapping = aes(y=Plasmid,x=Length,color="Plasmid Length"),shape=0,
               offset=Plasmidlen_offset,pwidth=Plasmidlen_pwidth,
               axis.params = list(axis="x", title = "Plasmid Lengths (bases)", title.size = 3, text.size=2,
                                  vjust=0.5, line.size=0, hjust = 1, text.angle=90),grid.params=list()) +
    scale_color_manual(values=c("blue","#0a7a68"),guide = guide_legend(order = 3,title.position = "top"))

  #qq3

  # Plot3 SNPs in Plasmid Core Genome layer
  qq4 <- qq3 + ggnewscale::new_scale_fill() +
    geom_fruit(data=snp_data_ggtree,geom=geom_point,
               mapping = aes(y=Plasmid,x=POS, color="SNP"),shape=1,
               offset=SNPs_offset,pwidth=SNPs_pwidth,#,GenePA_offset=0.1,SpeciesST_offset=0.1,Plasmidlen_offset=0.2,SNPs_offset=0.25,Hospital_offset=0.25,
               axis.params = list(axis="x", title = "SNPs in core genome (bases)", title.size = 3, text.size=2,
                                  vjust=0.5, line.size=0, hjust = 1, text.angle=90),grid.params=list()) +
    scale_color_manual(values=c("blue","#00AFBB"),
                       guide=guide_legend(ncol = 1, order=4,title.position = "top")) +
    labs(color="Scatter Points")

  #qq4
  # Plot5 Hospital Information layer
  qq5 <- qq4 + ggnewscale::new_scale_fill() +
    geom_fruit(data=meta,geom=geom_tile,
               mapping = aes(y=Plasmid,fill=factor(Hospital)),
               offset=Hospital_offset,pwidth=Hospital_pwidth, width=var_width,
               axis.params = list(axis="x", title = "Hospital", title.size = 3, text.size=2,
                                  vjust=0.5, line.size=0, hjust = 1, text.angle=90)) +
    scale_fill_manual(values = c('red', 'blue', 'green','orange', "#00AFBB","#000000"),
                      guide = guide_legend(ncol=1, order = 5,title.position = "top")) +
    labs(fill="Hospital")

  qq6 <- qq5

  # qq6 <- qq5 + ggnewscale::new_scale_fill() +
  #               geom_fruit(data=meta,geom=geom_tile,
  #               mapping = aes(y=Plasmid,fill=factor(Plasmid_Clusters_renamed)),
  #               offset=Hospital_offset,pwidth=Hospital_pwidth, width=var_width,
  #               axis.params = list(axis="x", title = "Cluster", title.size = 3, text.size=2,
  #                                  vjust=0.5, line.size=0, hjust = 1, text.angle=90)) +
  #               # scale_fill_manual(values = c('red', 'blue', 'green','orange', "#00AFBB","#000000"),
  #               #                   guide = guide_legend(ncol=1, order = 5,title.position = "top")) +
  #              labs(fill="Plasmid_Clusters_renamed")

  qq7 <- qq6 + vexpand(.4, -1) + vexpand(.1,1)
  #
  qq7
  
}

#Cluster 2 blaKPC-2 plasmids
#folderloc <- "/data02/Analysis/Projects/2_CPE_Transmission/AfterFeb102023/Multilayer_plots/PC1/gff/CoreGenome_AND_gubbins/"

folderloc <- "/data02/Analysis/Projects/2_CPE_Transmission/AfterFeb102023/Multilayer_plots/PC1/gff/CoreGenome_AND_gubbins_cd100/"

Plasmid_SNPtree_SNPlocations_Heatmap(
  location=folderloc,
  #treeInput=glue('postGubbins.final_tree.tre'),
  treeInput=glue('postGubbins.filtered_polymorphic_sites_fasttree.tre'),
  modifiedVCF=glue('{folderloc}modified.vcf'),
  genePresAbsCSV=glue('{folderloc}gene_presence_absence.csv'),
  genePresAbsTAB=glue('{folderloc}gene_presence_absence.Rtab'),
  plot_title="Plasmid Cluster 1 (PC1) carrying blaKPC-2 gene",
  GenePA_offset=0.1,SpeciesST_offset=0.1,Plasmidlen_offset=0.2,SNPs_offset=0.25,Hospital_offset=0.25,
  GenePA_pwidth=3,SpeciesST_pwidth=2,Plasmidlen_pwidth=1,SNPs_pwidth=1,Hospital_pwidth=0.07, var_width=0.02
  )

ggsave(filename = "/data02/Analysis/Projects/2_CPE_Transmission/AfterFeb102023/Multilayer_plots/PC1/gff/CoreGenome_AND_gubbins_cd100/PC1.jpeg", width = 17, height = 10, device='jpeg', dpi=700)



pdf(file = "/data02/Analysis/Projects/Misc/test_gubbins/Plot2_rootpoint_nodelabels.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 50) # The height of the plot in inches

Plasmid_SNPtree_SNPlocations_Heatmap(
  location=folderloc,
  treeInput=glue('postGubbins.final_tree.tre'),
  modifiedVCF=glue('{folderloc}modified.vcf'),
  genePresAbsCSV=glue('{folderloc}gene_presence_absence.csv'),
  genePresAbsTAB=glue('{folderloc}gene_presence_absence.Rtab'),
  plot_title="Plasmid Cluster 1 (PC1) carrying blaKPC-2",
  GenePA_offset=0.1,SpeciesST_offset=0.1,Plasmidlen_offset=0.2,SNPs_offset=1,Hospital_offset=0.25,
  GenePA_pwidth=3,SpeciesST_pwidth=2,Plasmidlen_pwidth=1,SNPs_pwidth=1,Hospital_pwidth=0.07, var_width=5
)


# Step 3: Run dev.off() to create the file!
dev.off()


