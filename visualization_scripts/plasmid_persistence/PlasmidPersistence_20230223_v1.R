library(lubridate)
library(ggplot2)
library(forcats)
library(glue)
library(dplyr)
library(ggrepel)

CP_Plasmids_meta_ClusterNumber_cutoff0.9_df <- data.table::fread("/data02/Analysis/Projects/2_CPE_Transmission/Reference_Excels/Closed_CP_Plasmids_metadata_07102022.csv", header= TRUE, sep = ",")

head(CP_Plasmids_meta_ClusterNumber_cutoff0.9_df)

Clusters_info_df <- CP_Plasmids_meta_ClusterNumber_cutoff0.9_df %>% 
  select(Plasmid_Clustering_KmerMethod_0.9cutoff, Inc_group_regrouped, Relaxase_type_corrected, Length,ENT_ID,PID) %>% 
  group_by(Plasmid_Clustering_KmerMethod_0.9cutoff) %>% 
  summarise(Inc_group_regrouped_collapse=paste(unique(Inc_group_regrouped),collapse='#'),
            Relaxase_type_corrected_collapse=paste(unique(Relaxase_type_corrected),collapse='#'),
            IsolateCount=n(),
            PatientCount = n_distinct(PID),
            PlasmidLength_Min = min(Length),
            PlasmidLength_Max = max(Length),
  ) %>% 
  rename(ClusterNum=Plasmid_Clustering_KmerMethod_0.9cutoff,
         Inc_Group=Inc_group_regrouped_collapse,
         Relaxase=Relaxase_type_corrected_collapse) 

Closed_CP_Plasmids_metadata_df <- data.table::fread("/data02/Analysis/Projects/2_CPE_Transmission/Reference_Excels/Closed_CP_Plasmids_metadata_07102022.csv", header= TRUE, sep = ",")

#View(Closed_CP_Plasmids_metadata_df)

#Closed_CP_Plasmids_metadata_df_subset <- Closed_CP_Plasmids_metadata_df %>% select(CapesID,Clonal_cluster)
Closed_CP_Plasmids_metadata_df_subset <- Closed_CP_Plasmids_metadata_df %>% select(ENT_ID,Clonal_cluster)


#head(Closed_CP_Plasmids_metadata_df_subset)
#View(Closed_CP_Plasmids_metadata_df)

epi <- data.table::fread("/data02/Analysis/Projects/2_CPE_Transmission/Plasmid_roary/Plasmids_1071_PCA_plot/from_Natascha/hospitalwardbedoverlaps/EpiOverlap_Bacterial_Combined_v5_updatedControls_forSharing.csv")


#head(epi)

#epi <- left_join(epi, Closed_CP_Plasmids_metadata_df_subset, by=c("Index_CaPES_ID"="CapesID")) %>% head()
epi_1 <- left_join(epi, Closed_CP_Plasmids_metadata_df_subset, by=c("Index_ENT_ID"="ENT_ID")) 


#head(epi_1)

#View(epi_1)

cc <- epi_1 %>% 
  filter(grepl('Cluster', Clonal_cluster)) %>% 
  select(Clonal_cluster) %>% arrange(Clonal_cluster) %>%  unique()  %>% c()

#remove df1
rm(cc_epi_df)
cc_epi_df <- data.frame()
#cc$Clonal_cluster

# Sum loop
for(ClusterNum in cc$Clonal_cluster) {
  #print(paste0(ClusterNum,"HI"))
  
  cc_epi_df = rbind(cc_epi_df, 
                    epi_1 %>% 
                      filter(Clonal_cluster == ClusterNum) %>% 
                      select(Clonal_cluster,Status) %>% 
                      mutate(Status = strsplit(as.character(Status), "#")) %>% 
                      tidyr::unnest(Status) %>% 
                      filter(trimws(Status) !="") %>% 
                      mutate(Short_Status=case_when(stringr::str_detect(Status,"Hospital Indirect") ~ "HI", 
                                                    stringr::str_detect(Status,"Hospital Direct" ) ~ "HD", 
                                                    stringr::str_detect(Status,"Ward Indirect" ) ~ "WI",  
                                                    stringr::str_detect(Status,"Ward Direct" ) ~ "WD",
                                                    TRUE ~ "")) %>% 
                      #head()
                      group_by(Clonal_cluster) %>% 
                      arrange(Short_Status) %>% 
                      summarize(Collapsed_Status = paste(unique(Short_Status), collapse = "|"))
  )
} 

Closed_CP_Plasmids_metadata_df2 <- left_join(Closed_CP_Plasmids_metadata_df, cc_epi_df, by=c("Clonal_cluster")) %>%
  mutate(Collapsed_Status = stringr::str_replace(Collapsed_Status, "^\\|", ""))

write.csv(Closed_CP_Plasmids_metadata_df2,"/data02/Analysis/Projects/2_CPE_Transmission/Plasmid_roary/Plasmids_1071_PCA_plot/Closed_CP_Plasmids_1071_complete_v6_Plasmid_meta_for_gubbins.csv")  

#View(Closed_CP_Plasmids_metadata_df2)

Phylodynamics_staticplot_withepi <- function(plot_title,clusterNum,pointsize){
  #View(meta)
  # Subsetting to a smaller set for ease of simplicity - selecting OXA-48 Plasmids records with n=41 rows
  subset_clusterNum <- Closed_CP_Plasmids_metadata_df2 %>% 
    #filter(Plasmid_Clustering_KmerMethod_0.9cutoff==clusterNum) %>%  # Although the variable number
    filter(Plasmid_Clustering_KmerMethod_0.9cutoff %in% clusterNum) %>% 
    mutate(Unique_Plasmid_Cluster = if_else(is.na(Clonal_cluster),paste0("Singleton_", row_number()), as.character(Clonal_cluster))) %>%  # assign unique id for NA values
    group_by(Unique_Plasmid_Cluster) %>% 
    mutate(group_id = cur_group_id()) #%>% View()
  
  subset_clusterNum$Date_of_culture <- as.Date(subset_clusterNum$Date_of_culture, format = "%d/%m/%Y")
  
  #View(subset_clusterNum)
  
  labelInfo <- subset_clusterNum %>%
    ungroup() %>% 
    mutate(Unique_Plasmid_Cluster = forcats::fct_reorder(Unique_Plasmid_Cluster, Date_of_culture, min)) %>% 
    select(Unique_Plasmid_Cluster,Date_of_culture,Collapsed_Status) %>% 
    group_by(Unique_Plasmid_Cluster) %>% 
    filter(Date_of_culture == max(Date_of_culture)) %>% 
    unique() 
  
  labelInfo$Collapsed_Status <- factor(labelInfo$Collapsed_Status)
  
  cc_hospCount <- subset_clusterNum %>%
    filter(grepl('Cluster', Clonal_cluster)) %>%
    #unique() %>%
    group_by(Clonal_cluster) %>%
    mutate(HospitalCount_byKmerCluster = n_distinct(Hospital)) %>%
    select(Clonal_cluster,HospitalCount_byKmerCluster) %>% unique()
  
  labelInfo2 <- left_join(labelInfo,cc_hospCount,by=c("Unique_Plasmid_Cluster"= "Clonal_cluster")) %>% 
    mutate(Collapsed_Status_withHospCnt_byKC = paste(Collapsed_Status, HospitalCount_byKmerCluster, sep = ";")) %>% 
    mutate(Collapsed_Status_withHospCnt_byKC = stringr::str_replace(Collapsed_Status_withHospCnt_byKC, "NA;NA",""))
  
  labelInfo2$Collapsed_Status_withHospCnt_byKC <- factor(labelInfo2$Collapsed_Status_withHospCnt_byKC)
  
  pd1 <- subset_clusterNum %>%
    ungroup() %>% 
    mutate(Unique_Plasmid_Cluster = forcats::fct_reorder(Unique_Plasmid_Cluster, Date_of_culture, min)) %>%
    ggplot(aes(Date_of_culture, Unique_Plasmid_Cluster)) +
    geom_line(aes(group=Unique_Plasmid_Cluster), color='grey') +
    geom_point(aes(color = Lab_Species,shape= factor(Hospital)),size=pointsize) +
    geom_text_repel(data = labelInfo2, aes(x = Date_of_culture, y = Unique_Plasmid_Cluster, 
                                           label = Collapsed_Status_withHospCnt_byKC, 
                                           #color = Collapsed_Status
    ), 
    size=3, min.segment.length = 0, max.overlaps = Inf, nudge_x = 20
    #size=3, min.segment.length = 0, max.overlaps = Inf, force = 5,  nudge_y = 20, direction = "x"
    ) +
    # geom_label_repel(data = labelInfo, aes(x = Date_of_culture, y = Unique_Plasmid_Cluster, 
    #                                        label = Collapsed_Status, 
    #                                        #color = Collapsed_Status
    # ), size=3, min.segment.length = 0, max.overlaps = Inf, nudge_x = -40
    # ) +
    labs(title = plot_title,
         #subtitle="of a facet barplot in ggplot",
         tag = "Epidemiological information\n\nHD - Hospital Direct\nHI - Hospital Indirect\nWD - Ward Direct\nWI - Ward Indirect") +
    ylab("Cluster") + xlab("Date") +
    theme_bw() +
    theme(
      rect = element_blank(),
      legend.position = "right",
      legend.key=element_rect(fill='gray96'),
      legend.title =element_text(size=10),
      text=element_text(size=12),
      axis.title.x = element_text(vjust = 0, size = 11),
      axis.title.y = element_text(vjust = 2, size = 11),
      #axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
      axis.text.y = element_blank(),
      plot.background = element_rect(fill = '#fbf9f4', color = '#fbf9f4'),
      panel.grid = element_line(color = "#b4aea9"),
      plot.tag.position = c(.865,.17),
      plot.tag = element_text(hjust =0, size=9), 
      panel.grid.major = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())
  
  pd1
}

Phylodynamics_staticplot_withepi(plot_title="Plasmid persistence with blaNDM-1 gene in Cluster 6",clusterNum=6,pointsize=1)