library(lubridate)
library(ggplot2)
library(forcats)
library(glue)
library(dplyr)
library(ggrepel)
library(ggtext)

library(showtext)
font_add_google(name = "Roboto Condensed", family = "roboto-condensed")

# Changes from previous v1: removing epi information from the clonal cluster points, adding species, ST, keeping hospital count and isolate counts

meta_df <- data.table::fread("/data02/Analysis/Projects/2_CPE_Transmission/Reference_Excels/Closed_CP_Plasmids_metadata_20230213.csv", header= TRUE, sep = ",")

meta_df <- meta_df %>% 
  mutate(ShortLab_Species=case_when(stringr::str_detect(Lab_Species,"Escherichia coli") ~ "EC", 
                                              stringr::str_detect(Lab_Species,"Klebsiella pneumoniae" ) ~ "KP", 
                                              stringr::str_detect(Lab_Species,"Enterobacter cloacae complex" ) ~ "ECloc_comp",
                                              stringr::str_detect(Lab_Species,"Enterobacter cloacae" ) ~ "ECloc",
                                              stringr::str_detect(Lab_Species,"Citrobacter freundii" ) ~ "CF",
                                              stringr::str_detect(Lab_Species,"Citrobacter species" ) ~ "Cit_sp",
                                              stringr::str_detect(Lab_Species,"Enterobacter aerogenes" ) ~ "EA",
                                              stringr::str_detect(Lab_Species,"Klebsiella oxytoca" ) ~ "KO",
                                              stringr::str_detect(Lab_Species,"Morganella morganii" ) ~ "MM",
                                              stringr::str_detect(Lab_Species,"Escherichia species" ) ~ "Esc_sp",
                                              TRUE ~ "Others"))

head(meta_df)

Plasmid_Persistence <- function(plot_title,
                                clusterNum,
                                CC_countlabel_x,CC_annotlabel_x,HGT_annotlabel_x,
                                CC_countlabel_y,CC_annotlabel_y,HGT_annotlabel_y)
  {
  
subset_clusterNum <- meta_df %>% 
  #filter(Plasmid_Clusters_renamed==clusterNum) %>%  # Although the variable number
  filter(Plasmid_Clusters_renamed %in% clusterNum) %>% 
  mutate(Unique_Plasmid_Cluster = if_else(is.na(Clonal_cluster),paste0("Singleton_", row_number()), as.character(Clonal_cluster))) %>%  # assign unique id for NA values
  group_by(Unique_Plasmid_Cluster) %>% 
  mutate(group_id = cur_group_id())

#View(subset_clusterNum)

subset_clusterNum$Date_of_culture <- as.Date(subset_clusterNum$Date_of_culture, format = "%d/%m/%Y")

labelInfo <- subset_clusterNum %>%
  filter(n()>1) %>% #This is not to double count if the clonal cluster has unclustered plasmid
  ungroup() %>% 
  mutate(Unique_Plasmid_Cluster = forcats::fct_reorder(Unique_Plasmid_Cluster, Date_of_culture, min)) %>% 
  select(Unique_Plasmid_Cluster,Date_of_culture,ShortLab_Species,ST) %>%
  group_by(Unique_Plasmid_Cluster) %>% 
  filter(Date_of_culture == max(Date_of_culture)) %>% 
  unique() 

cc_isolate_hospCount <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  group_by(Clonal_cluster) %>%
  mutate(IsolateCount_byCC = n_distinct(ENT_ID)) %>%
  mutate(HospitalCount_byCC = n_distinct(Hospital)) %>%
  select(Clonal_cluster,IsolateCount_byCC, HospitalCount_byCC) %>% 
  unique()

labelInfo2 <- left_join(labelInfo,cc_isolate_hospCount,by=c("Unique_Plasmid_Cluster"= "Clonal_cluster")) %>%
  mutate(CC_label = paste(ShortLab_Species,ST,IsolateCount_byCC,HospitalCount_byCC, sep = ";")) %>% 
  mutate(CC_label = stringr::str_replace(CC_label, ".*NA;NA",""))

labelInfo2$CC_label <- factor(labelInfo2$CC_label)

#View(labelInfo2)

# Annotate plot with Clonal Cluster count, Total Isolates in CCs and Plasmid Backbone and unique patients 

## Overall isolate and patient count in PCX cluster

total_isolatecount <- subset_clusterNum %>% ungroup() %>% select(Plasmid_name_filename) %>% n_distinct()
total_patientcount <- subset_clusterNum %>% ungroup() %>% select(PID) %>% unique() %>% n_distinct()

## Number of patients in Horizontal gene transfer (backbone)

uniq_subj_HGT <- subset_clusterNum %>%
  select(Date_of_culture,Clonal_cluster,ENT_ID,PID) %>% 
  arrange(Unique_Plasmid_Cluster, Date_of_culture) %>% slice(1) %>% 
  ungroup() %>% 
  select(PID) %>% 
  unique() %>% n_distinct()

## Number of isolates in Horizontal gene transfer (backbone)

total_isolates_HGT <- subset_clusterNum %>%
  select(Date_of_culture,Clonal_cluster,ENT_ID,PID) %>% 
  arrange(Unique_Plasmid_Cluster, Date_of_culture) %>% slice(1) %>% 
  ungroup() %>% 
  select(ENT_ID) %>% 
  unique() %>% n_distinct()

## Number of patients in Clonal Clusters

uniq_subj_CC <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  select(Date_of_culture,Clonal_cluster,ENT_ID,PID) %>% 
  filter(n()>1) %>% #This is not to double count if the clonal cluster has unclustered plasmid
  arrange(Unique_Plasmid_Cluster, Date_of_culture) %>% slice(2:n()) %>% 
  ungroup() %>% 
  select(PID) %>% 
  unique() %>% 
  n_distinct()

## Number of isolates in in Clonal Clusters

total_isolates_CC <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  select(Date_of_culture,Clonal_cluster,ENT_ID,PID) %>% 
  filter(n()>1) %>% #This is not to double count if the clonal cluster has unclustered plasmid
  arrange(Unique_Plasmid_Cluster, Date_of_culture) %>% slice(2:n()) %>% 
  ungroup() %>% 
  select(ENT_ID) %>% 
  unique() %>% 
  n_distinct()

# Number of isolates in clonal cluster

#Note: some CC's are 0 because the other plasmid is not in the current plasmid cluster and is unclustered plasmid.

cc_isolatecount_df <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  #group_by(Clonal_cluster) %>%
  mutate(IsolateCount_byCC = n_distinct(ENT_ID)) %>%
  select(Clonal_cluster,IsolateCount_byCC) %>% 
  unique() %>% 
  mutate(IsolateCount_exclFirstIsolate=(IsolateCount_byCC - 1)) # not counting the first isolate

cc_isolatecount <- sum(cc_isolatecount_df$IsolateCount_exclFirstIsolate)
cc_isolatecount_total <- sum(cc_isolatecount_df$IsolateCount_byCC)

# Number of clonalclusters

cc_count <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  #group_by(Clonal_cluster) %>%
  mutate(IsolateCount_byCC = n_distinct(ENT_ID)) %>%
  select(Clonal_cluster,IsolateCount_byCC) %>% 
  unique() %>% 
  mutate(IsolateCount_exclFirstIsolate=(IsolateCount_byCC - 1)) %>% 
  filter(IsolateCount_exclFirstIsolate>0) %>% 
  ungroup() %>% 
  select(Clonal_cluster) %>% 
  unique() %>% n_distinct()
  
#cc_count

# Number of unique patients in clonalcluster

cc_count_totalPIDs <- subset_clusterNum %>%
  filter(grepl('Cluster', Clonal_cluster)) %>%
  #group_by(Clonal_cluster) %>%
  mutate(IsolateCount_byCC = n_distinct(ENT_ID)) %>%
  select(Clonal_cluster,IsolateCount_byCC,PID) %>% 
  unique() %>% 
  mutate(IsolateCount_exclFirstIsolate=(IsolateCount_byCC - 1)) %>% 
  #filter(IsolateCount_exclFirstIsolate>0) %>% 
  ungroup() %>% 
  select(PID) %>% 
  unique() %>% n_distinct()

# Let make a table of the above generated results

Num <- c("1", "2", "3", "4")
Description <- c("Total Clonal Cluster Count:","Total Isolate Count:",  "Total Isolates in Clonal Clusters (excluded the firstmost isolate):", "Total Isolates involved in HGT:")

Count <- c(cc_count,total_isolatecount,  total_isolates_CC, total_isolates_HGT)
Unique_Patient_Count <- c(cc_count_totalPIDs, total_patientcount, uniq_subj_CC,uniq_subj_HGT)

df <- data.frame(Num, Description, Count,Unique_Patient_Count)
#df

# Let's create annotations using the values from above dataframe

CC_countlabel <- glue::glue("<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> Clonal Clusters: ",
                            "<span style='font-size:15pt;font-family:Roboto Condensed;color:#ff8c00'>",cc_count,"</span><br/>")

CC_annotlabel <- glue::glue("<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> Clonal Transmission: ",
                            "<span style='font-size:15pt;font-family:Roboto Condensed;color:#A020F0'>",total_isolates_CC,"</span>",
                            "<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> isolates (",
                            "<span style='font-size:10pt;font-family:Roboto Condensed;color:#FF0000'>",uniq_subj_CC,"</span>", 
                            "<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> unique subjects)</span><br/>")

HGT_annotlabel <- glue::glue("<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> Plasmid Transmission: ",
  "<span style='font-size:15pt;font-family:Roboto Condensed;color:#008080'>",total_isolates_HGT,"</span>", 
  "<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> isolates (",
  "<span style='font-size:10pt;font-family:Roboto Condensed;color:#FF0000'>",uniq_subj_HGT,"</span>", 
  "<span style='font-size:10pt;font-family:Roboto Condensed;color:#000000'> unique subjects)</span><br/>")


## Plasmid persistence plot

#piratepal(palette = "basel")


mybasel_palette_edited <- c("#0C5BB0FF",	#blue1
  "#EE0011FF", #red
  "#15983DFF", #green
  "#EC579AFF", #pink
  "#FA6B09FF", #orange
  "#149BEDFF", #blue2
  "#A1C720FF", #green2
  "#FEC10BFF", #yellow
  "#16A08CFF", #turquoise
  "#9A703EFF", #poop
  "maroon",
  "purple"
)
# temp_df <- subset_clusterNum %>%
#   ungroup() %>%
#   mutate(Unique_Plasmid_Cluster = forcats::fct_reorder(Unique_Plasmid_Cluster, Date_of_culture, min))

gg1 <- subset_clusterNum %>%
  ungroup() %>%
  mutate(Unique_Plasmid_Cluster = forcats::fct_reorder(Unique_Plasmid_Cluster, Date_of_culture, min)) %>%
  ggplot(aes(Date_of_culture, Unique_Plasmid_Cluster)) +
  geom_line(aes(group=Unique_Plasmid_Cluster), color="#E30B5C") + #'grey'
  geom_point(aes(color = Lab_Species,shape= factor(Hospital)),size=2, alpha = 0.5) +
  scale_shape_manual(values = c("1"= 0, "2"= 1,"3"= 2,"4"= 3, "5"=5, "6"= 15)) +
  #scale_color_paletteer_d("yarrr::basel", 12, type = "continuous") +
  #paletteer_dynamic("cartography::multi.pal", 20)
  scale_color_manual(values = mybasel_palette_edited) +
  scale_y_discrete(expand = expansion(add = 1.2)) +
  geom_text_repel(data = labelInfo2, aes(x = Date_of_culture, y = Unique_Plasmid_Cluster,
                                         label = CC_label, family = 'Roboto Condensed'
                                         #color = Collapsed_Status
  ),
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
  size=3.5, min.segment.length = 0, max.overlaps = Inf
  ) +
  # geom_text_repel(data = labelInfo2, aes(x = Date_of_culture, y = Unique_Plasmid_Cluster,
  #                                        label = CC_label, family = 'Roboto Condensed'
  #                                        #color = Collapsed_Status
  # ),
  # size=3.5, min.segment.length = 0, max.overlaps = Inf, nudge_x = 20
  # ) +
  labs(
      title = plot_title,
      #title = "Persistence of **Plasmid Cluster2** (PC2) carrying blaNDM-1 gene",
       #subtitle="of a facet barplot in ggplot",
      # tag = "Label Information\n\nSpecies; ST; Isolate Count;\nHospital Count",
       shape="Hospital", colour="Species",
       caption = "**Label Information:** Species; ST; Isolate Count; Hospital Count"
  ) +
  ylab("Plasmid/Clonal Clusters") + xlab("Date") +
  theme_bw(base_family="Roboto Condensed") +
  theme(
    #rect = element_blank(),
    plot.title = element_textbox_simple(colour = "black", size = 14, margin = margin(b = 20), family = "Roboto Condensed"),
    plot.title.position = "plot",
    legend.position = "right",
    legend.key=element_rect(fill='gray96'),
    legend.title =element_text(size=10),
    text=element_text(size=12),
    axis.title.x = element_text(vjust = 0, size = 11),
    axis.title.y = element_text(vjust = 2, size = 11),
    #axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_blank(),
    #plot.background = element_rect(fill = '#F0F5F5', color = '#F0F5F5'),
    #panel.background = element_rect(fill = "#F0F5F5", colour = "#F0F5F5"),
    #legend.background = element_rect(fill="#F0F5F5", colour ="#F0F5F5"),
    panel.grid = element_line(color = "#b4aea9"),
    plot.tag.position = c(.890,.25),
    plot.tag = element_text(hjust =0, size=9),
    panel.grid.major = element_blank(),
   # panel.grid.minor = element_line(color = 2, size = 0.25, linetype = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.caption = element_markdown(lineheight = 1.2),
    plot.margin = margin(1,1,1.5,1.2, "cm")
   
    ) +
  annotate(geom = "richtext", x = CC_countlabel_x, y = CC_countlabel_y,  label = CC_countlabel,fill = NA, family= "Roboto Condensed", label.color = NA) +
  annotate(geom = "richtext", x = CC_annotlabel_x, y = CC_annotlabel_y,  label = CC_annotlabel,fill = NA, family= "Roboto Condensed", label.color = NA) +
  annotate(geom = "richtext", x = HGT_annotlabel_x, y = HGT_annotlabel_y,  label = HGT_annotlabel,fill = NA, label.color = NA) +
  coord_cartesian(clip = "off")

gg1
}

p1 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster1** (PC1) carrying blaKPC-2 gene",
                    CC_countlabel_x=as.Date("2012-11-15"),CC_countlabel_y=330,
                    CC_annotlabel_x=as.Date("2013-02-01"),CC_annotlabel_y=316,
                    HGT_annotlabel_x=as.Date("2013-02-01"),HGT_annotlabel_y=300,
                    clusterNum="PC1")
p1

p2 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster2** (PC2) carrying blaNDM-1 gene",
                    CC_countlabel_x=as.Date("2011-10-15"),CC_countlabel_y=210,
                    CC_annotlabel_x=as.Date("2012-02-01"),CC_annotlabel_y=200,
                    HGT_annotlabel_x=as.Date("2012-02-01"),HGT_annotlabel_y=190,
                    clusterNum="PC2")

# p3 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster3** (PC3) carrying blaOXA-181 gene",
#                     CC_countlabel_x=as.Date("2011-11-20"),CC_countlabel_y=16,
#                     CC_annotlabel_x=as.Date("2012-02-01"),CC_annotlabel_y=15,
#                     HGT_annotlabel_x=as.Date("2012-01-25"),HGT_annotlabel_y=14,
#                     clusterNum="PC3")

p3 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster3** (PC3) carrying blaOXA-181 gene",
                          CC_countlabel_x=as.Date("2012-09-27"),CC_countlabel_y=16,
                          CC_annotlabel_x=as.Date("2012-12-01"),CC_annotlabel_y=15,
                          HGT_annotlabel_x=as.Date("2012-11-30"),HGT_annotlabel_y=14,
                          clusterNum="PC3")
#p3
p4 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster4** (PC4) carrying blaOXA-48 gene",
                    CC_countlabel_x=as.Date("2012-11-27"),CC_countlabel_y=23,
                    CC_annotlabel_x=as.Date("2013-02-01"),CC_annotlabel_y=22,
                    HGT_annotlabel_x=as.Date("2013-01-31"),HGT_annotlabel_y=21,
                    clusterNum="PC4")

p5 <- Plasmid_Persistence(plot_title="Persistence of **Plasmid Cluster5** (PC5) carrying blaNDM-1 gene",
                    CC_countlabel_x=as.Date("2012-12-01"),CC_countlabel_y=23,
                    CC_annotlabel_x=as.Date("2013-02-01"),CC_annotlabel_y=22,
                    HGT_annotlabel_x=as.Date("2013-01-31"),HGT_annotlabel_y=21,
                    clusterNum="PC5")

# ggpubr::ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")
# ggpubr::ggarrange(p3, p4, ncol=2, common.legend = TRUE, legend="bottom")
# ggpubr::ggarrange(p3, p4, p5, ncol=3, common.legend = TRUE, legend="bottom")

p1
#ggsave(filename = "PC1.jpeg", width = 12, height = 9, device='jpeg', dpi=700)
ggsave(filename = "PC1_edit.jpeg", width = 11, height = 12, device='jpeg', dpi=300)
p2
#ggsave(filename = "PC2.jpeg", width = 12, height = 9, device='jpeg', dpi=700)
ggsave(filename = "PC2_edit.jpeg", width = 11, height = 12, device='jpeg', dpi=300)
p3
#ggsave(filename = "PC3.jpeg", width = 10, height = 7, device='jpeg', dpi=700)
ggsave(filename = "PC3_edit.jpeg", width = 9, height = 9, device='jpeg', dpi=300)
p4
#ggsave(filename = "PC4.jpeg", width = 10, height = 7, device='jpeg', dpi=700)
ggsave(filename = "PC4_edit.jpeg", width = 9, height = 9, device='jpeg', dpi=300)
p5
#ggsave(filename = "PC5.jpeg", width = 10, height = 7, device='jpeg', dpi=700)
ggsave(filename = "PC5_edit.jpeg", width = 10, height = 9, device='jpeg', dpi=300)
