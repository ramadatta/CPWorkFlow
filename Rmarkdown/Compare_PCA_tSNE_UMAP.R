# load the Rtsne package
library(Rtsne)
library(ggplot2)
library(plotly)
library(dplyr)
library(magrittr)
library(ggrepel)
library(crosstalk) # For extracting data points from plot
library(DT) # data table

# load the gene presence/absence matrix
Gene_PA_df <- read.table("/data02/Analysis/Projects/2_CPE_Transmission/Plasmid_roary/Plasmids_1071_PCA_plot/gff/gene_presence_absence.Rtab", quote = "", check.names = FALSE, header = TRUE, sep = "\t")

# Transpose to samples in rows and genes in columns
Gene_PA_df_t  <- data.table::transpose(Gene_PA_df, keep.names = "Plasmid", make.names = "Gene")   

#View(Gene_PA_df_t)

# Load metadata file
Gene_PA_meta <- read.table("/data02/Analysis/Projects/2_CPE_Transmission/Plasmid_roary/Plasmids_1071_PCA_plot/Closed_CP_Plasmids_1071_complete_v3_Plasmid_meta.csv", check.names = FALSE, header = TRUE, sep = ",", fill = TRUE)
head(Gene_PA_meta)
# Combine the Matrix with Metadata
Gene_PA_df_with_metadata <- left_join(Gene_PA_df_t, Gene_PA_meta, by = c("Plasmid"))
View(Gene_PA_df_with_metadata)

Gene_PA_df_with_metadata$year=format(as.Date(Gene_PA_df_with_metadata$Date_of_culture, format="%d/%m/%Y"),"%Y")

# Convert first column to rownames
roaryTab.df  <- Gene_PA_df_with_metadata %>% remove_rownames %>% column_to_rownames(var="Plasmid")

############################################# Run PCA ##############################################################
# res.pca.scaleF  <- prcomp(roaryTab.df[1:(length(roaryTab.df)-17)], scale = FALSE)
# 
# options(ggrepel.max.overlaps = 0)
# # With the fviz_pca_ind function, we cannot have two categorical variables: one point shape and one fill color.
# # So we modify the output from fviz_pca_ind(), so you would need to take out the data from the results, and 
# # plot it again using ggplot2:
# 
# #######------------Multiple labels---------###########
# #basic_plot <- fviz_pca_ind(res.pca.scaleF, label="none")
# basic_plot <- fviz_pca_ind(res.pca.scaleF, label="var") # Adding var gives nicer legend symbols than none
# 
# #View(head(basic_plot$data))
# #tail(basic_plot$data)
# 
# # Combine the Matrix with Metadata
# basic_plot_with_meta <- left_join(basic_plot$data, Gene_PA_meta, by = c("name" = "Plasmid"))
# #View((basic_plot_with_meta))
# 
# # With confidence level
# 
# p3_rev <- ggplot(basic_plot_with_meta,
#                  aes(x=x,y=y,col=CP_Gene,shape=CP_Gene)) + geom_point(size=2, alpha = 0.3) +
#   theme_bw() +
#   stat_ellipse(level=0.95,linetype=2) +
#   ggtitle(paste0("PCA plot of all CP gene plasmids based on Inc group")) +
#   xlab("PC1") +
#   ylab("PC2")
# 
# p3_rev

#ggplotly(p3_rev)

############################################# PCA DONE ##############################################################

############################################# TSNE ##############################################################
## Run the t-SNE algorithm and store the results into an object called tsne_results
View(Gene_PA_df_with_metadata[1:(length(roaryTab.df)-14)]) 
View(Gene_PA_df_with_metadata[(length(roaryTab.df)-14):length(roaryTab.df)]) 

tsne.norm <- Rtsne(Gene_PA_df_with_metadata[1:(length(roaryTab.df)-14)], perplexity=10, check_duplicates = FALSE)

#plasmid_clusters <- Gene_PA_df_with_metadata[ ,2039] 
#info.norm <- Gene_PA_df_with_metadata[,2030:2045] 

col_start <- length(roaryTab.df)-14
col_end <- length(roaryTab.df)

info.norm <- Gene_PA_df_with_metadata[,col_start:col_end] 

#View(info.norm)

## Generate the t_SNE plot
par(mfrow=c(1,2)) # To plot two images side-by-side
plot(tsne.norm$Y, col = "blue", pch = 19, cex = 1) # Plotting the first image
#plot(tsne.norm$Y, col = "black", bg= plasmid_clusters, pch = 21, cex = 1) # Second plot: Color the plot by the real species type (bg= IR_species)

head(tsne.norm$Y)

info.norm %<>% mutate(tsne1 = tsne.norm$Y[, 1], tsne2 = tsne.norm$Y[,2])                                                               
head(info.norm)
ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = CP_Gene)) + 
  geom_point(alpha = 0.3) + theme_bw()

info.norm %>% 
  filter(CP_Gene=="blaKPC2") %>% 
ggplot(aes(x = tsne1, y = tsne2, colour = factor(Hospital))) + 
  geom_point(alpha = 0.3) + theme_bw()

info.norm %>% 
  filter(CP_Gene=="blaNDM1") %>% 
  ggplot(aes(x = tsne1, y = tsne2, colour = factor(Hospital))) + 
  geom_point(alpha = 0.3) + theme_bw()

# tSNE different clustering methods

# # Hierarchical clustering
# options(ggrepel.max.overlaps = Inf)
# hc.norm = hclust(dist(tsne.norm$Y))
# info.norm$hclust = factor(cutree(hc.norm, 10))
# hc.norm.cent = info.norm %>% group_by(Plasmid_cluster) %>% select(tsne1, 
#                                                          tsne2) %>% summarize_all(mean)
# hc_cl10 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = factor(Plasmid_cluster))) + 
#   geom_point(alpha = 0.6) + theme_bw() + geom_label_repel(aes(label = Plasmid_cluster), 
#                                                           data = hc.norm.cent) + guides(colour = "none") + ggtitle("KPC plasmids by plasmid clusters") 
# 
# ggplotly(hc_cl10)
# 
# 
# # Hierarchical clustering
# options(ggrepel.max.overlaps = Inf)
# hc.norm = hclust(dist(tsne.norm$Y))
# info.norm$hclust = factor(cutree(hc.norm, 10))
# hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1, 
#                                                          tsne2) %>% summarize_all(mean)
# hc_cl10 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) + 
#   geom_point(alpha = 0.6) + theme_bw() + geom_label_repel(aes(label = hclust), 
#                                                           data = hc.norm.cent) + guides(colour = "none") + 
#   ggtitle("KPC plasmids by plasmid clusters") 
# 
# hc.norm = hclust(dist(tsne.norm$Y), method = "ward.D")
# info.norm$hclust = factor(cutree(hc.norm, 10))
# hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1, 
#                                                          tsne2) %>% summarize_all(mean)
# hc_ward_cl10 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) + 
#   geom_point(alpha = 0.6) + theme_bw() + geom_label_repel(aes(label = hclust), 
#                                                           data = hc.norm.cent) + guides(colour = "none") + 
#   ggtitle("KPC plasmids by plasmid clusters") 
# 
# hc.norm = hclust(dist(tsne.norm$Y), method = "ward.D")
# info.norm$hclust = factor(cutree(hc.norm, 18))
# hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1, 
#                                                          tsne2) %>% summarize_all(mean)
# hc_ward_cl18 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) + 
#   geom_point(alpha = 0.6) + theme_bw() + geom_label_repel(aes(label = hclust), 
#                                                           data = hc.norm.cent) + guides(colour = "none") + 
#   ggtitle("KPC plasmids by plasmid clusters") 
# 
# # Kmeans
# km.norm = kmeans(tsne.norm$Y, 18, nstart = 100)
# info.norm$kmeans = factor(km.norm$cluster)
# km.cent = info.norm %>% group_by(kmeans) %>% select(tsne1, 
#                                                     tsne2) %>% summarize_all(mean)
# kmeans_cl18 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = kmeans)) + 
#   geom_point(alpha = 0.3) + theme_bw() + geom_label_repel(aes(label = kmeans), 
#                                                           data = km.cent) + guides(colour = "none") + ggtitle("kmeans 9 clusters")
# 
# #Mclust
# library(mclust)
# 
# mc.norm = Mclust(tsne.norm$Y, 18)
# info.norm$mclust = factor(mc.norm$classification)
# mc.cent = info.norm %>% group_by(mclust) %>% select(tsne1, 
#                                                     tsne2) %>% summarize_all(mean)
# mclust_18 <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = mclust)) + 
#   geom_point(alpha = 0.3) + theme_bw() + geom_label_repel(aes(label = mclust), 
#                                                           data = mc.cent) + guides(colour = "none") +  ggtitle("mclust 18 clusters")
# 
# #Density based clustering
# 
# library(fpc)
# ds.norm = dbscan(tsne.norm$Y, 2)
# info.norm$density = factor(ds.norm$cluster)
# ds.cent = info.norm %>% group_by(density) %>% select(tsne1, 
#                                                      tsne2) %>% summarize_all(mean)
# db_cl <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = density)) + 
#   geom_point(alpha = 0.3) + theme_bw() + geom_label_repel(aes(label = density), 
#                                                           data = ds.cent) + guides(colour = "none") +  ggtitle("density-based 18 clusters")
# 
# 
# library(igraph)
# library(FNN)
# k = 100
# knn.norm = get.knn(as.matrix(tsne.norm$Y), k = k)
# knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index), 
#                                  k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + 
#                                                                                       as.vector(knn.norm$nn.dist)))
# nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
# nw.norm = simplify(nw.norm)
# lc.norm = cluster_louvain(nw.norm)
# info.norm$louvain = as.factor(membership(lc.norm))
# lc.cent = info.norm %>% group_by(louvain) %>% select(tsne1, 
#                                                      tsne2) %>% summarize_all(mean)
# ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = louvain)) + 
#   geom_point(alpha = 0.3) + theme_bw() + geom_label_repel(aes(label = louvain), 
#                                                           data = lc.cent) + guides(colour = "none")
# 
# library(patchwork)
# 
# (hc_ward_cl18 + kmeans_cl18) / (mclust_18 + db_cl)


##UMAP
library(umap)

temp_umap <- Gene_PA_df_with_metadata[1:(length(roaryTab.df)-14)]

#test <- Gene_PA_df_with_metadata[(length(roaryTab.df)-14):length(roaryTab.df)]
#View(test)

umap_fit <- temp_umap %>%
  tibble::column_to_rownames("Plasmid") %>% 
  select(where(is.numeric)) %>% 
  scale() %>% 
  umap()

info.norm_umap <- Gene_PA_df_with_metadata[,col_start:col_end] 
test <- info.norm_umap
info.norm_umap %<>% mutate(umap1 = umap_fit$layout[, 1], umap2 = umap_fit$layout[, 2]) 


umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  tibble::rownames_to_column("Plasmid") 

head(umap_df)
info.norm_umap %<>% mutate(plasmid = umap_df$Plasmid, umap1 = umap_df$UMAP1, umap2 = umap_df$UMAP2)  

info.norm_umap %>% head()

ggplotly(ggplot(info.norm_umap, aes(x = umap1, y = umap2, colour = CP_Gene, text = paste("Plasmid_cluster:", Plasmid_cluster, "<br>Lab_Species:",Lab_Species, "<br>ST:",ST, "<br>Length:",Length))) + 
  geom_point(alpha = 0.3) + theme_bw(), tooltip = "text" )
head(Gene_PA_meta)

d <- highlight_key(info.norm_umap , ~plasmid)   

dataplot <- ggplotly(ggplot(d, 
                     aes(x = umap1, y = umap2, colour = CP_Gene, text = paste("Plasmid_cluster:", Plasmid_cluster, "<br>Lab_Species:",Lab_Species, "<br>ST:",ST, "<br>Length:",Length))) + 
                geom_point(alpha = 0.3) + theme_bw(), tooltip = "text" ) %>%
  highlight("plotly_selected", dynamic = TRUE)

#options(persistent = FALSE)
# 
dataplot_dt <- datatable(d, extensions = 'Buttons', options = list(
  dom = 'Blfrtip',
  buttons = c('copy', 'csv', 'excel', 'pdf'),
  lengthMenu = list(c(10,30, 50, -1), 
                    c('10', '30', '50', 'All')),
  paging = T) )
  

bscols(widths = c(8,1),  list(dataplot, dataplot_dt))

# info.norm_umap %>% 
#   filter(CP_Gene=="blaKPC2") %>% 
#   ggplot(aes(x = umap1, y = umap2, colour = factor(Hospital))) + 
#   geom_point(alpha = 0.3) + theme_bw()

# umap_hos <- info.norm_umap %>% 
#   filter(CP_Gene=="blaNDM1") %>% 
#   ggplot(aes(x = umap1, y = umap2, colour = factor(Hospital))) + 
#   geom_point(alpha = 0.3) + theme_bw()

# umap_len <- info.norm_umap %>% 
#   filter(CP_Gene=="blaNDM1") %>% 
#   ggplot(aes(x = umap1, y = umap2, colour = Length)) + 
#   geom_point(alpha = 0.3) + theme_bw()




