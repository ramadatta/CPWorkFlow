install.packages("scatterpie")
library(ggplot2)
library(scatterpie)
library(RColorBrewer)

my_df <- read.table("/data02/Analysis/for_Colleagues/for_zeqin/scatter_pie_plot/data.txt",sep = "\t", header = TRUE)
head(my_df)

xvals <- unique(my_df$Species)
xvals

yvals <- unique(my_df$Gene)
yvals

str(my_df)


my_df2 <- my_df %>% 
  mutate(Species_num = as.numeric(as.factor(Species)), 
         Gene_num = as.numeric(as.factor(Gene)))

my_df2

ggplot() + geom_scatterpie(data=my_df2, aes(x=Species_num, y=Gene_num, r =0.25), cols = colnames(my_df)[2:3]) + 
  scale_x_continuous(breaks=c(1,2,3),labels= c("Ecloc","KP","Ecoli")) +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8), labels=yvals) + 
  labs(x="Species", y="Gene") + 
  theme(panel.background  = element_blank()) +
  coord_fixed() 

# Make it more colorful 
Color<-brewer.pal(4, "Set2")

ggplot() + geom_scatterpie(data=my_df2, aes(x=Species_num, y=Gene_num, r =0.25), cols = colnames(my_df)[2:3]) +
  scale_fill_manual(name = "Sequence Type", values = Color) +
  scale_x_continuous(breaks=c(1,2,3),labels= c("Ecloc","KP","Ecoli")) +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8), labels=yvals) + 
  labs(x="Species", y="Gene") + 
  theme(panel.background  = element_blank()) +
  coord_fixed() 
