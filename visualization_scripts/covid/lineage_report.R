library(tidyverse)
library(dplyr)
library(ggplot2)

covid.df <- read.csv("/home/prakki/Documents/LeaRn/Covid/lineage_report_forR.csv")
head(covid.df)


pdf("/home/prakki/Documents/LeaRn/Covid/covid_report.pdf", width=10, height=7)

# Strain counts

positions <- c("Covid Negative", "No template control", "Not assigned" , "Delta (B.1.617.2-like)", 
               "Omicron (BA.1-like)" , "Omicron (BA.2-like)")

covid.df %>% ggplot(aes(x=factor(scorpio_call))) +
  geom_bar(stat="count", width=0.1, fill="#00AFBB") +
  geom_text(stat='count',aes(label=..count..),vjust=-1,size=3) +
  theme_bw() + ggtitle("Strain Counts") + # for the main title
  xlab("Strain") + # for the x axis label
  ylab("Count") + scale_x_discrete(limits = positions)


# Strains vs Lineage

covid.df %>% ggplot(aes(x=factor(lineage))) +
  geom_bar(stat="count", width=0.1, fill="#00AFBB") +
  geom_text(stat='count',aes(label=..count..),vjust=-1,size=3) +
  theme_bw() + ggtitle("Lineage Counts by Strain") + # for the main title
  xlab("Strain and Lineage") + # for the x axis label
  ylab("Count") + facet_wrap(~scorpio_call, scales = "free_x") +
  theme(strip.text = element_text(face = "plain", color = "black", size = 10),
        strip.background = element_rect(fill = "#FFB310")) 


# Coverage vs Percent Genome Coverage >100x

covid.df %>%  ggplot(aes(Coverage, Prcnt_of_Genome_Coverage_with_gt100x, color = scorpio_call)) +
  geom_point(size=2) +
  labs(
    x = "Coverage of the genome",
    y = "Percentage of genome with coverage > 100x",
    title = "Coverage vs Percent Genome Coverage",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
   ) +
   annotate("text",label = stringr::str_wrap("* Samples with sequencing depth of 250-300x had at least 90% of the genome covered >100x", 20),
            x = c(300), y = c(0.5), color = "#D55E00", size = 4) + # #D55E00,#009E73
  annotate("text",label = stringr::str_wrap("* Most samples with sequencing depth < 100x were unassigned variants", 20),
           x = c(300), y = c(0.25), color = "#D55E00", size = 4) # #D55E00,#009E73

# Coverage vs Nanodrop (ng/ul)

covid.df %>%  ggplot(aes(Nanodrop.ng_per_ul., Coverage, color = scorpio_call)) +
  geom_point(size=2) +
  labs(
    x = "Nanodrop.ng_per_ul.",
    y = "Coverage of the genome",
    title = "Coverage vs Nanodrop (ng/ul)",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) +
  annotate("text",label = stringr::str_wrap("A weak positive correlation between nanodrop and Sequencing Coverage", 20),
           x = c(20), y = c(100), color = "#D55E00", size = 5) # #D55E00,#009E73


# Ct value vs Coverage

covid.df %>%  ggplot(aes( CtValue, Coverage, color = scorpio_call)) +
  geom_point(size=2) +
  labs(
    x = "Ct value",
    y = "Coverage of the genome",
    title = "Ct value vs Coverage",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) +
  annotate("text",label = stringr::str_wrap("A weak negative correlation between Ctvalue and Sequencing coverage", 20),
           x = c(27), y = c(120), color = "#D55E00", size = 5) # #D55E00,#009E73

## Strains vs SNPs

covid.df %>%  ggplot(aes(fct_reorder(scorpio_call, SNPs), SNPs, color = scorpio_call)) +
  #geom_violin(width=0.2,lwd=1, trim=TRUE) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(
    x = "Strain",
    y = "SNPs",
    title = "Strains vs SNPs",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) +
  annotate("text",label = stringr::str_wrap("Low number of SNPs may be due to high Missing bases (Ns)", 20),
           x = c("Delta (B.1.617.2-like)"), y = c(18), color = "#D55E00", size = 5) # #D55E00,#009E73


# Strains vs Missingness

covid.df %>%  ggplot(aes(fct_reorder(scorpio_call, Ns), Ns, color = scorpio_call)) +
  #geom_violin(width=0.2,lwd=1, trim=TRUE) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(
    x = "Strain",
    y = "Missing bases (Ns)",
    title = "Strains vs Missing bases (Ns)",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) +
  annotate("text",label = stringr::str_wrap("Due to higher missingness, strains and lineage were not were assigned", 20),
           x = c("No template control"), y = c(25000), color = "#D55E00", size = 4) # #D55E00,#009E73

# ambiguity_score vs Ns

covid.df %>%  ggplot(aes(Ns, ambiguity_score, color = scorpio_call)) +
  geom_point(size=2) +
  labs(
    x = "Missing bases (Ns)",
    y = "Ambiguity score",
    title = "Missing bases (Ns) vs Ambiguity score",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) +
  annotate("text",label = stringr::str_wrap("Due to higher missingness, ambiguity score (confidence) become low", 20),
           x = c(20000), y = c(0.9), color = "#D55E00", size = 5) # #D55E00,#009E73


# HS_qubit (ng/ul) vs Coverage
head(covid.df)

covid.df %>%  ggplot(aes(HS_qubit.ng_per_ul., Coverage, color = scorpio_call)) +
  geom_point(size=2) +
  geom_smooth(method = "lm", se = FALSE)
  labs(
    x = "HS_qubit (ng/ul)",
    y = "Coverage",
    title = "HS_qubit (ng/ul) vs Coverage",
    subtitle = "",
  ) +
  scale_color_manual(values = c("Covid Negative" = "#7eb37a", "No template control" = "#7eb37a", "Not assigned" = "#FFB310", "Delta (B.1.617.2-like)" = "#00AFBB", 
                                "Omicron (BA.1-like)" = "navyblue", "Omicron (BA.2-like)" = "red")) +
  #scale_fill_manual(values = c25) + #Paired and Set3
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      size = 7
    ),
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0),
    strip.text.x.top = element_text(angle = 90)
  ) 
  
  
  qubit_cov.lm <- lm(HS_qubit.ng_per_ul. ~ Coverage, covid.df)
  #qubit_cov.lm
  ggplot(covid.df, aes(x = HS_qubit.ng_per_ul., y = Coverage)) + 
    geom_point() + 
    geom_abline(slope = coef(qubit_cov.lm )[["Coverage"]], 
                intercept = coef(qubit_cov.lm )[["(Intercept)"]])
  
  
dev.off()
