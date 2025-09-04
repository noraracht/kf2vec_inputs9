

require(ggplot2); require(scales); require(reshape2); 
#install.packages("dplyr")
require(dplyr)
#require(Hmisc)
library("readxl")
library(RColorBrewer)
library("ggsci")
#install.packages("ggrepel")
library("ggrepel")
library(ggpubr)


library(stringr)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()


library(RColorBrewer)
my_palette = c(brewer.pal(9, "RdBu")[c(1,2, 3, 7, 9)])
my_palette_lst =as.list(strsplit(my_palette, " "))


dark2_colors <-brewer.pal(n = 8, name = "Dark2")
dark2_colors

set1_colors <-brewer.pal(n = 8, name = "Set1")
set1_colors

accent_colors <-brewer.pal(n = 8, name = "Accent")
accent_colors


amber_df=read.csv(file.path("results", "TOL_full_genome_queries_claded_fullbackbone", "all.pl_err"),sep="\t",header=T)

# Comparison of full genome 500 queries on different models

# Claded, full, classified
df_CFC =read.csv(file.path("results", "TOL_full_genome_queries_claded_fullbackbone", "all.pl_err"),sep=" ",header=FALSE)
df_CFC ["claded"] = "Claded"
df_CFC ["backbone"] = "Full"
df_CFC ["clsf"] = "cmp"

# Claded, full, TRUE
df_CFT =read.csv(file.path("results", "TOL_full_genome_queries_claded_fullbackbone_TRUE", "all.pl_err"),sep=" ",header=FALSE)
df_CFT ["claded"] = "Claded"
df_CFT ["backbone"] = "Full"
df_CFT ["clsf"] = "tr"


# Claded, high quality, classified
df_CHC =read.csv(file.path("results", "TOL_full_genome_queries_claded_hiqualbackbone", "all.pl_err"),sep=" ",header=FALSE)
df_CHC ["claded"] = "Claded"
df_CHC ["backbone"] = "HighQ"
df_CHC ["clsf"] = "cmp"

# Claded, high quality, TRUE
df_CHT =read.csv(file.path("results", "TOL_full_genome_queries_claded_hiqualbackbone_TRUE", "all.pl_err"),sep=" ",header=FALSE)
df_CHT ["claded"] = "Claded"
df_CHT ["backbone"] = "HighQ"
df_CHT ["clsf"] = "tr"



# Uncladed, full, na
df_UF =read.csv(file.path("results", "TOL_full_genome_queries_uncladed_fullbackbone", "all.pl_err"),sep=" ",header=FALSE)
df_UF ["claded"] = "Uncladed"
df_UF ["backbone"] = "Full"
df_UF ["clsf"] = "na"

# Uncladed, high quality, TRUE
df_UH =read.csv(file.path("results", "TOL_full_genome_queries_uncladed_hiqualbackbone", "all.pl_err"),sep=" ",header=FALSE)
df_UH ["claded"] = "Uncladed"
df_UH ["backbone"] = "HighQ"
df_UH ["clsf"] = "na"




df_inter = rbind(df_CFC, df_CFT, df_CHC, df_CHT, df_UF, df_UH)

head(df_inter)

df_inter
# Unify format for all dataframes
df_inter$V1 <- gsub('.pl_err', '', df_inter$V1)
df_inter$V1 <- gsub('apples_input_di_mtrx_query_', '', df_inter$V1)
df_inter <- df_inter[,colSums(is.na(df_inter))<nrow(df_inter)]
colnames(df_inter)[colnames(df_inter)=="V3"] <- "V2"
colnames(df_inter)[colnames(df_inter)=="V4"] <- "V3"
colnames(df_inter)[colnames(df_inter)=="V5"] <- "V4"
head(df_inter)

df_inter = df_inter[!(is.na(df_inter$V2)), ]


nrow(df_inter)
head(df_inter)
df_inter

df_fin = df_inter
df_fin$condition <- paste(df_fin$claded, "-", df_fin$backbone, "-", df_fin$clsf)

head(df_fin)

q = read.csv('/Users/nora/Documents/ml_metagenomics/tol_quality_scores/quality_comparison_hor.csv')
head(q)

#nov = read.csv('/Users/nora/Documents/ml_metagenomics/pendant.txt',sep="\t",h=F)
#nov$V2 = nov$V2*2/100
#head(nov)

nov = read.csv('/Users/nora/Documents/ml_metagenomics/closest_dev_set_2col.txt',sep="\t",h=FALSE)
#head(nov)
nov$V2 = nov$V2*1/100
head(nov)

names(nov)[2] = "nov"
df_fin = merge(df_fin,nov,by="V1")
head(df_fin)

qdf_fin = merge(q,df_fin, by.x="Assembly", by.y = "V1") 

head(qdf_fin)

ggplot(aes(x=condition, y=V3),
       data=df_fin)+
  #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
  geom_boxplot()+
  #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
  #stat_summary(aes(fill=
  #                   ifelse(claded =="Claded","Claded","Uncladed")),geom="bar", )+
 #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", fun = "mean", vjust = 2.5, colour = "black" ) +
  
  #facet_wrap(.~chunked)
  stat_summary()+
  #geom_bar( stat="identity", fill="skyblue", alpha=0.7)
  stat_summary(geom="point")+
  #geom_violin()+
  theme_classic()
  

  
  quantile(nov$nov,(0:10)/10)
  qdf_fin 
  nrow(qdf_fin)

  qdf_fin <- qdf_fin %>%
    mutate(facet_label = sub("-.*-", "", condition)) %>%
    mutate(facet_label = recode_factor(facet_label,
                                       "Claded  cmp" = "Cladded",
                                       "Claded  tr"  = "Cladded (true)",
                                       "Uncladed  na"= "Uncladded"))
  
  # Set the order explicitly
  qdf_fin$facet_label <- factor(qdf_fin$facet_label, 
                                levels = c("Uncladded", "Cladded","Cladded (true)"))


  unique(qdf_fin$facet_label)
  
  ggplot(aes(y=V3, x=cut(nov,breaks=c(0,0.01,0.05,0.15,0.2,0.5,2)), 
          color = grepl("HighQ",condition )),
        # color = condition ),
         data=qdf_fin)+
    #data=df_fin[df_fin$chunked %in% c("Unchunked"),])+
    #geom_boxplot(outlier.shape=NA)+
    #geom_boxplot(fill="skyblue1", outlier.shape=NA)+
    stat_summary(aes(group=condition), geom="line" )+
    #stat_summary(aes(group=condition),geom="bar" ,position = position_dodge(0.9),color="black")+
    #stat_summary(position = position_dodge(0.9),color="black")+
    #facet_wrap(~sub("-.*-","",condition))+
    facet_wrap(~facet_label)+
    stat_summary()+
    theme_classic()+
    scale_fill_brewer(palette = "Paired",name =NULL)+
  scale_color_manual(labels = c("Full","High-quality"), name = NULL,
                      values=c("#a6cee3","#1f78b4","#33a02c","red"))+
    #theme( axis.text.x = element_text(angle = 10, hjust = 1))+
    theme(axis.text.x = element_text( angle = 30, hjust = 1, vjust = 1))+
    #theme(axis.text.x = element_text( size = 8))+
    theme(legend.position = c(0.1,0.88),
          legend.margin = margin(0, 0, 0, 0))+
    scale_y_continuous("Placement error")+scale_x_discrete("Query novelty")
  
  ggsave("Full-HighQ-D1-line.pdf",width = 6.5,height = 3.5)
  getwd()
  
    #geom_text(data = df_fin$V3, aes(label = sprintf("%.1f", round(V3, digits = 2)), y = V3 - 0.45))+
    #geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 2)))), stat = "summary", 
    #          fun = "mean", vjust = 2.5) +
    
    #facet_wrap(.~chunked)+



