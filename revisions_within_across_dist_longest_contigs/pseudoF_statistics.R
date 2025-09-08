

require(ggplot2); require(scales); require(reshape2); 
#install.packages("dplyr")
require("tidyverse")
require(dplyr)
#require(Hmisc)
library("readxl")
library(RColorBrewer)
library("ggsci")
#install.packages("ggrepel")
library("ggrepel")
library(ggpubr)


library(stringr)

library(dplyr)

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


dist_emb=read.csv('k7_v39_8k_s28_clade_All_TOL_Contigs_Global_logest/distances_kf2vec_emb_R.txt',sep="\t",header=TRUE)
dist_pl=read.csv('k7_v39_8k_s28_clade_All_TOL_Contigs_Global_logest/distances_kf2vec_placement_R.txt',sep="\t",header=TRUE)
dist_fastme=read.csv('fastme/distances_kf2vec_fastme_R.txt',sep="\t",header=TRUE)

dist_d2star=read.csv('cafe_dist/pairwise_distances_cafe_d2star_phylip_R.txt',sep="\t",header=TRUE)
dist_d2shepp=read.csv('cafe_dist/pairwise_distances_cafe_d2shepp_phylip_R.txt',sep="\t",header=TRUE)
dist_cvtree=read.csv('cafe_dist/pairwise_distances_cafe_cvtree_phylip_R.txt',sep="\t",header=TRUE)
dist_cosine=read.csv('cafe_dist/pairwise_distances_cafe_cosine_phylip_R.txt',sep="\t",header=TRUE)

dist_cophylog=read.csv('cafe_dist/pairwise_distances_cafe_co-phylog_phylip_R.txt',sep="\t",header=TRUE)
dist_eu=read.csv('cafe_dist/pairwise_distances_cafe_eu_phylip_R.txt',sep="\t",header=TRUE)
dist_js=read.csv('cafe_dist/pairwise_distances_cafe_js_phylip_R.txt',sep="\t",header=TRUE)
dist_ma=read.csv('cafe_dist/pairwise_distances_cafe_ma_phylip_R.txt',sep="\t",header=TRUE)

#dist_all <- rbind(dist_emb, dist_pl, dist_fastme)
#colnames(dist_all)

dfs <- list(dist_emb, dist_pl, dist_fastme, dist_d2star, dist_d2shepp, dist_cvtree, dist_cosine, dist_cophylog, dist_eu, dist_js, dist_ma)


# 1. Find common (Contig1, Contig2) pairs across all dataframes
common_pairs <- Reduce(function(x, y) inner_join(x, y, by = c("Contig1","Contig2")),
                       lapply(dfs, function(df) df %>% select(Contig1, Contig2)))

# 2. Sample 30% of the common pairs
set.seed(123)  # reproducible
sampled_pairs <- sample_frac(common_pairs, 1)

# 3. Subset all dataframes to the sampled pairs
dfs_sub <- lapply(dfs, function(df) semi_join(df, sampled_pairs, by = c("Contig1","Contig2")))


combined_df <- bind_rows(dfs_sub)

# Example: access subsetted dataframes
#df1_sub <- dfs_sub[[1]]
#df2_sub <- dfs_sub[[2]]


ggplot(aes(x=Condition, y=Distance_scaled, color = Type),
       data=combined_df)+
  #stat_summary(geom="crossbar")+
  #stat_summary(geom="point")+
  #geom_violin()+
  stat_summary()+
  #geom_boxplot() +
  theme_classic()
  



summarise

combined_df %>% 
  distinct(Condition, Type)

combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD=median(Distance)) %>% 
 # head
  pivot_wider(names_from = "Type",values_from = "mD" ) %>%
  mutate(pF=Across/Within)

head(combined_df)

df_summary


# Table included in a paper
df_summary <- combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD=mean(Distance))%>%
  pivot_wider(names_from = "Type", values_from = "mD" ) %>%
  mutate(pF=Across/Within)%>%
  mutate(Across = if_else(grepl("kf2vec", Condition), Across / 100, Across))%>%
mutate(Within = if_else(grepl("kf2vec", Condition), Within / 100, Within))%>%
# Round numbers for display
mutate(
  Across   = round(Across, 3),
  Within   = round(Within, 3),
  pF = round(pF, 2)
)
#df_summary$pF <- format(df_summary$pF, scientific = FALSE, nsmall = 6)
#print(df_summary)

#df_summary <- df_summary %>%
#  mutate(Across = if_else(grepl("kf2vec", Condition), Across / 100, Across))

print(as.data.frame(df_summary))




# Table included in a paper (median)
df_summary <- combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD=median(Distance))%>%
  pivot_wider(names_from = "Type", values_from = "mD" ) %>%
  mutate(pF=Across/Within)%>%
  mutate(Across = if_else(grepl("kf2vec", Condition), Across / 100, Across))%>%
  mutate(Within = if_else(grepl("kf2vec", Condition), Within / 100, Within))%>%
  # Round numbers for display
  mutate(
    Across   = round(Across, 3),
    Within   = round(Within, 3),
    pF = round(pF, 2)
  )
#df_summary$pF <- format(df_summary$pF, scientific = FALSE, nsmall = 6)
#print(df_summary)

#df_summary <- df_summary %>%
#  mutate(Across = if_else(grepl("kf2vec", Condition), Across / 100, Across))

print(as.data.frame(df_summary))




pF = combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD = round(mean(Distance, na.rm = TRUE), 4)) %>%
  pivot_wider(names_from = "Type", values_from = "mD") %>%
  mutate(pF = Across/Within)

combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD = round(mean(Distance, na.rm = TRUE), 4)) #%>%
  pivot_wider(names_from = "Type", values_from = "mD") %>%
  mutate(pF = Across/Within)


combined_df %>%
  group_by(Condition, Type) %>%
  dplyr::summarise(mD=mean(Distance))%>%  # compute mean
  print(n = 22)   


stat_summary(fun.data = "mean_cl_boot",size = 1, alpha = 0.7)+
  coord_cartesian(ylim=c(0,10))
#theme_classic()

combined_df2 = combined_df
combined_df2$Contig2 = combined_df$Contig1
combined_df2$Contig1 = combined_df$Contig2


sibling = rbind(combined_df,combined_df2) %>%
  group_by(Contig1,Condition) %>%
  #filter(Condition=="cafe_phylip_js") %>%
  dplyr::summarise(cl=Type[which.min(Distance)],clw=sum(Type=="Within")>0) %>% 
  filter(clw==TRUE) %>%
  group_by(Condition,cl) %>%
  dplyr::summarise(n=n())%>%  
  pivot_wider(names_from = "cl", values_from = "n") %>%
  mutate(errorrate = Across/(Across+Within),n=Across+Within)

sibling%>%
  print(n = 22) 

#pal <- colorRampPalette(brewer.pal(8, "Set2"))  # make it flexible
#my_colors <- pal(12) 
set2_colors <- brewer.pal(8, "Set2")
set3_colors <- brewer.pal(12, "Set3")  # 12 colors, distinct from Set2
my_colors <- c(set2_colors, set3_colors[8])


dark2_colors <- brewer.pal(8, "Dark2")  # max 8 colors

merge(sibling,pF,by="Condition") %>%
  mutate(pFi=1/pF) %>%
  select(Condition,errorrate,pFi) %>%
  filter(Condition != "cafe_phylip_co-phylog") %>%
  filter(Condition != "kf2vec_after_placement") %>%
  pivot_longer(cols=2:3) %>%
  #ggplot(aes(x=reorder(Condition,value),y=value,fill=name))+
  #  geom_bar(stat="identity",position = position_dodge())
  ggplot(aes(x=name,y=value,fill=Condition))+
  geom_bar(stat="identity",position = position_dodge(),color="black")+
  #scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = my_colors) +
  theme_classic()


# Paper plot
merge(sibling,pF,by="Condition") %>%
  mutate(pFi=1/pF) %>%
  select(Condition,errorrate,pFi) %>%
  filter(Condition != "cafe_phylip_co-phylog") %>%
  filter(Condition != "kf2vec_after_placement") %>%
  pivot_longer(cols=2:3) %>%
#ggplot(aes(x=reorder(Condition,value),y=value,fill=name))+
#  geom_bar(stat="identity",position = position_dodge())
  filter(name == "errorrate") %>%
ggplot(aes(x=reorder(Condition, -value),y=value,fill=Condition))+
  geom_bar(stat="identity",position = position_dodge(),color="black")+
  #scale_fill_brewer(palette = "Set3")+
  #scale_fill_manual(values = my_colors)+ 
  scale_fill_manual(values = c(dark2_colors[7], dark2_colors[6], dark2_colors[6], dark2_colors[6], dark2_colors[7], dark2_colors[7], dark2_colors[7], dark2_colors[5], dark2_colors[5]))+ 
  scale_x_discrete( labels = c('CVtree', 'JS', 'Ma', 'Eu',expression(italic("D")[italic("2")]^italic("*")), 'Cosine',  expression(italic("D")[italic("2")]^italic("S")), 'kf2vec', expression(atop("kf2vec/", "FastME2")) )) +
  theme_classic()+
  xlab(NULL) +
  scale_y_continuous("Closest contig mismatch",labels = percent)+
  theme(legend.position = "none")
  #theme(
  #  axis.text.x = element_text(lineheight = 0.3)  # only affects x-axis labels
  #)# 
ggsave("contigs_across_within_bar_D3.pdf",width=4.8,height = 4.0)



# Annotated

merge(sibling,pF,by="Condition") %>%
  mutate(pFi = 1/pF) %>%
  select(Condition, errorrate, pFi) %>%
  filter(Condition != "cafe_phylip_co-phylog") %>%
  filter(Condition != "kf2vec_after_placement") %>%
  pivot_longer(cols=2:3) %>%
  filter(name == "errorrate") %>%
  ggplot(aes(x=reorder(Condition, -value), y=value, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +
  geom_text(aes(label = scales::percent(value, accuracy = 0.1)), 
            vjust = -0.3, size=3.5) +  # annotate with percent values
  scale_fill_manual(values = c(
    dark2_colors[7], dark2_colors[6], dark2_colors[6], dark2_colors[6],
    dark2_colors[7], dark2_colors[7], dark2_colors[7],
    dark2_colors[5], dark2_colors[5]
  )) +
  scale_x_discrete(labels = c(
    'CVtree', 'JS', 'Ma', 'Eu',
    expression(italic("D")[italic("2")]^italic("*")),
    'Cosine',
    expression(italic("D")[italic("2")]^italic("S")),
    'kf2vec',
    expression(atop("kf2vec/", "FastME2"))
  )) +
  theme_classic() +
  xlab(NULL) +
  scale_y_continuous("Closest contig mismatch", labels = scales::percent) +
  theme(legend.position = "none")



#scale_x_discrete(labels = c('Cosine', 'Eu', 'Ma',"CVtree", expression(italic("D")[italic("2")]^italic("*")), expression(italic("D")[italic("2")]^italic("S")), 'kf2vec' ))+
  
getwd()

rbind(combined_df,combined_df2) %>%
    filter(Contig1=="G000007705.part_NC_005085.1",Condition=="cafe_phylip_js")

  dplur::summarise(mD = round(mean(Distance, na.rm = TRUE), 4)) %>%
  pivot_wider(names_from = "Type", values_from = "mD") %>%
  mutate(pF = Across/Within)
  
  
# sed '$!N;s/\n/,/' kmer_pl_3layer_lowerLrRate.out | awk 'BEGIN { FS = "[, ]" } ; {print $9, $21}' > kmer_pl_3layer_lowerLrRate.csv
