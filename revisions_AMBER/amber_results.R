

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
require(plyr)
require(tidyr)



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


################################################################
# Plot AMBER results for gsa short
#folder = "amber_closest_subset_5kfilt"

#folder = "amber_closest_subset_hiqual" # also includes >10k cut off

#folder = "amber_closest_subset_5kfilt_above5k"
#folder = "amber_closest_subset_5kfilt_above10k"
folder = "amber_closest_subset_5kfilt_above10k_revisions"
#folder = "amber_closest_subset_5kfilt_above20k"

#folder = "amber_closest_distance_long_reads_S0_hiqual"
#folder = "amber_closest_distance_long_reads_S0_above10k"

#folder = "amber_closest_subset_combo"
#folder = "amber_closest_subset_hiqual_ALL"

amber_df=read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "results.tsv"),sep="\t",header=T)
amber_df <- subset(amber_df, Tool != "closest_subset_v10_revtax")

nrow(amber_df)
head(amber_df)

colnames(amber_df)
unique(amber_df$Tool)
head(amber_df$accuracy_bp..unfiltered.ool)

amber_df$rank  <- factor(amber_df$rank, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

amber_df$Tool

ggplot(aes(x=rank, y=accuracy_seq, fill = Tool,color = Tool, group = Tool),
       data=amber_df[(!grepl("losest",amber_df$Tool) |grepl("v10",amber_df$Tool)) & 
                       !amber_df$Tool %in% c("Closest_tree_dist_phyl","tax2tree", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard", "closest_subset_v2", "closest_subset_v3"),])+
      #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line",  alpha=0.7)+
  #stat_summary(geom="bar" , alpha = 0.7,position = position_dodge())+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  coord_cartesian(ylim=c(0.1, 1.0))+
  #scale_color_brewer(name = "",palette="Set1")
  #scale_color_brewer(name = "",palette="Dark2")
  scale_fill_brewer(name = "",palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  scale_color_brewer(name = "",palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="", label = c("Superkindom", "Phylum", "Class","Order", "Family", "Genus", "Species"))  +
  ylab("Accuracy")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust=0.95))+
theme(legend.position = c(0.3,0.4),legend.margin=margin(0,-2,0,0)
      #,
      #axis.text.x = element_text(size = 8)
)

nrow(amber_df)

head(amber_df)

ggplot(aes(x=rank, y=accuracy_seq, fill = reorder(Tool,accuracy_seq)),
       data=amber_df[(!grepl("losest",amber_df$Tool) |grepl("v10",amber_df$Tool)) & 
                       !amber_df$Tool %in% c("Closest_tree_dist_phyl","tax2tree", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard", "closest_subset_v2", "closest_subset_v3"),])+
  #stat_summary(geom="line",  alpha=0.7)+
  stat_summary(geom="bar",  position = position_dodge(width=0.8),color="black",width=0.75)+
  theme_classic()+
  coord_cartesian(ylim=c(0.1, 1.0))+
  #scale_color_brewer(name = "",palette="Set1")
  #scale_color_brewer(name = "",palette="Dark2")
  
  #scale_fill_brewer(name = "", palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_color_brewer(name = "", palette="Dark2", labels=c("closest_subset_v10" = "kf2d"))+
  
  #scale_fill_manual(values = c(dark2_colors[2], dark2_colors[4], dark2_colors[5], dark2_colors[3],
  #                             dark2_colors[1], dark2_colors[6],dark2_colors[7], dark2_colors[8]),
  #                  name = "", labels=c("closest_subset_v10" = "kf2d"))+
  #scale_color_manual(values = c(dark2_colors[2], dark2_colors[4], dark2_colors[5], dark2_colors[3],
  #                              dark2_colors[1], dark2_colors[6],dark2_colors[7], dark2_colors[8]),
  #                  name = "", labels=c("closest_subset_v10" = "kf2d"))+
  
 # scale_color_manual(values = c(my_set1_colors[8], my_set1_colors[4], my_set1_colors[2], my_set1_colors[3], my_set1_colors[1],
#                                my_set1_colors[5], my_set1_colors[6],my_set1_colors[7] ),
#                     labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2d", "Kraken 2.0.8-beta" ), name = "")+
#  scale_fill_manual(values = c(my_set1_colors[8], my_set1_colors[4], my_set1_colors[2], my_set1_colors[3], my_set1_colors[1],
#                                my_set1_colors[5], my_set1_colors[6],my_set1_colors[7] ),
 #                    labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2d", "Kraken 2.0.8-beta" ), name = "")+
 
  scale_color_manual(values = c(#9BA8AEFF", dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                  dark2_colors[2], dark2_colors[6],dark2_colors[8], dark2_colors[5]),
                     labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2vec", "kraken_customDB"="Kraken*", "Kraken 2.0.8-beta" ), name = NULL)+
  scale_fill_manual(values = c("#9BA8AEFF", dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                 dark2_colors[2], dark2_colors[6],dark2_colors[8], dark2_colors[5] ),
                    labels= c("DIAMOND 0.9.28","MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "closest_subset_v10" = "kf2vec", "kraken_customDB"="Kraken*", "Kraken 2.0.8-beta" ), name = NULL)+
  
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="", label = c("Superkindom", "Phylum", "Class","Order", "Family", "Genus", "Species"))  +
  ylab("Accuracy")+
  theme(axis.text.x = element_text(angle = 25, vjust = 1.1, hjust=1))+
  theme(#legend.position="bottom",
    legend.position=c(.48,-0.18),
        legend.margin=margin(0, 0,0,0), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.5, "lines")
        #axis.text.x = element_text(size = 4)
        )+
#guides(color=guide_legend(ncol=2))
guides(fill=guide_legend(ncol=3))


ggsave("CAMI-D5-accuracy-line-revisions.pdf",width=3.8,height = 4.0)
#ggsave("/Users/nora/Documents/CV/IMSI/CAMI-D5-accuracy-line-poster.pdf",width=4.8,height = 4.0)
getwd()



unique(amber_df$Tool)

stats_data = amber_df[amber_df$Tool=="closest_subset_v10",]
stats_data2 = amber_df[amber_df$Tool=="LSHVec cami2",]
stats_data$accuracy_seq
stats_data2$accuracy_seq
summary(c(stats_data$accuracy_seq[-length(stats_data$accuracy_seq)]))
summary(c(stats_data2$accuracy_seq[-length(stats_data2$accuracy_seq)]))


  ggplot(aes(x=rank, y=accuracy_seq, color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()+
    coord_cartesian(ylim=c(0.55, 1.0))
  #ylim(0.75, 1.05)
  scale_color_brewer(name = "",palette="Set1")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  

  ggplot(aes(x=rank, y=accuracy_seq..unfiltered., color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
    #ylim(0.75, 1.05)
    scale_color_brewer(name = "",palette="Dark2")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  
  
  colnames(amber_df)
  
  ggplot(aes(x=rank, y=f1_score_per_seq..unfiltered., color = Tool, group = Tool),
         data=amber_df[!amber_df$Tool %in% c("Closest_tree_dist_phyl", "Closest_tree_dist_cls", "Closest_tree_dist_ordr", "Closest_tree_dist_fam", "Closest_tree_dist_gen", "Gold standard"),])+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    stat_summary(geom="line" , alpha = 0.7)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
    scale_color_brewer(name = "",palette="Dark2")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50))
  #scale_x_discrete(name="Contig length (KB)", label = c("(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", bquote( '(2560-'  *10^4* ']' )))  +
  scale_x_discrete(name="Contig length (KB)", label = c("(5-10]", "(10-20]", "(20-40]","(40-80]", "(80-160]", "(160-320]", "(320-640]", "(640-1280]", "(1280-2560]", "(2560-10000]"))  +
    ylab("Placement error")+
    theme(axis.text.x = element_text(angle = 10, vjust = 0.7, hjust=0.5))
  


  

  
#df1 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/Closest_tree_dist_sp/metrics_per_bin.tsv"),sep="\t",header=T)
df2 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/DIAMOND 0.9.28/metrics_per_bin.tsv"),sep="\t",header=T)
df3 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/Kraken 2.0.8-beta/metrics_per_bin.tsv"),sep="\t",header=T)
df4 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/LSHVec cami2/metrics_per_bin.tsv"),sep="\t",header=T)
df5 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/MEGAN 6.15.2/metrics_per_bin.tsv"),sep="\t",header=T)
df6 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/PhyloPythiaS+ 1.4/metrics_per_bin.tsv"),sep="\t",header=T)
#df7 = read.csv( file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset/metrics_per_bin.tsv"),sep="\t",header=T)
#df8 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v2/metrics_per_bin.tsv"),sep="\t",header=T)
#df9 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/tax2tree/metrics_per_bin.tsv"),sep="\t",header=T)
#df10 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v3/metrics_per_bin.tsv"),sep="\t",header=T)
#df11 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v4/metrics_per_bin.tsv"),sep="\t",header=T)
#df12 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v5/metrics_per_bin.tsv"),sep="\t",header=T)
#df13 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v6/metrics_per_bin.tsv"),sep="\t",header=T)
#df14 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v7/metrics_per_bin.tsv"),sep="\t",header=T)
#df15 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v8/metrics_per_bin.tsv"),sep="\t",header=T)
df16 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v9/metrics_per_bin.tsv"),sep="\t",header=T)
df17 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v10/metrics_per_bin.tsv"),sep="\t",header=T)
df18 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/kraken_customDB/metrics_per_bin.tsv"),sep="\t",header=T)
#df19 = read.csv(file.path("/Users/nora/Documents/ml_metagenomics/tol_variable_k_resuls", folder, "taxonomic/closest_subset_v10_revtax/metrics_per_bin.tsv"),sep="\t",header=T)

#df1["Tool"] = "Closest_tree_dist_sp"
df2["Tool"] = "DIAMOND 0.9.28"
df3["Tool"] = "Kraken 2.0.8-beta"
df4["Tool"] = "LSHVec cami2"
df5["Tool"] = "MEGAN 6.15.2"
df6["Tool"] = "PhyloPythiaS+ 1.4"
#df7["Tool"] = "closest_subset"
#df8["Tool"] = "closest_subset_v2"
#df9["Tool"] = "tax2tree"
#df10["Tool"] = "closest_subset_v3"
#df11["Tool"] = "closest_subset_v4"
#df12["Tool"] = "closest_subset_v5"
#df13["Tool"] = "closest_subset_v6"
#df14["Tool"] = "closest_subset_v7"
#df15["Tool"] = "closest_subset_v8"
df16["Tool"] = "closest_subset_v9"
df17["Tool"] = "closest_subset_v10"
df18["Tool"] = "kraken_customDB"
#df19["Tool"] = "closest_subset_v10_revtax"

require(tidyr)

#install.packages("tidyverse")
require(tidyverse)

metrics_per_bin_df = rbind(df2, df3, df4, df5, df6, df16, df17, df18)
colnames (metrics_per_bin_df)
nrow(metrics_per_bin_df)
metrics_per_bin_df$Taxonomic.rank  <- factor(metrics_per_bin_df$Taxonomic.rank, levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))

metrics_per_bin_df%>%group_by(Taxonomic.rank, Tool, sample_id, Filtered) %>% mutate(Sample_Purity=mean(Purity..seq., na.rm = TRUE))

nrow(metrics_per_bin_df)

require(tidyverse)

my_set1_colors = RColorBrewer::brewer.pal(8, "Set1")
my_set1_colors[1]

unique(metrics_per_bin_df$Tool)




metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) %>%
ggplot(aes(x=purity/t,y=comp/t, color = reorder(Tool, purity), shape=Taxonomic.rank,
           group = Tool))+
  geom_path(size=0.1)+
  geom_point(size=2,fill=NA)+
  facet_wrap(.~ifelse(Taxonomic.rank=="species","species","higher ranks"),scales="free")+
  #scale_shape_manual(values=c(21,25,5,15,17,16,19), name = "", labels=c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))+
  scale_shape_manual(values=c(21,25,5,15,17,16,19), name = "")+
   #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values = c(dark2_colors[7], dark2_colors[5], dark2_colors[3], dark2_colors[4],
                                dark2_colors[2], dark2_colors[6],"#9BA8AEFF", dark2_colors[5]),
                     #breaks = tool_levels,
                     labels= c("MEGAN 6.15.2", "LSHVec cami2", "PhyloPythiaS+ 1.4", "kf2vec", "Kraken*", "Kraken 2.0.8-beta" ), name = "")+
  #coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  xlab ("Purity")+
  ylab ("Completeness")+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  guides(
    shape = guide_legend(order = 1, title = NULL, byrow = T),
    color = guide_legend(order = 2, title = NULL, byrow = T)
    
  ) +
  theme(
    legend.spacing.y = unit(1, "lines"),
    legend.margin = margin(0, 0, 0, 0), 
    legend.box = "vertical", 
    legend.title=element_blank()
  )
  #theme(legend.position = c(0.85,0.3), legend.margin=margin(0,0,0,0),  
  #     axis.text.x = element_text(size = 1)
  #)

getwd()
ggsave("CAMI-D5-purity_complt-line-revision.pdf",width=6.2,height = 4.0)



m <- metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) 
m$x = m$purity/m$t
m$y = m$comp/m$t

m[m$Taxonomic.rank=="species",]


metrics_per_bin_df %>%
  filter((!grepl("losest",Tool) | grepl("v10",Tool)) & Filtered == "False") %>%
  filter(Tool != "tax2tree") %>%
  filter(Tool != "DIAMOND 0.9.28") %>%
  #mutate(bin=cut(True.size..bp.,c(0,10^4,10^7,10^8.3,10^15))) %>%
  group_by(Taxonomic.rank,Tool) %>%
  dplyr::summarise(comp=sum(Completeness..seq.* True.size..bp.,na.rm = T),
                   purity=sum(Purity..seq.* True.size..bp.,na.rm = T),
                   t=sum( True.size..bp.),
                   n=n()) %>%
  ggplot(aes(x=purity/t,y=comp/t, color = reorder(Tool,purity), shape=Taxonomic.rank,
             group = Tool))+
  geom_path(size=0.1)+
  geom_point(size=2,fill=NA)+
  facet_wrap(.~ifelse(Taxonomic.rank=="species","species","higher ranks"),scales="free")+
  scale_shape_manual(values=c(21,25,5,15,17,16,19))+
  scale_color_brewer(palette = "Dark2")+
  #coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()

ggplot(aes(x=Taxonomic.rank, y=Purity..seq., color = Tool,fill=Tool, group = Tool),
       data=metrics_per_bin_df[
         (!grepl("losest",metrics_per_bin_df$Tool)|grepl("v10",metrics_per_bin_df$Tool)) &
           metrics_per_bin_df$Tool != "tax2tree"
         & metrics_per_bin_df$Filtered == "False",])+
  #data=metrics_per_bin_df)+
  #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  stat_summary()+
  stat_summary(geom="line", alpha=0.7)+
  #stat_summary(geom="bar" , alpha = 0.7, position = position_dodge())+
  coord_cartesian(ylim=c(0.,1))+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()
  #ylim(0.996, 1.0005)
  #scale_color_brewer(name = "",palette="Dark2")
  
  ggplot(aes(x=Taxonomic.rank, y=Completeness..seq., fill=Tool,color = Tool, group = Tool),
         data=metrics_per_bin_df[ (!grepl("losest",metrics_per_bin_df$Tool)|grepl("v10",metrics_per_bin_df$Tool)) &
                                    metrics_per_bin_df$Tool != "tax2tree" &
                                    metrics_per_bin_df$Filtered == "False",])+
    #data=metrics_per_bin_df)+
    #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
    #geom_violin()+
    stat_summary()+
    #stat_summary(geom="line", aes(group=conditions, color = conditions))
    #stat_summary(geom="line" , alpha = 0.7)+
    stat_summary(geom="line", alpha=0.7)+
    #stat_summary(geom="bar" , alpha = 0.7, position = position_dodge())+
    #geom_smooth(method=lm)+
    #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
    theme_classic()
  scale_color_brewer(name = "",palette="Dark2")

colnames(metrics_per_bin_df)

metrics_per_bin_df$Purity..seq.

metrics_per_bin_df["gp"] = 0

ggplot(aes(x=Purity..seq., y=Completeness..seq., color = Tool, group = gp),
       data=metrics_per_bin_df[metrics_per_bin_df$Filtered == "False" & metrics_per_bin_df$Taxonomic.rank == "species",])+
  #facet_grid(.~Taxonomic.rank)+
  #data=metrics_per_bin_df)+
  #data=amber_df[amber_df$Tool %in% c("closest_subset"),])+
  #geom_violin()+
  #stat_summary(fun.data = "mean_cl_boot", linewidth = 2, size = 0.1)
  stat_summary(fun = mean, geom = "point")

  #stat_summary(geom="line", aes(group=conditions, color = conditions))
  stat_summary(geom="line" , alpha = 0.7)+
  #geom_smooth(method=lm)+
  #stat_summary( geom="crossbar",fun.data = function(x) median_hilow_(x,ci=0.8))+
  theme_classic()+
  scale_color_brewer(name = "",palette="Dark2")



################################################################
################################################################
################################################################

