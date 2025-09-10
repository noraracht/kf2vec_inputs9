

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
library(data.table)

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

# Some entries are NAs
#dist_emb=read.csv('k7_v39_8k_s28_clade_All_TOL_Contigs_Global_logest/distances_kf2vec_emb_R.txt',sep="\t",header=TRUE)
#dist_pl=read.csv('k7_v39_8k_s28_clade_All_TOL_Contigs_Global_logest/distances_kf2vec_placement_R.txt',sep="\t",header=TRUE)
#dist_fastme=read.csv('fastme/distances_kf2vec_fastme_R.txt',sep="\t",header=TRUE)

#dist_d2star=read.csv('cafe_dist/pairwise_distances_cafe_d2star_phylip_R.txt',sep="\t",header=TRUE)
#dist_d2shepp=read.csv('cafe_dist/pairwise_distances_cafe_d2shepp_phylip_R.txt',sep="\t",header=TRUE)
#dist_cvtree=read.csv('cafe_dist/pairwise_distances_cafe_cvtree_phylip_R.txt',sep="\t",header=TRUE)
#dist_cosine=read.csv('cafe_dist/pairwise_distances_cafe_cosine_phylip_R.txt',sep="\t",header=TRUE)

#dist_cophylog=read.csv('cafe_dist/pairwise_distances_cafe_co-phylog_phylip_R.txt',sep="\t",header=TRUE)
#dist_eu=read.csv('cafe_dist/pairwise_distances_cafe_eu_phylip_R.txt',sep="\t",header=TRUE)
#dist_js=read.csv('cafe_dist/pairwise_distances_cafe_js_phylip_R.txt',sep="\t",header=TRUE)
#dist_ma=read.csv('cafe_dist/pairwise_distances_cafe_ma_phylip_R.txt',sep="\t",header=TRUE)



# Filtered subset include only contig pairs that are shared between samples
#dist_emb=read.csv('filtered/dist_emb_filtered.csv',sep="\t",header=TRUE)
#dist_pl=read.csv('filtered/dist_pl_filtered.csv',sep="\t",header=TRUE)
#dist_fastme=read.csv('filtered/dist_fastme_filtered.csv',sep="\t",header=TRUE)

#dist_d2star=read.csv('filtered/dist_d2star_filtered.csv',sep="\t",header=TRUE)
#dist_d2shepp=read.csv('filtered/dist_d2shepp_filtered.csv',sep="\t",header=TRUE)
#dist_cvtree=read.csv('filtered/dist_cvtree_filtered.csv',sep="\t",header=TRUE)
#dist_cosine=read.csv('filtered/dist_cosine_filtered.csv',sep="\t",header=TRUE)

#dist_cophylog=read.csv('filtered/dist_cophylog_filtered.csv',sep="\t",header=TRUE)
#dist_eu=read.csv('filtered/dist_eu_filtered.csv',sep="\t",header=TRUE)
#dist_js=read.csv('filtered/dist_js_filtered.csv',sep="\t",header=TRUE)
#dist_ma=read.csv('filtered/dist_ma_filtered.csv',sep="\t",header=TRUE)



 
#dfs <- list(dist_emb, dist_pl, dist_fastme, dist_d2star, dist_d2shepp, dist_cvtree, dist_cosine, dist_eu, dist_js, dist_ma)

#combined_df <- rbind(dist_emb, dist_ma, dist_js, dist_eu, dist_cophylog,  dist_cosine, dist_cvtree, dist_d2shepp, dist_d2star,  dist_fastme, dist_pl )



combined_df <- fread("pairwise_distances_kf2vec_cafe_R.txt", sep="\t", data.table = FALSE)

#combined_df=read.csv("pairwise_distances_kf2vec_cafe_R.txt",sep="\t", header=TRUE)

unique(combined_df$Condition)

#combined_df[combined_df$Condition=="cafe_phylip_js","Distance"]=1-combined_df[combined_df$Condition=="cafe_phylip_js","Distance"]
#combined_df[combined_df$Condition == "cafe_phylip_js", "Distance"] <- sqrt(combined_df[combined_df$Condition == "cafe_phylip_js", "Distance"])


nrow(combined_df)


combined_df2 = combined_df
combined_df2$Contig2 = combined_df$Contig1
combined_df2$Contig1 = combined_df$Contig2
combined_df2$Sample1 = combined_df$Sample2
combined_df2$Sample2 = combined_df$Sample1


  ######################################################################
  #Filtered for sister contigs
  
  # Find contigs that have at least one sibling (same genome, not itself)
  intra_contigs <- combined_df %>%
    filter(Sample1 == Sample2 & Contig1 != Contig2) %>%
    select(Contig1, Contig2) %>%
    unlist() %>%
    unique()

  length(intra_contigs)
    
  # Keep only rows where both contigs have siblings
  filtered_df <- rbind(combined_df,combined_df2) %>%
    filter(Contig1 %in% intra_contigs & Contig2 %in% intra_contigs)
  
  # Check how many rows were kept
  nrow(filtered_df)
  nrow(combined_df)
  

  
  # Table included in a paper (median) 
  df_summary_filt <- filtered_df %>%
    group_by(Sample1,Condition,Type) %>%
    dplyr::summarise(md=mean(Distance),md=mean(Distance),n=n())  %>% 
    pivot_wider(names_from = "Type", values_from = c("md","n")) %>%
    filter(!is.na(n_Within)) %>% 
    mutate(Fs = md_Across/(md_Within)) %>% 
    group_by(Condition) %>%
    filter(!is.na(Fs))%>%
    #dplyr::summarise(Fsmean=mean(Fs),n=n())
    dplyr::summarise(Fsmedian=median(Fs),n=n())
  ggplot(aes(x=Condition,y=Fs))+
    stat_summary()
  
  print(df_summary_filt)
  
  
sibling_filtered <- filtered_df %>%
    group_by(Contig1, Condition) %>%
    summarise(
      cl = Type[which.min(Distance)],  # assign closest type
      clw = sum(Type == "Within") > 0, # has at least one Within pair
      .groups = "drop"
    ) %>%
    filter(clw == TRUE) %>%
    group_by(Condition, cl) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = cl, values_from = n, values_fill = 0) %>%
    mutate(
      n = Across + Within,
      errorrate = Across / (Across + Within)
    )

  print(sibling_filtered,n=50)
  
  plot_df <- merge(sibling_filtered, df_summary_filt, by = "Condition") %>%
    mutate(pFi = 1 / Fsmedian) %>%
    select(Condition, errorrate, pFi) %>%
    filter(!Condition %in% c("cafe_phylip_co-phylog", "kf2vec_after_placement")) %>%
    pivot_longer(cols = c(errorrate, pFi), names_to = "name", values_to = "value")
  

  ggplot(plot_df, aes(x = name, y = value, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    theme_classic() +
    ylab("Value") +
    xlab(NULL) +
    scale_fill_manual(values = my_colors)  # if you have a color vector
  
  
  
  # plot
  plot_df %>%
    filter(name == "errorrate") %>%  # only plot errorrate
    ggplot(aes(x = reorder(Condition, -value), y = value, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    theme_classic() +
    xlab(NULL) +
    scale_y_continuous("Closest contig mismatch", labels = percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # optional: rotate labels
 
  
  
  # Annotated paper plot (final version)
  
  merge(sibling_filtered, df_summary_filt, by = "Condition") %>%
    #select(Condition, errorrate, pFi) %>%
    filter(!Condition %in% c("cafe_phylip_co-phylog", "kf2vec_after_placement")) %>%
    mutate(`F statistics` = Fsmedian, `Closest contig mismatch`=errorrate) %>%
    pivot_longer(cols = c(`Closest contig mismatch`, `F statistics`), names_to = "name", values_to = "value")%>% 
    filter(!Condition %in% c("cafe_phylip_js","cafe_phylip_cvtree")) %>%
    #filter(name == "errorrate") %>%  # only errorrate
  ggplot(aes(x = reorder(Condition, 1/value), y = value, fill = Condition)) +
    facet_wrap(.~name, scales="free_y",nrow=2,strip.position = "left")+
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_text(aes(label = scales::percent(value, accuracy = 1)), 
            vjust = +2, size = 2.5,color="yellow") +  # annotate bars with percentage
   scale_fill_manual(values = c( dark2_colors[7], dark2_colors[6], dark2_colors[6], dark2_colors[7], dark2_colors[7], dark2_colors[5], dark2_colors[5]))+ 
  scale_x_discrete( labels = c( 'Ma', 'Eu', expression(italic("D")[italic("2")]^italic("*")),  'Cosine', expression(italic("D")[italic("2")]^italic("S")), 'kf2vec', "kf2vec+tree")) +
    theme_classic() +
    xlab(NULL) +
    scale_y_continuous("",labels = percent)+
    theme(legend.position = "none",axis.title = element_blank())  # optional: rotate labels
  ggsave("contigs_across_within_bar_D3_sister_filtered.pdf",width=4.8,height = 4.0)
  
  getwd()
# sed '$!N;s/\n/,/' kmer_pl_3layer_lowerLrRate.out | awk 'BEGIN { FS = "[, ]" } ; {print $9, $21}' > kmer_pl_3layer_lowerLrRate.csv

  
  # OBSOLETE
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
  