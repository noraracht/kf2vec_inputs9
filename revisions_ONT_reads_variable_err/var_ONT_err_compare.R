

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

library(dplyr)
library("tidyverse")

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


df =read.csv('ONT_best10TOLquer_varReadErr.csv',sep=",",header=TRUE)
nrow(df)
head(df)
# Replace 1 with 100 where read_err == 1
df$read_err[df$read_err == 1] <- 100

colnames(df)
df
#e1 = df_fin[df_fin$greedy == "ng" & (df_fin$clsf!="tr" |df_fin$chunked=="Chunked") & !grepl("Exp",df_fin$clsf),]

#e2 = e1[e1$chunks %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "1"),]
#unique(e2$cond)
#e3= e2[e2$cond=="Claded - Chunked - cmp",]
#e4= e2[e2$cond=="Claded - Chunked - tr",]
#summary(e3$V3)
#summary(e4$V3)

#e5 = e1[!e1$chunks %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "1", "full"),]
#head(e5)
#e6= e5[e5$cond=="Claded - Chunked - cmp",]
#summary(e6$V3)


ggplot(aes(x=read_err, y=pl_err, color = factor(interaction(model, clade_ind))),
       data=df)+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  stat_summary(geom="line")+
  #geom_violin()+
  theme_classic()+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  #coord_cartesian(xlim=c(0.85,1))
  
  #scale_x_discrete(name="Fragment length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  scale_colour_brewer(palette = "Dark2", name="")
  scale_linetype(name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  #facet_wrap(~claded)+
  #guides(col= guide_legend(title= "Conditions"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(legend.position = c(0.75,0.7), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
#scale_colour_brewer(palette = "Dark2", name="")

#ggsave("ONT-line.pdf",width=6.5,height = 4.5)

head(df[, c("clade_ind", "model")])
df$facet_label <- paste(df$clade_ind, df$model, sep = "_")
df$facet_label

# Explicit order + relabel
df$facet_label <- factor(df$facet_label,
                         levels = c("classified_Uncladed_Chunked", "classified_Claded_Chunked", "true_Claded_Chunked"),
                         labels = c("Uncladded",
                                    "Cladded",
                                    "Cladded (true)"))

unique(df$facet_label)

ggplot(aes(x=cut(length_adj/10^3,c(5,6,8,10,13,18,70),include.lowest = TRUE), y=pl_err, 
           color = factor(interaction(read_err))),
       data=df)+
  #geom_boxplot(alpha = 0.6)+
  #geom_point(aes(color = factor(k))) +
  #stat_summary(geom="crossbar")
  #stat_summary(geom="line")+
  #geom_violin()+
  theme_classic()+
  #facet_wrap(.~interaction(clade_ind,model,lex.order = F))+
  facet_wrap(~facet_label)+
  stat_summary(aes(group=interaction(model,clade_ind,read_err)), alpha = 0.7,geom="line")+
  stat_summary( alpha = 0.7,geom="point")+
  stat_summary( alpha = 0.7,geom="errorbar",width=0.1,linetype=1,size=0.3)+
  #coord_cartesian(xlim=c(0.85,1))
  scale_x_discrete(name="Read length (KB)")+ 
  ylab("Placement error")+
  #scale_x_discrete(name="Fragment length (KB)", label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10","20", "30", "40", "50", "100", "200", "300", "400", "500", "1000", "full"))  +
  ylab("Placement error")+
  scale_colour_brewer(palette = "Dark2", name="Read accuracy", label = c("0.85", "0.90", "0.95", "0.99", "1.00"))+
#scale_linetype(name="", labels = c("Cladded", "Cladded (true)", "Uncladded"))+
  #facet_wrap(~claded)+
  #guides(col= guide_legend(title= "Conditions"))+
  #scale_fill_discrete(name = "New Legend Title")
  #guides(fill=guide_legend(title="k"))
  theme(axis.text.x = element_text( angle = 30, hjust = 1, vjust = 1))+
  theme(legend.position = c(0.90,0.7), legend.margin=margin(0,0,0,0),
        #axis.text.x = element_text(size = 8)
  )
#scale_colour_brewer(palette = "Dark2", name="")

ggsave("variable-ONTerr-D1-line.pdf",width = 6.5,height = 3.5)
getwd()





# Compute fraction of equal rows per condition
result_df <- df %>%
  group_by(read_err, model) %>%
  summarize(
    fraction_equal = mean(top_class == clade)
  )

print(result_df)


df$cond_clfacc <- paste(as.character(df$read_err), df$model, df$clade_ind, sep = "_")


unique(df$cond_clfacc)

result_df <- df %>%
  group_by(interaction(read_err, model, clade_ind)) %>%
  summarize(
    fraction_equal = mean(top_class == clade)
  )
result_df


df_subset <- df[df$cond_clfacc == "85_Claded_Chunked_classified", ]
df_correct_class <- df_subset[df_subset$top_class == df_subset$clade, ]
round(nrow(df_correct_class)/nrow(df_subset), 2)

df_subset <- df[df$cond_clfacc == "90_Claded_Chunked_classified", ]
df_correct_class <- df_subset[df_subset$top_class == df_subset$clade, ]
round(nrow(df_correct_class)/nrow(df_subset), 2)

df_subset <- df[df$cond_clfacc == "95_Claded_Chunked_classified", ]
df_correct_class <- df_subset[df_subset$top_class == df_subset$clade, ]
round(nrow(df_correct_class)/nrow(df_subset), 2)

df_subset <- df[df$cond_clfacc == "99_Claded_Chunked_classified", ]
df_correct_class <- df_subset[df_subset$top_class == df_subset$clade, ]
round(nrow(df_correct_class)/nrow(df_subset), 2)

df_subset <- df[df$cond_clfacc == "100_Claded_Chunked_classified", ]
df_correct_class <- df_subset[df_subset$top_class == df_subset$clade, ]
round(nrow(df_correct_class)/nrow(df_subset), 2)

nrow(df_correct_class)
nrow(summary_df)

head(df)


# sed '$!N;s/\n/,/' kmer_pl_3layer_lowerLrRate.out | awk 'BEGIN { FS = "[, ]" } ; {print $9, $21}' > kmer_pl_3layer_lowerLrRate.csv
