setwd("~/Projects/TCR/results_1/proposal_figs/tcga_meta/")
library(data.table)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
#----------------------------
# figure 2, TCGA heatmap
#----------------------------
# read case and file meta
wgs_case <- read.table("WGS_case/chosen_case.txt",sep="\t",header=T)
wgs_file <- read.table("WGS_case/file_meta_case.tsv",sep="\t",header=T)

wxs_case <- read.table("WXS_case/chosen_case.txt",sep="\t",header=T)
wxs_file <- read.table("WXS_case/file_meta_case.tsv",sep="\t",header=T)

rna_case <- read.table("RNA_case/chosen_case.txt",sep="\t",header=T)
rna_file <- read.table("RNA_case/file_meta_case.tsv",sep="\t",header=T)

geno_case <- read.table("genoty_array_case/chosen_case.txt",sep="\t",header=T)
geno_file <- read.table("genoty_array_case/file_meta_case.tsv",sep="\t",header=T)

methy_case <- read.table("methy_array_case/chosen_case.txt",sep="\t",header=T)
methy_file <- read.table("methy_array_case/file_meta_case.tsv",sep="\t",header=T)

#---------------------------------------------
# extract common case have both sample sequenced
# get case_id union which have both blood and tumor seq in wgs or wxs

get_seq_case_id <- function(file_meta,case_file,overlap=2) {
  file_chosen <- file_meta[which(file_meta$Sample.ID %in% unique(case_file$sample_submitter_id)),c(1,2,7,8)]
  colnames(file_chosen)[3] <- "sample_submitter_id"

  merge_1 <- merge(case_file,file_chosen,by="sample_submitter_id",all=T)

  # remove NA
  merge_2<- na.omit(merge_1)
  file_df_merge <- unique(merge_2[,c(1,2,3,4,8)])

  # group case
  file_freq <- file_df_merge %>% count(project_id, case_submitter_id)

  common_case_name <- unique(file_freq$case_submitter_id[which(file_freq$n==overlap)])
  return(common_case_name)
}

#---------------------------------------------

wgs_common_case_name <- get_seq_case_id(wgs_file,wgs_case)
wxs_common_case_name <- get_seq_case_id(wxs_file,wxs_case)

wgs_wxs_case <- unique(c(wgs_common_case_name,wxs_common_case_name))

#write.table(wgs_wxs_case,"wgs_wxs_union_case_id.txt",sep="\t",quote = F,row.names = F,col.names = F)


#--------------------------------------------
# extract case ID from wgs_wxs_case for RNA,geno,methy
#--------------------------------------------
# RNA only have tumor
rna_case_id <- get_seq_case_id(rna_file,rna_case,overlap=1)

# genotype have both
geno_case_id <- get_seq_case_id(geno_file,geno_case,overlap=2)

# methy only have tumor
methy_case_id <- get_seq_case_id(methy_file,methy_case,overlap=1)

# filter by wgs_wxs_case
rna_case_chosen <- rna_case_id[which(rna_case_id %in% wgs_wxs_case)]
geno_case_chosen <- geno_case_id[which(geno_case_id %in% wgs_wxs_case)]
methy_case_chosen <- methy_case_id[which(methy_case_id %in% wgs_wxs_case)]

# extract project_id
wgs_chosen_project <- unique(wgs_case[which(wgs_case$case_submitter_id %in% wgs_common_case_name),c(1,3)])
wxs_chosen_project <- unique(wxs_case[which(wxs_case$case_submitter_id %in% wxs_common_case_name),c(1,3)])
rna_chosen_project <- unique(rna_case[which(rna_case$case_submitter_id %in% rna_case_chosen),c(1,3)])
geno_chosen_project <- unique(geno_case[which(geno_case$case_submitter_id %in% geno_case_chosen),c(1,3)])
methy_chosen_project <- unique(methy_case[which(methy_case$case_submitter_id %in% methy_case_chosen),c(1,3)])

wgs_wxs_chosen_project <- unique(rbind(wgs_chosen_project,wxs_chosen_project))

#count freq
wgs_feq <- data.frame(table(wgs_chosen_project$project_id))
colnames(wgs_feq)[2] <- "wgs_both"

wxs_feq <- data.frame(table(wxs_chosen_project$project_id))
colnames(wxs_feq)[2] <- "wxs_both"

rna_feq <- data.frame(table(rna_chosen_project$project_id))
colnames(rna_feq)[2] <- "rna_tumor"

geno_feq <- data.frame(table(geno_chosen_project$project_id))
colnames(geno_feq)[2] <- "geno_both"

methy_feq <- data.frame(table(methy_chosen_project$project_id))
colnames(methy_feq)[2] <- "methy_tumor"

wgs_wxs_feq <- data.frame(table(wgs_wxs_chosen_project$project_id))
colnames(wgs_wxs_feq)[2] <- "wgs_wxs_union"


# merge all
tmp_list <- list(wgs_feq,wxs_feq,wgs_wxs_feq,rna_feq,geno_feq,methy_feq)

tcga_count_for_heatmap <- tmp_list %>% reduce(full_join, by='Var1')

tcga_count_for_heatmap[is.na(tcga_count_for_heatmap)] <- 0

write.table(tcga_count_for_heatmap,"tcga_count_for_heatmap.txt",sep="\t",quote = F,row.names = F)

#-----------------
# plot heatmap
#-----------------
tcga_df_tmp <- read.table("tcga_count_for_heatmap.txt",sep="\t",header=T)
tcga2orien <- read.table("tcga_orien_cancer_name_for_heatmap.txt",sep="\t",header=T)
orien_df <- read.table("ORIEN_count.txt",sep="\t",header=T)

#------------------------
# add race and gander
#------------------------
clinical_meta <- as.data.frame(fread("../all_tcga_clinical_meta/clinical.tsv",select = c(2,3,12,15)))
wgs_wxs_case <- read.table("wgs_wxs_union_case_id.txt",sep="\t",header=F)
colnames(wgs_wxs_case) <- "case_submitter_id"
# subset clinical_meta by wgs_wxs_case_id

meta_chosen <- unique(merge(wgs_wxs_case,clinical_meta,by="case_submitter_id"))

meta_chosen_freq_gender <- data.frame(table(meta_chosen$project_id,meta_chosen$gender))
meta_chosen_freq_race <- data.frame(table(meta_chosen$project_id,meta_chosen$race))

meta_chosen_freq_race_gender <- data.frame(table(meta_chosen$project_id,
                                                 meta_chosen$race,
                                                 meta_chosen$gender))

meta_chosen_freq_race_gender$Var4 <- paste0(meta_chosen_freq_race_gender$Var2,"_",meta_chosen_freq_race_gender$Var3)
# long to width
gender_count <- spread(meta_chosen_freq_gender,Var2,Freq)
race_count <- spread(meta_chosen_freq_race,Var2,Freq)

gender_race_count <- spread(meta_chosen_freq_race_gender[,-c(2,3)],Var4,Freq)

meta_count_merge <-merge(gender_count,race_count,by="Var1",all=T)

meta_count_merge_2 <- merge(meta_count_merge,gender_race_count,by="Var1")
# add tcga_df count
tcga_df <- merge(tcga_df_tmp,meta_count_merge_2,by="Var1",all=T)
# replace tcga_name
tcga_2 <- merge(tcga_df,tcga2orien[,c(1,3)],by="Var1")

write.table(tcga_df,"tcga_count_for_heatmap_clinical_meta_gender_per_race.txt",sep="\t",quote = F,row.names = F)


# sum same name
tcga_2 <- tcga_2[,-1]
tcga_2$heatmap_name <- as.factor(tcga_2$heatmap_name)
tcga_3 <- aggregate(.~heatmap_name,tcga_2,sum)

# merge orien data
colnames(orien_df) <-c("heatmap_name","WXS_ORIEN","RNA_ORIEN")

heatmap_df <- merge(tcga_3,orien_df,by="heatmap_name",all=T)
heatmap_df[is.na(heatmap_df)] <- 0


row.names(heatmap_df) <- heatmap_df$heatmap_name
heatmap_df_1 <- heatmap_df[,-1]


# heatmap 1
# heatmap_df_1 <- log2(heatmap_df_1+1)
# heatmap_df_1 <- heatmap_df[,-1]

# change column order
heatmap_df_1 <- heatmap_df_1[,c(1:6,15:16,7:14)]
sample_col <- data.frame(Database=rep(c("TCGA","ORIEN","Gender","Race"),c(6,2,2,6)))
row.names(sample_col) <- colnames(heatmap_df_1)

# database color
my_colour <- list(Database = c(TCGA = brewer.pal(3, "Blues")[1], ORIEN = brewer.pal(3, "Blues")[3],
                               Gender=brewer.pal(5, "Purples")[2],
                               Race=brewer.pal(5, "Purples")[4]))
#heatmap_df_1$cancer <- rownames(heatmap_df_1)
#write.table(heatmap_df_1,"../fig2_heatmap_df.txt",quote = F,sep="\t",row.names = F)

#write.table(heatmap_df_1,"../fig2_heatmap_df_w_clinical.txt",quote = F,sep="\t",row.names = F)

pheatmap(heatmap_df_1,
         color = c("white",brewer.pal(6, "Reds")),
         breaks=c(0,1,150,500,1000,1500,2000,2500),
         gaps_col = 8,
         cluster_rows = F,
         cluster_cols = F,
         annotation_colors = my_colour,
         annotation_col = sample_col,
         )

# heatmap 2
heatmap_df_2 <- data.frame(t(heatmap_df_1))

sample_row <- data.frame(Database=rep(c("TCGA","ORIEN"),c(6,2)))
row.names(sample_row) <- row.names(heatmap_df_2)

# database color
my_colour <- list(Database = c(TCGA = brewer.pal(3, "Dark2")[1], ORIEN = brewer.pal(3, "Dark2")[2]))

pheatmap(heatmap_df_2,
         color = brewer.pal(8, "Blues"),
         #breaks=c(0,3,5,7,9,10),
         gaps_col = 6,
         cluster_rows = F,
         cluster_cols = F,
         annotation_colors = my_colour,
         annotation_row = sample_row,
)

#-----------------
# figure 3, cell boxplot
#-----------------
mm_estimate <- read.table("../../ovarian/o_estimation_matrix.csv",sep=",",header=T)
mm_estimate <- as.data.frame(t(mm_estimate))

mm_estimate$SL_tumor <- rownames(mm_estimate)
mm_estimate[1,] <- gsub(" ","_",mm_estimate[1,])
colnames(mm_estimate) <- mm_estimate[1,]
mm_estimate <- mm_estimate[-1,]
colnames(mm_estimate)[120] <- "SL_tumor"


mm_meta <- read.table("../../../metadata/mm_tumor_meta.txt",sep="\t",header=T)
colnames(mm_meta)[2] <- "CH"
mm_estimate_merge_tmp <- merge(mm_meta[,c("SL_tumor","CH")],mm_estimate,by="SL_tumor")
# remove all 0
mm_estimate_merge_tmp[,-c(1,2)] <- mm_estimate_merge_tmp[,-c(1,2)] %>% mutate_if(is.character,as.numeric)
mm_estimate_merge <- mm_estimate_merge_tmp[, colSums(mm_estimate_merge_tmp != 0) > 0]

mm_estimate_long <- gather(mm_estimate_merge,cell_p,value,
                           B_cell_TIMER:Cancer_associated_fibroblast_MCPCOUNTER,factor_key = T)


mm_estimate_long$value <- as.numeric(mm_estimate_long$value)
colnames(mm_estimate_long)[2] <-"CH"
mm_estimate_long$CH <- as.factor(mm_estimate_long$CH)

mm_estimate_long_eso <- mm_estimate_long[which(mm_estimate_long$cell_p=="Eosinophil_CIBERSORT"),]
mm_estimate_long_eso$cancer <- "Ovarian"

#---------------------------------
mm_eso <- mm_estimate_long_eso
o_eso <- mm_estimate_long_eso
o_eso$CH <- ifelse(o_eso$CH==1,"CH","NO")

box_plot_df <- rbind(mm_eso,o_eso)

write.table(box_plot_df,"../fig3_boxplot_df.txt",sep="\t",row.names = F,quote = F)
sig_v <- box_plot_df %>%
  group_by(cancer) %>%
  wilcox_test(value ~ CH) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")
sig_v$p <- round(sig_v$p,digits=4)

sig_v <-  sig_v %>% add_xy_position(x = "CH")


ggplot(box_plot_df, aes(x=cancer, y=value,fill=CH)) +
  geom_boxplot(width=0.5,position = position_dodge(width=0.6))+
  #facet_wrap(~cancer)+
  theme_classic()+
  scale_y_continuous(trans='sqrt')+
  scale_fill_manual(values=c("firebrick","skyblue"))+
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_blank(),
    axis.text=element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank())

#---------------------------------
# figure 4 bar plot of sig vdjc count
#---------------------------------
# fisher test
b_fisher <- read.table("../../breast/b_vdjc_fisher.txt",sep="\t",header=T)
l_fisher <- read.table("../../lung/l_vdjc_fisher.txt",sep="\t",header=T)
o_fisher <- read.table("../../ovarian/o_vdjc_fisher.txt",sep="\t",header=T)
mm_fisher <- read.table("../../mm/mm_vdjc_fisher.txt",sep="\t",header=T)

# significant vdjc
b_sig <- b_fisher$VDJC[which(b_fisher$pval < 0.05)]
l_sig <- l_fisher$VDJC[which(l_fisher$pval < 0.05)]
o_sig <- o_fisher$VDJC[which(o_fisher$pval < 0.05)]
mm_sig <- mm_fisher$VDJC[which(mm_fisher$pval < 0.05)]

bar_df <- data.frame(cancer=c("mm","b","o","l"),
                     count=c(length(unique(mm_sig)),
                          length(unique(b_sig)),
                          length(unique(o_sig)),
                          length(unique(l_sig))))

bar_df$cancer <- factor(bar_df$cancer,levels=c("mm","b","o","l"))
bar_df$normal <- c(bar_df$count[1]/648,bar_df$count[2]/291,bar_df$count[3]/127,bar_df$count[4]/69)
ggplot(bar_df,aes(x=cancer,y=normal,fill=cancer))+
  geom_bar(stat="identity",width = 0.8)+
  theme_classic()+
  scale_fill_manual(values=brewer.pal(4, "Spectral"))+
  theme(plot.title = element_blank(),
        axis.text=element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())





