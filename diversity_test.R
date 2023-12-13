setwd("~/Projects/TCR/")

library(vegan)
sample_list <- read.table("metadata/blo_sample_list",sep="\t",header=F)


# combine file
all_report <- data.frame()
for (i in 1:length(sample_list$V1)) {
  print(i)
  in_name <- paste0("blo_trust4/",sample_list$V1[i],"_report.tsv")
  in_tmp <- read.table(in_name,sep="\t",header=T,comment.char="")

  in_tmp$sample_id <- sample_list$V1[i]
  in_tmp$VDJC <- paste(in_tmp$V,in_tmp$D,in_tmp$J,in_tmp$C,sep="_")

  all_report <- rbind(all_report,in_tmp)

}

write.table(all_report,"blo_trust4_out_all_report.tsv",sep="\t",row.names = F,quote=F)

# make to OncoDiversity input file

onco_in_cdr3nt <- data.frame(patient=all_report$sample_id,
                             CDR3=all_report$CDR3nt,
                             count=all_report$X.count)

cdr3aa_tmp <-aggregate(X.count~CDR3aa+sample_id,data=all_report[,c(1,4,11)],FUN="sum")
VDJC_tmp <-aggregate(X.count~VDJC+sample_id,data=all_report[,c(1,11,12)],FUN="sum")

onco_in_cdr3aa <- data.frame(patient=cdr3aa_tmp$sample_id,
                             CDR3=cdr3aa_tmp$CDR3aa,
                             count=cdr3aa_tmp$X.count)
onco_in_VDJC <- data.frame(patient=VDJC_tmp$sample_id,
                             CDR3=VDJC_tmp$VDJC,
                             count=VDJC_tmp$X.count)

test_sample <- head(sample_list$V1)
write.table(onco_in_cdr3nt[which(onco_in_cdr3nt$patient %in% test_sample[1]),],"onco_cdr3nt_input_test_1.csv",sep=",",row.names = F,quote=F)
write.table(onco_in_cdr3nt[which(onco_in_cdr3nt$patient %in% test_sample[c(1,2)]),],"onco_cdr3nt_input_test_2.csv",sep=",",row.names = F,quote=F)

write.table(onco_in_cdr3nt,"blo_onco_cdr3nt_input.csv",sep=",",row.names = F,quote=F)
write.table(onco_in_cdr3aa,"blo_onco_cdr3aa_input.csv",sep=",",row.names = F,quote=F)
write.table(onco_in_VDJC,"blo_onco_VDJC_input.csv",sep=",",row.names = F,quote=F)


#----------------------------------------------------------
# GDI
#----------------------------------------------------------
# gdi_in_df <- all_report[,c(1,3,11,12)]
#
# cdr3nt_out <- data.frame()
#
# for (s in 1:nrow(sample_list)) {
#   print(s)
#   cdr3_tmp <- gdi_in_df[which(gdi_in_df$sample_id==sample_list$V1[s]),]
#   cdr3_tmp$frequency <- cdr3_tmp$X.count/sum(cdr3_tmp$X.count)
#
#   vdjc_tmp <-aggregate(X.count~VDJC,data=cdr3_tmp[,c(1,4)],FUN="sum")
#   vdjc_tmp$frequency <- vdjc_tmp$X.count/sum(vdjc_tmp$X.count)
#
#   # shannon index
#   shan_cdr3 <- diversity(cdr3_tmp$frequency,index = "shannon")
#   shan_vdjc <- diversity(vdjc_tmp$frequency,index = "shannon")
#
#   # simpson index
#   simp_cdr3 <- diversity(cdr3_tmp$frequency,index = "simpson")
#   simp_vdjc <- diversity(vdjc_tmp$frequency,index = "simpson")
#
#   # inverse simp
#   invsimp_cdr3 <- diversity(cdr3_tmp$frequency,index = "invsimpson")
#   invsimp_vdjc <- diversity(vdjc_tmp$frequency,index = "invsimpson")
#
#   # 0.01D (low q)
#   gdi_001_cdr3  <- sum(cdr3_tmp$frequency^0.01)^(1/(1-0.01))
#   gdi_001_vdjc  <- sum(vdjc_tmp$frequency^0.01)^(1/(1-0.01))
#
#   # 100D (high q)
#   gdi_100_cdr3  <- sum(cdr3_tmp$frequency^100)^(1/(1-100))
#   gdi_100_vdjc  <- sum(vdjc_tmp$frequency^100)^(1/(1-100))
#
#   # diff_D
#   diff_d_cdr3 <- gdi_001_cdr3-gdi_100_cdr3
#   diff_d_vdjc <- gdi_001_vdjc-gdi_100_vdjc
#
#   # GDI plot q from 0.01~100
#   q <- seq(0.01,100,0.1)
#   gdi_df <- data.frame()
#   for (i in 1:length(q)) {
#     print(i)
#     gdi_cdr3_tmp <- sum(cdr3_tmp$frequency^q[i])^(1/(1-q[i]))
#     gdi_vdjc_tmp <- sum(vdjc_tmp$frequency^q[i])^(1/(1-q[i]))
#
#     gdi_df_tmp <- data.frame(q=q[i],
#                              gdi_cdr3=gdi_cdr3_tmp,
#                              gdi_vdjc=gdi_vdjc_tmp)
#     gdi_df <- rbind(gdi_df,gdi_df_tmp)
#   }
#
#
#
# }
#
# plot(log(gdi_df$q),gdi_df$gdi_vdjc,pch=16,cex=1)
#
# ggplot(gdi_df,aes(x=log(q),y=gdi_cdr3))+
#   geom_smooth(method = "loess",span=0.1)
#
# ggplot(gdi_df,aes(x=log(q),y=gdi_vdjc))+
#   geom_smooth(method = "loess",span=0.1)
#
# # IP_p


#----------------------------------------------------------
# correlation test
#----------------------------------------------------------
all_report <- read.table("blo_trust4_out_all_report.tsv",sep="\t",header=T)
gdi_cdr3nt <- read.table("diversity_results/blo_onco_cdr3nt.csv",sep=",",header=T)
gdi_cdr3aa <- read.table("diversity_results/blo_onco_cdr3aa.csv",sep=",",header=T)
gdi_vdjc <- read.table("diversity_results/blo_onco_vdjc.csv",sep=",",header=T)

sample_meta <- read.table("metadata/blo_tumor_meta.txt",sep="\t",header=T)
#cdr3nt_count <- read.table("cdr3nt_count",sep="\t",header=F)

freq_tmp <- as.data.frame(table(all_report$CDR3nt))
summary(freq_tmp$Freq)

# rename columns
colnames(cdr3nt_count) <- c("cdr3nt_count","SL_tumor")
colnames(gdi_cdr3nt) <- c("SL_tumor","lowQ_cdr3nt","highQ_cdr3nt",
                          "deltaqD_cdr3nt","IPq_cdr3nt","IPslope_cdr3nt",
                          "qDatIPq_cdr3nt","Shan_cdr3nt","Simp_cdr3nt")

colnames(gdi_cdr3aa) <- c("SL_tumor","lowQ_cdr3aa","highQ_cdr3aa",
                          "deltaqD_cdr3aa","IPq_cdr3aa","IPslope_cdr3aa",
                          "qDatIPq_cdr3aa","Shan_cdr3aa","Simp_cdr3aa")

colnames(gdi_vdjc) <- c("SL_tumor","lowQ_vdjc","highQ_vdjc",
                          "deltaqD_vdjc","IPq_vdjc","IPslope_vdjc",
                          "qDatIPq_vdjc","Shan_vdjc","Simp_vdjc")

# merge to one

mtmp1 <- merge(sample_meta,gdi_cdr3nt,by="SL_tumor")
mtmp2 <- merge(mtmp1,gdi_cdr3aa,by="SL_tumor")
mtmp3 <- merge(mtmp2,gdi_vdjc,by="SL_tumor")

all_var <- merge(mtmp3,cdr3nt_count,by="SL_tumor")

write.table(mtmp3,"results_1/blo_all_var.txt",sep="\t",row.names = F,quote = F)

#-------------------------------------------------------------
# fisher test
#-------------------------------------------------------------
# background (CH:48,non-CH:439)
#table(sample_meta$CH_status)

cdr3nt_list <- unique(all_report[,3,drop=F])

PT <- nrow(sample_meta)
PH <- 48


cdr3nt_list$pval<- apply(cdr3nt_list, MARGIN = 1, fisher_function_by_row)


fisher_test_by_row <- function(row) {
  tmp_df <- unique(all_report[which(all_report$CDR3nt==row[1]),c("CDR3nt","sample_id")])
  colnames(tmp_df)[2]<-"SL_tumor"
  # add meta
  tmp_2 <- merge(tmp_df,sample_meta[,c(2,6)],by="SL_tumor")
  LT <- nrow(tmp_2)
  LH <- length(which(tmp_2$CH_status==1))

  contigency_m <- matrix(c(LH,LT-LH,PH-LH,PT-LT-(PH-LH)),nrow=2,
                         dimnames= list(c("CH", "Not CH"),
                                        c("target_sample", "rest_sample")))
  ft <- fisher.test(contigency_m)
  return(ft$p.value)
}


#-------------------------------------------------------------
# correlation and PCA
#-------------------------------------------------------------
# code CH
all_var$CH <- ifelse(all_var$CH=="CH",1,0)

correlation_m <- cor(all_var[,-c(1:3,5)])

# correlation heatmap
# lower triangle
correlation_m[upper.tri(correlation_m)] <-NA

lower_tri <- correlation_m
library(reshape2)
library(ggplot2)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
cormat <- reorder_cormat(correlation_m)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()


CH_cor <- data.frame(correlation_m[,1])

# PCA
library(factoextra)
library(ggbiplot)


res.pca <- prcomp(all_var[,-c(1:5)],scale. = TRUE,center=T)

g <- ggbiplot(res.pca,
              groups = as.factor(all_var$CH),
              ellipse = TRUE,
              circle = F)+
  theme_classic()







