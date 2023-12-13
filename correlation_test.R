setwd("~/Projects/TCR/cell_enrichment/")


mm_estimate <- read.table("blo_timer_estimation_matrix.csv",sep=",",header=T)

mm_estimate <- as.data.frame(t(mm_estimate))

mm_estimate$SL_tumor <- rownames(mm_estimate)
mm_estimate[1,] <- gsub(" ","_",mm_estimate[1,])
colnames(mm_estimate) <- mm_estimate[1,]
mm_estimate <- mm_estimate[-1,]
colnames(mm_estimate)[120] <- "SL_tumor"


mm_meta_in <- read.table("../metadata/blo_tumor_meta.txt",sep="\t",header=T)

# separete cancer (BRE - Breast Cancer, GYN - Ovarian Cancer, THO - Lung Cancer,)
mm_meta <- mm_meta_in[which(mm_meta_in$Cancer_type=="GYN - Ovarian Cancer"),]


mm_estimate_merge_tmp <- merge(mm_meta[,c("SL_tumor","CH_status")],mm_estimate,by="SL_tumor")
# remove all 0
mm_estimate_merge_tmp[,-c(1,2)] <- mm_estimate_merge_tmp[,-c(1,2)] %>% mutate_if(is.character,as.numeric)
mm_estimate_merge <- mm_estimate_merge_tmp[, colSums(mm_estimate_merge_tmp != 0) > 0]

mm_estimate_long <- gather(mm_estimate_merge,cell_p,value,
                           B_cell_TIMER:Cancer_associated_fibroblast_MCPCOUNTER,factor_key = T)


mm_estimate_long$value <- as.numeric(mm_estimate_long$value)
colnames(mm_estimate_long)[2] <-"CH"
mm_estimate_long$CH <- as.factor(mm_estimate_long$CH)



# wilcox test
stat.test <- mm_estimate_long %>%
  group_by(cell_p) %>%
  wilcox_test(value ~ CH) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

wilcox_out_final <- stat.test

write.table(wilcox_out_final,"../results_1/o_wilcox_out.txt",sep="\t",row.names = F,quote = F)




#-------------------------------------------------------------
# correlation and PCA
#-------------------------------------------------------------
# code CH
mm_estimate_merge$CH <- ifelse(mm_estimate_merge$CH=="CH",1,0)

mm_estimate_merge[,-1]<- lapply(mm_estimate_merge[,-1],as.numeric)

correlation_m <- cor(mm_estimate_merge[,-c(1)])

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

#boxplot
library(ggpubr)
library(rstatix)
library(tidyverse)
# target_v <- wilcox_out_final[which(wilcox_out_final$pval<=0.05),]
#
# plot_m <- mm_estimate_merge[,c(2,which(colnames(mm_estimate_merge)%in%target_v$variable))]


mm_estimate_long <- gather(mm_estimate_merge,cell_p,value,
                           B_cell_TIMER:Cancer_associated_fibroblast_MCPCOUNTER,factor_key = T)


mm_estimate_long$value <- as.numeric(mm_estimate_long$value)
colnames(mm_estimate_long)[2] <-"CH"
mm_estimate_long$CH <- as.factor(mm_estimate_long$CH)

# wilcox test
stat.test <- mm_estimate_long %>%
  group_by(cell_p) %>%
  wilcox_test(value ~ CH) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

wilcox_out_final <- cbind(stat.test,CH_cor[-1,])

write.table(wilcox_out_final,"blo_CH_wilcox_test.txt",sep="\t",row.names = F,quote = F)

# sign_v
sig_v <- stat.test[which(stat.test$p<=0.05),]
plot_df <- mm_estimate_long[which(mm_estimate_long$cell_p %in% sig_v$cell_p),]

sig_v <- plot_df %>%
  group_by(cell_p) %>%
  wilcox_test(value ~ CH) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")
sig_v$p <- round(sig_v$p,digits=4)

sig_v <-  sig_v %>% add_xy_position(x = "CH")


gtmp <- ggplot(plot_df,aes(x=CH,y=value,color=CH))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap(~cell_p,nrow=3)+
  stat_pvalue_manual(sig_v,hide.ns = T,bracket.nudge.y = 0,label="{p}{p.signif}")+
  scale_y_continuous(trans='exp')

pdf("mm_significant_variables.pdf",width=16,height = 9)
print(gtmp)
dev.off()






# PCA
library(factoextra)
library(ggbiplot)


res.pca <- prcomp(mm_estimate_merge[,-c(1,2)],scale. = TRUE,center=T)

g <- ggbiplot(res.pca,
              groups = as.factor(mm_estimate_merge$CH),
              ellipse = TRUE,
              circle = F)+
  theme_classic()
