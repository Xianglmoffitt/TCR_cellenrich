setwd("~/Projects/TCR/")

# fisher test
b_fisher <- read.table("results_1/breast/b_cdr3aa_fisher.txt",sep="\t",header=T)
l_fisher <- read.table("results_1/lung/l_cdr3aa_fisher.txt",sep="\t",header=T)
o_fisher <- read.table("results_1/ovarian/o_cdr3aa_fisher.txt",sep="\t",header=T)
mm_fisher <- read.table("results_1/mm/mm_cdr3aa_fisher.txt",sep="\t",header=T)

# significant vdjc
b_sig <- b_fisher$CDR3aa[which(b_fisher$pval < 0.05)]
l_sig <- l_fisher$CDR3aa[which(l_fisher$pval < 0.05)]
o_sig <- o_fisher$CDR3aa[which(o_fisher$pval < 0.05)]
mm_sig <- mm_fisher$CDR3aa[which(mm_fisher$pval < 0.05)]

# significant vdjc df

vdjc_df <- data.frame(cdr3aa=unique(c(b_sig,l_sig,o_sig,mm_sig)),
                      b=0,
                      l=0,
                      o=0,
                      mm=0)

vdjc_df$b[which(vdjc_df$cdr3aa %in% b_sig)] <- 1
vdjc_df$l[which(vdjc_df$cdr3aa %in% l_sig)] <- 1
vdjc_df$o[which(vdjc_df$cdr3aa %in% o_sig)] <- 1
vdjc_df$mm[which(vdjc_df$cdr3aa %in% mm_sig)] <- 1

vdjc_df$sum <- rowSums(vdjc_df[,c(2:5)])

vdjc_df_sort <- vdjc_df[order(-vdjc_df$sum),]

write.table(vdjc_df_sort,"results_1/sig_cdr3aa_cancer.txt",sep="\t",row.names = F,quote = F)



#---------------------------------------
# wilcox cell plot
#---------------------------------------
mm_estimate <- read.table("results_1/mm/mm_estimation_matrix.csv",sep=",",header=T)

mm_estimate <- as.data.frame(t(mm_estimate))

mm_estimate$SL_tumor <- rownames(mm_estimate)
mm_estimate[1,] <- gsub(" ","_",mm_estimate[1,])
colnames(mm_estimate) <- mm_estimate[1,]
mm_estimate <- mm_estimate[-1,]
colnames(mm_estimate)[120] <- "SL_tumor"


mm_meta <- read.table("metadata/mm_tumor_meta.txt",sep="\t",header=T)

mm_estimate_merge_tmp <- merge(mm_meta[,c("SL_tumor","CH")],mm_estimate,by="SL_tumor")
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


#write.table(stat.test,"results_1/ovarian/o_wilcox_out.txt",sep="\t",row.names = F,quote = F)

# boxplot of wilcox test
#boxplot
library(ggpubr)
library(rstatix)
library(tidyverse)
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
  stat_pvalue_manual(sig_v,hide.ns = T,bracket.nudge.y = 0.3,label="{p}{p.signif}")+
  scale_y_continuous(trans='sqrt')

pdf("results_1/mm/mm_sig_cell_boxplot.pdf",width=16,height = 9)
print(gtmp)
dev.off()





