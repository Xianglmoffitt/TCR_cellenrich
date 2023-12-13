setwd("~/Projects/TCR/")

all_report <- read.table("~/Projects/TCR/blo_trust4_out_all_report.tsv",sep="\t",header=T)
sample_meta <- read.table("~/Projects/TCR/metadata/blo_tumor_meta.txt",sep="\t",header=T)

colnames(sample_meta)[2] <- "CH"

sample_meta$CH <- ifelse(sample_meta$CH=="1","CH","NO")

# separete cancer (BRE - Breast Cancer, GYN - Ovarian Cancer, THO - Lung Cancer,)

sample_meta <- sample_meta[which(sample_meta$Cancer_type=="THO - Lung Cancer"),]

#tmp_df <- unique(all_report[which(all_report$VDJC %in% target_list$V1),c("VDJC","sample_id")])
tmp_df <- unique(all_report[,c("VDJC","sample_id")])


colnames(tmp_df)[2]<-"SL_tumor"


tmp_2 <- merge(tmp_df,sample_meta[,c(2,6)],by="SL_tumor")
tmp_3 <- tmp_2 %>% count(VDJC, CH, sort = TRUE)
tmp_wide <- spread(tmp_3, CH, n)

tmp_wide[is.na(tmp_wide)] <- 0

tmp_wide$group <- paste0(tmp_wide$CH,"_",tmp_wide$NO)

group_list <- unique(tmp_wide$group)
PT <- nrow(sample_meta)
PH <- length(which(sample_meta$CH=="CH"))

fisher_out <- data.frame()
for (i in 1:length(group_list)) {
  g_tmp <- as.integer(unlist(strsplit(group_list[i],"_")))


  LT <- g_tmp[1]+g_tmp[2]
  LH <- g_tmp[1]

  contigency_m <- matrix(c(LH,LT-LH,PH-LH,PT-LT-(PH-LH)),nrow=2,
                         dimnames= list(c("CH", "Not CH"),
                                        c("target_sample", "rest_sample")))
  ft <- fisher.test(contigency_m)
  fisher_out_tmp <- data.frame(group=group_list[i],
                               pval=ft$p.value)
  fisher_out <- rbind(fisher_out,fisher_out_tmp)
}

# merge with nt name
fisher_out_2 <- merge(tmp_wide,fisher_out,by="group")

write.table(fisher_out_2, file="results_1/l_vdjc_fisher.txt", sep="\t",row.names=FALSE,quote=F)


# p.adj
p_in <- read.table("results_1/mm_vdjc_fisher.txt",sep="\t",header=T)

group_out <- unique(p_in[,c("group","pval")])

group_out$padj <- p.adjust(group_out$pval,method = "BH")

write.table(group_out, file="results_1/mm_vdjc_fisher_group_padj.txt", sep="\t",row.names=FALSE,quote=F)


