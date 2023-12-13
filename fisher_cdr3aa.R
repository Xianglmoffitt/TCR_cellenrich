setwd("~/Projects/TCR/")

all_report <- read.table("trust4_out_all_report.tsv",sep="\t",header=T)
sample_meta <- read.table("metadata/mm_tumor_meta.txt",sep="\t",header=T)
CDR3aa_list <- unique(all_report[,c("CDR3nt","sample_id")])

colnames(sample_meta)[2] <- "CH"
CH_sample <- sample_meta$SL_tumor[which(sample_meta$CH!="CH")]
target_freq <- data.frame(table(CDR3aa_list$CDR3nt))

# filter target
CDR3aa_1 <- target_freq$Var1[which(target_freq$Freq==1)]

CDR3aa_list_filtered_1 <- CDR3aa_list[which(CDR3aa_list$CDR3nt %in% CDR3aa_1),]

CDR3aa_list_filtered_2 <- CDR3aa_list_filtered_1[which(CDR3aa_list_filtered_1$sample_id %in% CH_sample),]

remove_list <- unique(CDR3aa_list_filtered_2$CDR3nt)

# final list
target_list <- unique(CDR3aa_list$CDR3nt[!(CDR3aa_list$CDR3nt %in%remove_list)])
write.table(target_list,"fisher/blo_cdr3nt_target_list.txt",sep="\t",quote=F,row.names = F,col.names = F)

target_df <- read.table("fisher/blo_cdr3aa_target_list.txt",sep="\t",header=F)

target_list <- target_df$V1
PT <- nrow(sample_meta)
PH <- length(which(sample_meta$CH=="CH"))

fisher_out <- data.frame()
for (i in 1:length(target_list)) {
  if(i%%5000==0) {
    print(i)
  }


  tmp_df <- unique(all_report[which(all_report$CDR3aa==target_list[i]),c("CDR3aa","sample_id")])
  colnames(tmp_df)[2]<-"SL_tumor"
  # add meta
  tmp_2 <- merge(tmp_df,sample_meta[,c(2,6)],by="SL_tumor")
  LT <- nrow(tmp_2)
  LH <- length(which(tmp_2$CH==1))
  #LH <- length(which(tmp_2$CH=="CH"))

  contigency_m <- matrix(c(LH,LT-LH,PH-LH,PT-LT-(PH-LH)),nrow=2,
                         dimnames= list(c("CH", "Not CH"),
                                        c("target_sample", "rest_sample")))
  ft <- fisher.test(contigency_m)
  fisher_out_tmp <- data.frame(target=target_list[i],
                               pval=ft$p.value)
  fisher_out <- rbind(fisher_out,fisher_out_tmp)
}

write.table(fisher_out,"fisher/fisher_out/blo_cdr3nt_fisher.txt",sep="\t",quote=F,row.names = F)




tmp_df <- unique(all_report[which(all_report$VDJC=="IGHV4-31*03_IGHD1-1*01_IGHJ4*02_IGHG1"),c("VDJC","sample_id")])
colnames(tmp_df)[2]<-"SL_tumor"
# add meta
tmp_2 <- merge(tmp_df,sample_meta[,c(3,5)],by="SL_tumor")
LT <- nrow(tmp_2)
#LH <- length(which(tmp_2$CH==1))
LH <- length(which(tmp_2$CH=="CH"))

matrix(c(LH,LT-LH,PH-LH,PT-LT-(PH-LH)),nrow=2,
                       dimnames= list(c("CH", "Not CH"),
                                      c("target_sample", "rest_sample")))

