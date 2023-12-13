#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

all_report <- read.table("~/Projects/TCR/trust4_out_all_report.tsv",sep="\t",header=T)
sample_meta <- read.table("~/Projects/TCR/metadata/mm_tumor_meta.txt",sep="\t",header=T)

fisher_test_by_row <- function(row) {
  tmp_df <- unique(all_report[which(all_report$CDR3nt==row[1]),c("CDR3nt","sample_id")])
  colnames(tmp_df)[2]<-"SL_tumor"
  # add meta
  tmp_2 <- merge(tmp_df,sample_meta[,c(3,5)],by="SL_tumor")
  LT <- nrow(tmp_2)
  LH <- length(which(tmp_2$CH=="CH"))

  contigency_m <- matrix(c(LH,LT-LH,PH-LH,PT-LT-(PH-LH)),nrow=2,
                         dimnames= list(c("CH", "Not CH"),
                                        c("target_sample", "rest_sample")))
  ft <- fisher.test(contigency_m)
  return(ft$p.value)
}


target_list <- read.table(args[1],sep="\t",header=F)

PT <- nrow(sample_meta)
PH <- length(which(sample_meta$CH=="CH"))


pval<- apply(target_list, MARGIN = 1, fisher_test_by_row)

df_out <- data.frame(target = target_list,
                     pval=pval)

write.table(df_out, file=args[2], row.names=FALSE,quote=F)

