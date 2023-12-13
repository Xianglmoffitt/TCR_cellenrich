setwd("~/Projects/TCR/")

in_tmp <- read.table("cell_enrichment/blo_expression_matrix_TPM.txt",sep="\t",header=T)
mm_meta_in <- read.table("metadata/blo_tumor_meta.txt",sep="\t",header=T)

# separete cancer (BRE - Breast Cancer, GYN - Ovarian Cancer, THO - Lung Cancer,)
mm_meta <- mm_meta_in[which(mm_meta_in$Cancer_type=="GYN - Ovarian Cancer"),]

out_tmp <- in_tmp[,c(1,which(colnames(in_tmp) %in% mm_meta$SL_tumor))]


write.table(out_tmp,"l_expression_matrix_TPM.txt",sep="\t",row.names = F,quote = F)
