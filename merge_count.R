setwd("~/Projects/TCR/cell_enrichment/")

library(GenomicFeatures)


list <- read.table("blo_count_list",sep = "\t",header=F)


txdb <- makeTxDbFromGFF("gencode.v38lift37.annotation.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))

exonic.gene.sizes$gene_id <- rownames(exonic.gene.sizes)

colnames(exonic.gene.sizes)[1] <- "exo_size"

gtf.gr <- rtracklayer::import("gencode.v38lift37.annotation.gtf")
gtf.df <- as.data.frame(gtf.gr)
genes <- unique(gtf.df[which(gtf.df$gene_type=="protein_coding"&
                             gtf.df$type=="gene"),c("gene_id","gene_name")])

gene_exo_size <- merge(genes,exonic.gene.sizes,by="gene_id")


#write.table(gene_exo_size,"gene_exo_size_hg19.txt",sep="\t",quote=F,row.names = F)

output <- data.frame(gene = tmp_2$gene_name)
output_normal <- data.frame(gene=tmp_2$gene_name)
for (i in 1:nrow(list)) {
  print(i)
  tmp_in_1 <- read.table(paste0("blo_count/",list$V1[i],"_fc"),header=T,sep="\t")
  tmp_in <- tmp_in_1[,c(1,7)]
  colnames(tmp_in)<- c("gene_id",list$V1[i])
  tmp_2 <- merge(tmp_in,gene_exo_size,by="gene_id")
  # raw_count
  output <- cbind(output,tmp_2[,c(2),drop=F])

  # nomralized TPM
  tmp_2$rpk <- tmp_2[,2]/tmp_2[,4]
  scale_f <- sum(tmp_2$rpk)/1000000

  tmp_2$tmp <- tmp_2$rpk/scale_f
  colnames(tmp_2)[c(2,6)] <- c("tmp",colnames(tmp_in)[2])
  output_normal <- cbind(output_normal,round(tmp_2[,c(6),drop=F],digit=3))
}

write.table(output,"blo_expression_matrix.txt",sep="\t",quote=F,row.names = F)
write.table(output_normal,"blo_expression_matrix_TPM.txt",sep="\t",quote=F,row.names = F)




