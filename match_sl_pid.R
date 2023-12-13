setwd("~/Projects/TCR/")

main_in <- read.table("metadata/Avatar_breast_lung_ovarian_504_CHIPstatus_withGene_withCancertype.txt",header=T,sep="\t")
normal_pid <- read.table("metadata/blo_normal_pid",header=F,sep="\t")
tumor_pid <- read.table("metadata/blo_tumor_pid",header=F,sep="\t")


#-----------
# mm
#-----------
main_in$SL_normal <- gsub("_normal","",main_in$patient_id)

normal_tumor_pid<- merge(normal_pid,tumor_pid,by="V2")

colnames(normal_tumor_pid) <-c("PID","SL_normal","SL_tumor")

ch_out <- merge(main_in,normal_tumor_pid,by="SL_normal")

#write.table(ch_out,"metadata/tumor_normal_pid_sl_CH.txt",sep="\t",quote=F,row.names = F)


#-----------
# blo
#-----------

colnames(main_in)[1] <- "SL_normal"

normal_tumor_pid<- merge(normal_pid,tumor_pid,by="V2")

colnames(normal_tumor_pid) <-c("PID","SL_normal","SL_tumor")

ch_out <- merge(main_in,normal_tumor_pid,by="SL_normal")

write.table(ch_out,"metadata/blo_tumor_meta.txt",sep="\t",quote=F,row.names = F)


