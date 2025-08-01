# Load libraries
library(Seurat)
library(tidyverse)


main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_path<-file.path(main_path,'rds')




seurat_integrated<-readRDS(file.path(rds_path,"clustered_granuled_integrated_doublet_rm_filtered_all_patients_qc_cca.rds"))
seurat_list<-SplitObject(seurat_integrated,split.by = "orig.ident")

patients=c("OO100","OO77")
for (patient in patients){
organoids = paste(patient,c("no", "vitro","vivo"), sep = "-")

for (organoid in organoids){
seurat<-seurat_list[[organoid]]

BC <- gsub(".*?_","",rownames(seurat@meta.data))

#Write the BC file to scReadCounts
out_dir=file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scReadCounts/",patient,organoid)
dir.create(out_dir,recursive = T)

write.table(BC,file.path(out_dir,"barcodes.tsv"),sep="\t",quote=F,row.names = F)
}
}
