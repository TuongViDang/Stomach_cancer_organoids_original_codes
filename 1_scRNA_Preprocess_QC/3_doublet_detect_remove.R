
# Single-cell RNA-seq analysis - detect and remove doublet 
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(patchwork)
library(DoubletFinder)

#Specify the path
main_path='/group/poetsch_projects/Projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
plot_path=file.path(main_path,'plots')
plot_path_doublet<-file.path(main_path,'plots','doublet')
txt_files<-file.path(main_path,'txt_files')
rds_path<-file.path(main_path,'rds')

#Specify argument
CA <- commandArgs(trailingOnly = TRUE)
#Load the data
filtered_all_patients_qc<-readRDS(file.path(rds_path,paste0(argu,'.rds')))



## Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_all_patients_qc, split.by = "orig.ident")


Find_Doublet<-function(object,expected_doublet){
  object <- SCTransform(object)
  object <- RunPCA(object)
  object <- RunUMAP(object,dims=1:30)
  
  sweep.res.x <- paramSweep_v3(object, PCs = 1:10, sct = TRUE)
  sweep.x <- summarizeSweep(sweep.res.x, GT = FALSE)
  bcmvn_x <- find.pK(sweep.x)
  pK <- as.numeric(levels(bcmvn_x$pK)[which.max(bcmvn_x$BCmetric)])
  nExp <- round(ncol(object) * expected_doublet)
  
  object <- doubletFinder_v3(object, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:50, sct = TRUE)
  
  #Plot to vidualize
  DF.name=colnames(object@meta.data)[grepl("DF.classification", colnames(object@meta.data))]
  cowplot::plot_grid(ncol = 2, DimPlot(object, group.by = "orig.ident") + NoAxes(),
                     DimPlot(object, group.by = DF.name) + NoAxes())
  
sample_name<-unique(object$orig.ident)
ggsave(filename = paste0(sample_name,'_db_detection.pdf'),path = plot_path_doublet,width = 10000,height = 5000,dpi = 1500,units = 'px')  
  #remove doublet
  object@meta.data<-object@meta.data%>%rename(doublet_detect=DF.name)
  object<-subset(x = object,subset = (doublet_detect=='Singlet'))
  object
}

#Specify expected proportion of doublet (10X) based on number of barcode
expected_doublet<-c(rep(0.032,3),0.008,0.04,0.032,0.08,0.072,0.048,0.056,0.08,0.048,0.072,0.04,0.056,0.04,0.032)
expected_doublet

#Do the detection and removal
db_split_seurat <- purrr::map2(split_seurat,expected_doublet,Find_Doublet)

#Save the processed lists of seurat object
saveRDS(db_split_seurat,file.path(rds_path,paste0('doublet_rm_',argu,'.rds')))
