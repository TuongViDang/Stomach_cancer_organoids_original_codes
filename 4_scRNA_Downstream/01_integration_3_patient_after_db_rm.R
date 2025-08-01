
# Single-cell RNA-seq analysis - integration
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(patchwork)

#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_path<-file.path(main_path,'rds')
rds_3_path <- file.path(main_path,'rds', 'three_patients')
dir.create(rds_3_path)

#Load the data

argu <- "doublet_rm_filtered_all_patients_qc.rds"
db_filtered_all_patients_qc<-readRDS(file.path(rds_path,argu))

#cell cycle markers
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


## Split seurat lists 
split_seurat <- db_filtered_all_patients_qc[9:17]
options(future.globals.maxSize = 1.5 * 1024^3)  # Set to 1.5 GB
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.genes, s.features=s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]] 
                                   #,method = "glmGamPoi",vars.to.regress = c('mitoRatio','phase')
  )
}
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

# Find best buddies - can take a while to run
# integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        # normalization.method = "SCT", 
                                        # anchor.features = integ_features,
                                        # reduction="cca"
                                        # )

# Integrate across conditions
# seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                  #  normalization.method = "SCT")
# Save integrated seurat object

# saveRDS(seurat_integrated, file.path(rds_3_path,  "three_patient_integrated_cca.rds"))


##Ref as Untreated samples
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        reduction="cca",
                                        reference = c(1,4,7)
                                        )

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# Save integrated seurat object
saveRDS(seurat_integrated, file.path(rds_3_path,"three_patient_integrated_cca_ref.rds"))