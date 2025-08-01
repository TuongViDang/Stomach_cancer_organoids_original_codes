

# Load libraries
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(patchwork)
library(RColorBrewer)
library(factoextra)
library(pheatmap)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(ggvenn)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)

#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_path_3 <-file.path(main_path,'rds',"three_patients")
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_3_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_3_path,"Transcriptomic_dissimilarity_by_patient")
dir.create(plot_result_path)

#Load data
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))
seurat_integrated


#Before integration
DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
seurat_integrated <- SCTransform(seurat_integrated, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)
seurat_integrated  <- RunPCA(seurat_integrated, npcs = 30, )
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

seurat_integrated$patient<- str_extract(seurat_integrated$orig.ident,"^[^[_-]]*")
seurat_integrated$treatment<-str_extract(seurat_integrated$orig.ident,"[^-]*$")
meta <- seurat_integrated@meta.data
meta <- meta %>%
  mutate(sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ), 
   treatment = case_when(
    treatment == "no" ~ "Untreated",
    treatment == "vitro" ~ "Ex vivo",
    treatment == "vivo" ~ "In vivo"
   ))
seurat_integrated@meta.data <- meta


#Transform to SCE for bulk analysis
DefaultAssay(seurat_integrated)<-'RNA'
counts_matrix<- GetAssayData(object = seurat_integrated, slot = "counts")
coldata<-seurat_integrated@meta.data
#umap_embeddings<-Embeddings(object = seurat_integrated, reduction = "umap")

sce <- SingleCellExperiment(assays = list(counts=counts_matrix), colData=coldata)
#reduceDim(sce,"umap")<-umap_embeddings


summed <- aggregateAcrossCells(sce, id=colData(sce)[,"sample"])

#Draw the transcriptomic dissimilarity 
bulk_df<-as.data.frame(assay(summed,"counts"))

dist_matrix<- 1- cor(bulk_df,method="spearman")

dist_matrix_bis <- dist_matrix * 100
clean_dmat <- as.data.frame(rbind(dist_matrix_bis[3,1:2],
                      dist_matrix_bis[6,4:5],
                      dist_matrix_bis[9,7:8]))

rownames(clean_dmat)<-c("OO100","OO77","OO99")
colnames(clean_dmat)<-c("untreated_vs_in_vitro","untreated_vs_in_vivo")
clean_dmat$patients <-rownames(clean_dmat)
clean_dmat_long <- clean_dmat%>%pivot_longer(cols=c(1:2),names_to="treatment",values_to="distance")
clean_dmat_long$patient <- factor(clean_dmat_long$patient , levels = cc("OO77","OO99","OO100"))

colors = c( "goldenrod3", "darkcyan","indianred1","tomato", "violet")

ggplot(clean_dmat_long, aes(x=treatment, y=distance, group=patients)) +
  geom_point(aes(colour = treatment),size= 2, position=position_dodge(width=0)) +
  geom_line( aes(colour = patients), linewidth= 2, alpha=0.5, position=position_dodge(width=0)) +
  scale_x_discrete(labels = c("untreated_vs_in_vitro" = " Ex vivo", "untreated_vs_in_vivo" = "In vivo")) +
  scale_colour_manual(values = colors)+
  ylim(0,8)+
  ylab("Transcriptomic dissimilarity") +
  #+geom_text(aes(label = patients), vjust = -0.3, position = position_dodge(width = 0.55),size= 5) +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x.bottom=element_blank(),
        axis.text.x.bottom=element_text(size= 16),
        axis.title.y.left=element_text(size= 19),
         axis.text.y.left=element_text(size= 16))

  ggsave("distance_vitro_vivo.pdf",path= plot_result_path,width= 650,height= 650,dpi=150,units="px")
