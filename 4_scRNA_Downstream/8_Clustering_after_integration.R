
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(patchwork)
library(RColorBrewer)

#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path<-file.path(main_path,'rds',"three_patients")
plot_3_path<-file.path(main_path,'plots_3_patients')
dir.create(plot_3_path)
plot_path_results_bf <- file.path(plot_3_path,'before_integration')
plot_path_results_aft <- file.path(plot_3_path,'after_integration')

Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)


functions_path <- '/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/downstream_analysis'

source(file.path(functions_path,'functions.R'))


#Load input file

seurat_integrated <- readRDS(file.path(rds_3_path,  paste0("three_patient_integrated_cca_ref.rds")))

seurat_integrated


#Before integration
DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 8 * 1024^3)  
seurat_integrated <- SCTransform(seurat_integrated, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = F)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

#change the name of sample
meta <- seurat_integrated@meta.data
meta <- meta %>%
  mutate(Sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ))
seurat_integrated@meta.data <- meta

seurat_integrated$Sample = factor(seurat_integrated$Sample, levels = c("OO77 Untreated",  "OO77 Ex vivo", "OO77 In vivo",
                                  "OO99 Untreated", "OO99 Ex vivo", "OO99 In vivo",
                                  "OO100 Untreated", "OO100 Ex vivo", "OO100 In vivo" ))

#Dim plot by Sample before integration
Idents(seurat_integrated) <- "Sample"
color_palette <- c(brewer.pal(12, "Set3")[1:8],  brewer.pal(8, "Set2")[4])

DimPlot(seurat_integrated, reduction = "umap",label=TRUE,label.size = 5, cols= color_palette, repel=T, raster=T)  + theme_void() +
        theme(legend.text = element_text( size = 12))

ggsave(file.path(plot_path_results_bf , 'bf_integration_by_patients.pdf'),width = 10,height = 8)




DefaultAssay(seurat_integrated)<-'integrated'
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = F)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

Idents(seurat_integrated) = "Sample"
plot <- DimPlot(seurat_integrated, reduction = "umap", cols = color_palette , raster=T, pt.size=1.5) + theme_void()
ggsave("After_integration_by_sample.pdf",path = plot_path_results_aft,width = 7.5,height = 4.7)

data = plot$data 
colnames(data)[3] = "Sample"

write.table(data, file.path(Figures_data, "Figure_Sup_7A_after_integration_by_sample.txt"), sep = "\t", quote = F)





#DimPlot(seurat_integrated, reduction = 'umap',group.by = 'Phase') + theme_void()
#ggsave("After_intergation_by_phase.pdf",path = plot_path_results_aft,width = 1000,height = 650,dpi = 150,units = 'px')



#Plot S score
FeaturePlot(seurat_integrated, features = c("S.Score"),raster=T, pt.size = 1.5,  min.cutoff = 0) +
  ggtitle("S Score") +
  theme_void() + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, color = "black"))+
  scale_color_distiller(palette = "Blues", direction = 1)
ggsave(filename=paste0("Featureplot_S_score.pdf"),path= plot_path_results_aft,width= 7,height= 4.7)

#Plot G2M score
FeaturePlot(seurat_integrated, features = c("G2M.Score"),raster=T, pt.size = 1.5,   min.cutoff = 0)  +
 ggtitle("G2M Score")+ 
 theme_void()+
 theme(plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, color = "black"))+
  scale_color_distiller(palette = "Greens",  direction = 1) 
ggsave(filename=paste0("Featureplot_G2M_score.pdf"),path= plot_path_results_aft,width= 7,height= 4.7)

#save data
data = seurat_integrated@meta.data %>% select(S.Score, G2M.Score)

umap_coords <- Embeddings(seurat_integrated, reduction = "umap")

all(rownames(umap_coords) == rownames(data))
data = cbind(umap_coords, data)

write.table(data, file.path(Figures_data, "Figure_Sup_7C_S.Score_G2M.Score.txt"), sep = "\t", quote = F)


#Plot cytokines
cytokines_genes = c("CXCL1", "CXCL2", "CXCL3", "CXCL8")

plot_data = list()
for (cytokine in cytokines_genes){
plot <- FeaturePlot(seurat_integrated, reduction = 'umap',features = cytokine, min.cutoff = -1, raster=T, pt.size = 1.5) +
        theme_void() +
        theme(plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, color = "black"))
plot
#ggsave(paste0("Featureplot_", cytokine,".pdf"),path = plot_path_results_aft,width = 7.2,height = 4.7)
plot_data[[cytokine]] = plot$data
}

for (cytokine in cytokines_genes[2:4]){
  data = cbind(plot_data[[cytokines_genes[1]]], plot_data[[cytokine]][4])
}

write.table(data, file.path(Figures_data, "Figure_Sup_7D_FeaturePlot_cytokines_expression.txt"), sep = "\t", quote = F)


#DO CLUSTERING
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:30)
print("Clustering finished!")
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c( 0.3,0.35, 0.4))
seurat_integrated$seurat_clusters<-seurat_integrated$integrated_snn_res.0.35


Idents(seurat_integrated) <- "seurat_clusters"


#set color

#set3_palette <- brewer.pal(12, "Set3")[3:12]
custom_palette <- c(  "#80B1D3", "limegreen", "#FCCDE5" ,
                  "#FDB462","#D9D9D9", "#CCEBC5" ,
                 "#BC80BD" , "khaki1",  "brown1")   

#rename clusters 
seurat_integrated@meta.data <- seurat_integrated@meta.data %>% mutate(seurat_clusters_new = case_when(seurat_clusters == 2 ~ 1,
                                                                                                      seurat_clusters == 3 ~ 2,
                                                                                                      seurat_clusters == 4 ~ 3,
                                                                                                      seurat_clusters == 5 ~ 4,
                                                                                                      seurat_clusters == 6 ~ 8,
                                                                                                      seurat_clusters == 0 ~ 6,
                                                                                                      seurat_clusters == 7 ~ 7,
                                                                                                      seurat_clusters == 1 ~ 5,
                                                                                                      seurat_clusters == 8 ~ 9))


seurat_integrated@meta.data$seurat_clusters_new  = factor(seurat_integrated@meta.data$seurat_clusters_new, levels = 1:9)
seurat_integrated@meta.data$seurat_clusters = seurat_integrated@meta.data$seurat_clusters_new
seurat_integrated@meta.data = seurat_integrated@meta.data %>% select(-seurat_clusters_new)


#Seurat object with transcriptomic cluster identity
seurat_integrated <- readRDS(file.path(rds_3_path, "clustered_0.35_integrated_seurat_ref.rds"))



Idents(seurat_integrated) = "seurat_clusters"
plot <- DimPlot(seurat_integrated,reduction='umap',label=TRUE, label.size = 8, cols= custom_palette, raster =T, pt.size = 1.5) + theme_void()+
       theme(legend.text = element_text(size = 17))
plot
ggsave('After_integration_by_cluster_0.35.pdf',path=plot_path_results_aft,width = 7,height = 4.4)

#Save Figures data
data = plot$data 
colnames(data)[3] = "Transcriptional_cluster"
write.table(data, file.path(Figures_data, "Figures_5B_after_intergration_by_cluster.txt"), sep = "'t", quote = F)




plot_integrated_clusters(seurat_integrated)
ggsave("Patient_contribution_each_cluster_0.35.pdf",path=plot_path_results_aft,width= 11.5,height= 6.5)



plot_cluster_contribution_each_sample(seurat_integrated)
ggsave("cluster_contribution_each_sample_0.35.pdf",path=plot_path_results_aft,width=10,height=6)

data = seurat_integrated@meta.data
data = data %>% select(Sample, seurat_clusters)

write.table(data, file.path(Figures_data, "Figure_Sup_7B_cluster_distribution_per_sample.txt"), sep = "\t", quote= F)

#save file

col_keep <- grep("^pANN|integrated", colnames(seurat_integrated@meta.data), invert = TRUE)

seurat_integrated@meta.data <- seurat_integrated@meta.data[, col_keep]

#saveRDS(seurat_integrated,file=file.path(rds_3_path, 'clustered_0.35_integrated_seurat_ref.rds'))

#reload the file
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))

