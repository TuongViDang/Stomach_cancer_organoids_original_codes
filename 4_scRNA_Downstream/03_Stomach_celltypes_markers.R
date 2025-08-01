# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(patchwork)
library(stringr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(HGNChelper)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(SingleCellExperiment)
library(scater)

#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path<-file.path(main_path,'rds',"three_patients")
plot_3_path<-file.path(main_path,'plots_3_patients')
dir.create(plot_3_path)
plot_path_results <- file.path(plot_3_path,"celltype_markers")
dir.create(plot_path_results)


#Figures data path
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)


#Load data
seurat_integrated<-readRDS(file.path(rds_3_path, "three_patient_integrated_cca_ref.rds"))

seurat_integrated$patient<- str_extract(seurat_integrated$orig.ident,"^[^[_-]]*")
seurat_integrated$treatment<-str_extract(seurat_integrated$orig.ident,"[^-]*$")

#Normalize and scale data
DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 10 * 1024^3)
seurat_integrated <- SCTransform(seurat_integrated,  vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = F)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

meta <- seurat_integrated@meta.data
meta <- meta %>%
  mutate(sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ))
seurat_integrated@meta.data <- meta


#color
set3_palette <- brewer.pal(12, "Set3")[c(1:8,10)]

Idents(seurat_integrated) <- "sample"

DimPlot(seurat_integrated, reduction = "umap",label=TRUE, cols= set3_palette, repel=T, raster = T, pt.size = 1.25) + theme_void()
ggsave('Bf_integration_by_patients.pdf',path = plot_path_results,width = 1000,height = 800,dpi = 150,units = 'px')





## Markers
Markers = c("EPCAM","PECAM1", "FBN1", "CD3E",
            "CDH1","CD34", "COL3A1", "PIM2")

# Generate FeaturePlots

plots <- FeaturePlot(seurat_integrated, features = Markers, ncol = 4, raster = TRUE, pt.size = 1.25)

# Increase legend text size for each plot
plots <- lapply(plots, function(p) {
  p + theme(legend.text = element_text(size = 17))  # Adjust size as needed
})


combined_plot <- wrap_plots(plots, ncol = 4) 
combined_plot
ggsave('Cell_type_markers.pdf',path = plot_path_results,width = 14,height = 6)

plot_data_list <- lapply(plots, function(p) p$data)

# Optionally, name each list element with the marker name
names(plot_data_list) <- Markers

df = plot_data_list[[1]] 
for (marker in Markers[2:8]){
 df = df %>% cbind(plot_data_list[[marker]][marker]) 
}
                                              
colnames(df)[3] = "sample"     
write.table(df, file.path(Figures_data, "Figure_Sup_1B_FeaturePlot_stomach_cell_types_markers.txt" ), sep = "\t", quote = F)