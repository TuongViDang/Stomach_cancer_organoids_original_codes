## Load libraries
library(Seurat)
library(tidyverse)
library(HGNChelper)


#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path<-file.path(main_path,'rds','three_patients')
txt_files_3 <-file.path(main_path,'txt_files','three_patients')
dir.create(txt_files_3)


#Load file
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))

#normalize and scale RNA assay
DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 8 * 1024^3)  
seurat_integrated <- NormalizeData(seurat_integrated)
seurat_integrated <- FindVariableFeatures(seurat_integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_integrated)
seurat_integrated <- ScaleData(seurat_integrated, features = all.genes)



Idents(seurat_integrated) = "seurat_clusters"

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  annotations <- read.csv(file.path(main_path,"txt_files/annotation.txt"))
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) 
}
clusters <- levels(seurat_integrated$seurat_clusters)
clusters
Markers_all_clusters_list<-purrr::map(clusters, get_conserved)

#add cluster names for each element in list
length(Markers_all_clusters_list)
names(Markers_all_clusters_list)<-clusters

#bind row to create dataframe and add a column name 'cluster'
Markers_all_clusters  <- dplyr::bind_rows(Markers_all_clusters_list,.id = 'cluster')

Markers_all_clusters %>% filter(cluster == 9) %>% head()
#save the table
write.table(Markers_all_clusters,file.path(txt_files_3,paste0('markers_0.35_ref.txt')),sep = '\t',quote = F)
