# Load libraries
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(ggvenn)
library(msigdbr)
library(ggpubr)
library(igraph)




#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
txt_files<-file.path(main_path,"txt_files")
rds_path<-file.path(main_path,'rds')
rds_path_3 <-file.path(main_path,'rds',"three_patients")
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_3_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_3_path,"NFKb_transcriptomic_clusters")
dir.create(plot_result_path)




#Load input file
#Specify the input argument
processed_file = file.path(rds_path,"doublet_rm_filtered_all_patients_qc.rds")
seurat_all <- readRDS(processed_file)


seurats <- seurat_all[c("OO100-no","OO100-vitro","OO100-vivo","OO77-no","OO77-vitro","OO77-vivo","OO99-no","OO99-vitro","OO99-vivo")]
#cell cycle markers
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Prepare Hallmark pathways genes list (version 2024)
Hallmark_files = list.files("/group/poetsch_projects/poetsch_sc/Driver_predict/cancer_gene_lists/HALLMARK", full.names = T )


Hallmark_df_build <- function(file){
  lines <- readLines(file)
  start_line <- grep("^GENE_SYMBOLS\t", lines)
  end_line <- grep("^FOUNDER_NAMES\t", lines)
  gene_symbols_lines <- lines[start_line:end_line-1]
  gene_symbols_combined <- paste(gene_symbols_lines, collapse = " ")
  gene_symbols <- unlist(strsplit(gsub("GENE_SYMBOLS\t|FOUNDER_NAMES\t", "", gene_symbols_combined), ","))
  gene_symbols <- trimws(gene_symbols)
  gene_symbols
}

Hallmark_genes_lists <- lapply(Hallmark_files, Hallmark_df_build )
Pathway_names = sapply(Hallmark_files, function(file){ str_extract(basename(file), "(?<=_)[^.]+")})
names(Hallmark_genes_lists) = Pathway_names 



seurats_list = list()
for (i in seq_along(seurats)) {
seurat = seurats[[i]]
DefaultAssay(seurat) = "RNA"
seurat <- NormalizeData(seurat, verbose = TRUE)
#Add cell cycle score
seurat <- CellCycleScoring(seurat, g2m.features=g2m.genes, s.features=s.genes)

# Add Module Scores for these pathways

seurat <- AddModuleScore(seurat, features = Hallmark_genes_lists )

seurats_list[[i]] = seurat
}

seurat_merged = seurats_list[[1]]
for (i in 2:9){
   seurat_merged = merge(seurat_merged ,  seurats_list[[i]])
}



meta <- seurat_merged@meta.data
meta <- meta %>%
  mutate(sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ),
  patient = 
   treatment = sub("^[^ ]* ","",sample))


colnames(meta)[14:21] = Pathway_names 
seurat_merged@meta.data <- meta[,c("orig.ident","sample","treatment", "TNFA_SIGNALING_VIA_NFKB")]


#Seurat object with transcriptomic cluster identity
seurat_integrated <- readRDS(file.path(rds_3_path, "clustered_0.35_integrated_seurat_ref.rds"))

#check identical cell CB
all(rownames(meta) == rownames(seurat_merged@meta.data))

#Add information of Nfkb
meta_new <- cbind(seurat_merged@meta.data , seurat_clusters = seurat_integrated$seurat_clusters)
seurat_integrated@meta.data = meta_new

cluster_sizes_NFKB = meta_new %>% group_by(treatment, seurat_clusters) %>% 
                                  summarize(count = n(),
                                  mean_NFKb = mean(TNFA_SIGNALING_VIA_NFKB))



context_df_dict <- readRDS(file.path(rds_3_path ,"Liana_context_df_dict.rds"))





#Untreated 
cluster_sizes_NFKB_Untreated = cluster_sizes_NFKB %>% filter( treatment  == "Untreated")

untreateds = c(1,4,7)
dict_sub <- context_df_dict[untreateds]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))
top_100$source = factor(top_100$source, levels = 1:9)
top_100$target = factor(top_100$target, levels = 1:9)
df = as.matrix(table(top_100$source,top_100$target))


# Create the igraph object from the interaction matrix
g_Untreated <- graph_from_adjacency_matrix(as.matrix(df), mode = "undirected", weighted = TRUE, diag = FALSE)


# Assign node attributes
V(g_Untreated)$nfkb_score <- cluster_sizes_NFKB_Untreated$mean_NFKb

# Use continuous color scale for NFkB score
nfkb_min <- min(cluster_sizes_NFKB$mean_NFKb)
nfkb_max <- max(cluster_sizes_NFKB$mean_NFKb)

# Use a color ramp from light to dark red
color_palette <- colorRampPalette(brewer.pal(9, "Reds"))(100)

# Scale NFkB score 
color_index_untreated <- as.numeric(cut(V(g_Untreated)$nfkb_score, breaks=seq(nfkb_min, nfkb_max, length.out=101), include.lowest = TRUE))

#check color with score
print(data.frame(cluster=V(g_Untreated)$name,
                  nfkb=V(g_Untreated)$nfkb_score,
                  color_index= color_index_untreated ))

# Assign colors
V(g_Untreated)$color <- color_palette[color_index_untreated ]

#Layout as coordinate alculate from edge weight
#layout_Untreated <- layout_with_fr(g_Untreated, weights = E(g_Untreated)$weight)

# Plot and save
pdf(file.path(plot_result_path, "NFkB_network_untreated.pdf"), width = 13, height = 8)
set.seed(123) 
plot(g_Untreated,
     #layout = layout_Untreated,
     vertex.label = V(g_Untreated)$name,
     vertex.label.color = "black",
     vertex.label.cex = 4,
     edge.width = E(g_Untreated)$weight / 1.5,s
     vertex.color = V(g_Untreated)$color,
     vertex.size = 25,   # constant size for all nodes
     main = "Untreated 3 patients")

# Add color legend for continuous NFkB activity
legend_gradient <- rev(color_palette)  # Reverse the color gradient
legend_image <- as.raster(matrix(legend_gradient, ncol=1))
rasterImage(legend_image, 1.2, 0.2, 1.1, 0.8)
text(1.25, 0.2, sprintf("%.2f", nfkb_min), adj=0)
text(1.25, 0.8, sprintf("%.2f", nfkb_max), adj=0)

dev.off()




### Ex vivo
cluster_sizes_NFKB_Ex_vivo = cluster_sizes_NFKB %>% filter( treatment  == "Ex vivo")

Ex_vivo = c(2,5,8)
dict_sub <- context_df_dict[Ex_vivo]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))
top_100$source = factor(top_100$source, levels = 1:9)
top_100$target = factor(top_100$target, levels = 1:9)
df = as.matrix(table(top_100$source,top_100$target))


# Create the igraph object from the interaction matrix
g_Ex_vivo <- graph_from_adjacency_matrix(as.matrix(df), mode = "undirected", weighted = TRUE, diag = FALSE)

# Assign node attributes
V(g_Ex_vivo)$nfkb_score <- cluster_sizes_NFKB_Ex_vivo$mean_NFKb



# Scale NFkB score to 1-100 (use Untreated as reference)
color_index_ex_vivo <- as.numeric(cut(V(g_Ex_vivo)$nfkb_score, 
                          breaks=seq(nfkb_min, nfkb_max, length.out=101), include.lowest = TRUE))

#check color with score
print(data.frame(cluster=V(g_Ex_vivo)$name,
                  nfkb=V(g_Ex_vivo)$nfkb_score,
                  color_index=color_index_ex_vivo))

# Assign colors
V(g_Ex_vivo)$color <- color_palette[color_index_ex_vivo]

#Layout as coordinate alculate from edge weight
#layout_Ex_vivo <- layout_with_fr(g_Ex_vivo, weights = E(g_Ex_vivo)$weight)

# Plot and save
pdf(file.path(plot_result_path, "NFkB_network_Ex_vivo.pdf"), width = 13, height = 8)
set.seed(123) 
plot(g_Ex_vivo,
   # layout = layout_Ex_vivo,
     vertex.label = V(g_Ex_vivo)$name,
     vertex.label.color = "black",
     vertex.label.cex = 4,
     edge.width = E(g_Ex_vivo)$weight / 1.5,
     vertex.color = V(g_Ex_vivo)$color,
     vertex.size = 25,   # constant size for all nodes
     main = "Ex vivo 3 patients")

# Add color legend for continuous NFkB activity
legend_gradient <- rev(color_palette)  # Reverse the color gradient
legend_image <- as.raster(matrix(legend_gradient, ncol=1))
rasterImage(legend_image, 1.2, 0.2, 1.1, 0.8)
text(1.25, 0.2, sprintf("%.2f", nfkb_min), adj=0, cex = 1.5)
text(1.25, 0.8, sprintf("%.2f", nfkb_max), adj=0, cex = 1.5)

dev.off()


#In vivo

cluster_sizes_NFKB_In_vivo = cluster_sizes_NFKB %>% filter( treatment  == "In vivo")

In_vivo = c(3,6,9)
dict_sub <- context_df_dict[In_vivo]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))
top_100$source = factor(top_200$source, levels = 1:9)
top_100$target = factor(top_200$target, levels = 1:9)
df = as.matrix(table(top_100$source,top_100$target))


# Create the igraph object from the interaction matrix
g_In_vivo <- graph_from_adjacency_matrix(as.matrix(df), mode = "undirected", weighted = TRUE, diag = FALSE)

# Assign node attributes
V(g_In_vivo)$nfkb_score <- cluster_sizes_NFKB_In_vivo$mean_NFKb



# Scale NFkB score to 1-100 (use Untreated as reference)
# Cut based on the Untreated range
color_index_in_vivo <- as.numeric(cut(V(g_In_vivo)$nfkb_score , 
                                      breaks=seq(nfkb_min, nfkb_max, length.out=101),
                                      include.lowest=TRUE))
#check color with score
print(data.frame(cluster=V(g_In_vivo)$name,
                  nfkb=V(g_In_vivo)$nfkb_score,
                  color_index= color_index_in_vivo ))

# Assign colors
V(g_In_vivo)$color <- color_palette[color_index_in_vivo]

#Layout as coordinate alculate from edge weight
#layout_In_vivo <- layout_with_fr(g_In_vivo, weights = E(g_In_vivo)$weight)

# Plot and save
pdf(file.path(plot_result_path, "NFkB_network_In_vivo.pdf"), width = 13, height = 8)
set.seed(123) 
plot(g_In_vivo,
   # layout = layout_In_vivo,
     vertex.label = V(g_In_vivo)$name,
     vertex.label.color = "black",
     vertex.label.cex = 4,
     edge.width = E(g_In_vivo)$weight/1.5 ,
     vertex.color = V(g_In_vivo)$color,
     vertex.size = 25,   # constant size for all nodes
     main = "In vivo 3 patients")

# Add color legend for continuous NFkB activity
legend_gradient <- rev(color_palette)  # Reverse the color gradient
legend_image <- as.raster(matrix(legend_gradient, ncol=1))
rasterImage(legend_image, 1.2, 0.2, 1.1, 0.8)
text(1.25, 0.2, sprintf("%.2f", nfkb_min), adj=0, cex = 1.5)
text(1.25, 0.8, sprintf("%.2f", nfkb_max), adj=0, cex = 1.5)

dev.off()



#Line Plot pathway score 
meta 
seurat_integrated = meta_99_no_vivo %>% group_by(Clone, treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB))