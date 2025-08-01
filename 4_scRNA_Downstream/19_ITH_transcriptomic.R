library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(GeneNMF)
library(ggrastr)

#Specify the path
Integrate_dir="/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/ITH"
plot_dir <- file.path(Integrate_dir,"plot" )
dir.create(plot_path, recursive = T)
rds_path <- '/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/rds'

#Figures data path
main_path <-'/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)



argu <- "integrated_doublet_rm_filtered_all_patients_qc_cca"


seurat <- readRDS(file.path(rds_path,paste0(argu,'.rds')))
DefaultAssay(seurat) = "RNA"
options(future.globals.maxSize = 10 * 1024^3)  # 8 GB

seurat <- SCTransform(seurat , vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat <- RunPCA(seurat, npcs = 30, verbose = F)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30, verbose = F)

#change the name of sample
meta <- seurat@meta.data
meta <- meta %>%
  mutate(sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ))
seurat@meta.data <- meta



pca_embeddings <- Embeddings(seurat, reduction = "pca")[, 1:30]
pca_df <- data.frame(pca_embeddings, Sample = seurat@meta.data$sample)

library(stats)  # For dist function


sum_distance_per_sample <- function(pca_df) {
  # Ensure the Sample column exists
  if (!"Sample" %in% colnames(pca_df)) {
    stop("The data frame must contain a 'Sample' column.")
  }
  
  # Split data by sample
  samples <- unique(pca_df$Sample)
  norm_distances <- numeric(length(samples))
  names(norm_distances) <- samples
  
  for (sample in samples) {
    # Subset PCA data for the current sample
    sample_data <- pca_df[pca_df$Sample == sample, 1:30]
    
    # Compute pairwise Euclidean distances
    dists <- as.matrix(dist(sample_data, method = "euclidean"))
    
    # Sum all distances in the matrix (upper triangle to avoid double counting)
    total_distance <- sum(dists) / 2
    
    # Normalize by the number of cells in the sample
    num_cells <- nrow(sample_data)
    
    if (num_cells > 1) {
      norm_distances[sample] <- total_distance / (num_cells * (num_cells-1)* 1/2)
    } else {
      norm_distances[sample] <- 0  # If only one cell, distance is zero
    }
  }
  
  return(norm_distances)
}


norm_distances <- sum_distance_per_sample(pca_df)
norm_distances_df <- as.data.frame (norm_distances[9:17])

norm_distances_df <- norm_distances_df %>% tibble::rownames_to_column("sample")
colnames(norm_distances_df )[2] = "ITH"

norm_distances_df <- norm_distances_df %>% 
  mutate(
    patient = sub(" .*", "", sample),   # Extract text before the first space
    treatment = sub(".*? ", "", sample) # Extract text after the first space
  )

All <- norm_distances_df
All$patient = factor(All$patient, levels = c("OO77", "OO99", "OO100"))
# Create a new column for labeling only tumor samples
All$label <- ifelse(All$treatment == "Untreated", All$patient, NA)
#Ex vivo
All_ex <- All %>% filter(treatment != "In vivo")
All_ex$treatment = factor(All_ex$treatment, levels = c("Untreated","Ex vivo"))

#t test
t_test <- t.test(ITH ~ treatment, data = All_ex, paired = TRUE, alternative = "greater")
p_value <-  t_test$p.value


#plot

colors = c("deepskyblue2","tomato", "goldenrod3", "darkcyan","indianred1")
ggplot(All_ex, aes(x = treatment, y =  ITH , group = patient)) +
  geom_point(aes(colour = treatment), size = 2, position = position_dodge(width = 0)) +
  geom_line(aes(group = patient, colour = patient), linewidth = 2, alpha = 0.5, position = position_dodge(width = 0)) +
  scale_y_continuous(limits = c(0,max(All_ex$ITH + 60)),breaks = c(0,50,100, 150))+
  scale_colour_manual(values = colors)+
  ylab("Transcriptomic ITH") +
  #geom_text(aes(label =  label), vjust = -0.5, position = position_dodge(width = 0), size = 6) +
  theme_classic() +
  theme( legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
         panel.spacing = unit(1, "lines"))  +
   # Add the p-value as annotation
  annotate("text", x = 1.5, y = 140, label = paste("p =", format.pval(p_value, digits = 3)), size = 8)

  ggsave("transcriptomic_ITH_ex_vivo.pdf",path=plot_dir,width= 650 ,height= 650,dpi=150,units="px")

#In vivo
All_in <- All %>% filter(treatment != "Ex vivo")
All_in$treatment = factor(All_in$treatment, levels = c("Untreated","In vivo"))

#t test
t_test <- t.test(ITH ~ treatment, data = All_in, paired = TRUE, alternative = "greater")
p_value <-  t_test$p.value


colors = c("deepskyblue2","plum", "goldenrod3", "darkcyan","indianred1")
#Plot
ggplot(All_in, aes(x = treatment, y = ITH , group = patient)) +
  geom_point(aes(colour = treatment), size = 2, position = position_dodge(width = 0)) +
  geom_line(aes(group = patient, colour = patient), linewidth = 2, alpha = 0.5, position = position_dodge(width = 0)) +
  scale_y_continuous(limits = c(-5,max(All_in$ITH + 50)), breaks = c(0, 50,100))+
  scale_colour_manual(values = colors)+
  ylab("Transcriptomic ITH") +
  #geom_text(aes(label = label), vjust = -0.5, position = position_dodge(width = 0.5), size = 6) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
         panel.spacing = unit(1, "lines"))  +
  # Add the p-value as annotation
  annotate("text", x = 1.5, y = 135, label = paste("p =", format.pval(p_value, digits = 3)), size = 8)

  ggsave("Transcriptomic_ITH_in_vivo.pdf",path=plot_dir,width= 650,height= 650,dpi=150,units="px")



#Save Figure data
All = All %>% select(sample, patient, treatment, ITH)
write.table( All, file.path(Figures_data, "Figure_Sup_8A_Transcriptomic_ITH.txt"), sep = "\t", quote = F)