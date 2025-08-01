# Load libraries
library(tidyverse)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(RCurl)
library(cowplot)
library(patchwork)
library(vegan)
library(pheatmap)
library(ggpubr)
library(magick)
library(variancePartition)
library(edgeR)


Integrate_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA")
plot_dir = file.path(Integrate_dir,"plots","Clone_Transcriptomic")
dir.create(plot_dir, recursive = T)


main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)

#calculate variance explained by Clone for each sample
variance_explained <- function(sample, seurat_patient){
    seurat <- subset(seurat_patient, subset = treatment == sample)
    expr_mat <- seurat$RNA$counts
    metadata <- seurat@meta.data %>% select("Clone")
    # Normalize counts using voom
    dge <- DGEList(counts = expr_mat)
    dge <- calcNormFactors(dge)
    expr_voomed <- voom(dge)$E  # Log2-normalized expression

    # Assuming expr_voomed is your log2-transformed expression data
    gene_variances <- apply(expr_voomed, 1, var)  # Variance for each gene (rows = genes, columns = samples)

    # Order the genes by variance (highest variance first)
    ordered_genes_var <- sort(gene_variances, decreasing = TRUE)

    # Select the top N most variated genes (e.g., top 2000 genes)
    top_genes <- names(ordered_genes_var)[1:2000]
    # Subset the expression data to include only the top variated genes
    expr_voomed_top <- expr_voomed[top_genes, ]

    #Formular: both patient and treatment as fixed effect
    formula <- ~ as.factor(Clone) 
    varPart <- fitExtractVarPartModel(expr_voomed_top, formula, metadata)
    colnames(varPart) = c( "Clone","residuals")
    varPart$Clone
}


variance_patient <- function(patient){
    rds_dir = file.path(Integrate_dir,patient, "tree30_change_name","rds")
    seurat_patient <- readRDS( file.path(rds_dir,paste0(patient,"three_samples_final_Clones.rds")))
    samples = c("Untreated","Ex vivo","In vivo")
    results_list <- map(samples,~ variance_explained(.x, seurat_patient = seurat_patient))
    names(results_list) = samples
    result_df = data.frame(Patient = rep(patient, 6000),
                           Treatment = rep(c("Untreated", "Ex vivo", "In vivo"), each = 2000),
                           Variance = unlist(results_list))
    result_df 
}

patients = c("OO77", "OO99", "OO100")

All_df = map_dfr(patients, variance_patient)
All_df = All_df %>% filter(Treatment != "Ex vivo")

All_df$Treatment = factor(All_df$Treatment, levels = c("Untreated",  "In vivo"))
All_df$Patient = factor(All_df$Patient , levels = patients)


cols = c("cadetblue2","plum")
All_df %>%
  ggplot(aes(x = Treatment, y = Variance, fill = Treatment)) +
    geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
    geom_boxplot( width = 0.1, outlier.size = 0.1)+
    geom_violin(
    trim = TRUE, 
    width = 0.5, 
    alpha = 0.5, 
    color = "black"
  ) +  
    scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(-0.05,0.43), breaks = seq(0,0.4,0.2))+
    ylab("% Variance Explained by genetic clone") +
    theme_classic() +
    theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17, color = "black"),
        strip.text = element_text(size = 17, color = "black"),
        panel.spacing = unit(1.5, "lines") ,
         strip.background = element_blank(),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 17,color = "black" ))  +
    facet_wrap(vars(Patient),  nrow = 3)
ggsave(file.path(plot_dir, "variance_clone_transcriptomic.pdf"), width = 4, height = 7)

rownames(All_df) <- NULL
write.table(All_df, file.path(Figures_data, "Figure_Sup_8B_Variance_Clone.txt"), sep = "\t", quote= F)


#Perform statistic

stat_results <- All_df %>%
  group_by(Patient) %>%
  summarise(
    ttest_ex_vivo = list(t.test(Variance[Treatment == "Untreated"], 
                                Variance[Treatment == "Ex vivo"], 
                                alternative = "greater",
                                var.equal = TRUE)),
    ttest_in_vivo = list(t.test(Variance[Treatment == "Untreated"], 
                                Variance[Treatment == "In vivo"], 
                                alternative = "greater",
                                var.equal = TRUE))
  ) %>%
  mutate(
    p_value_ex_vivo = sapply(ttest_ex_vivo, function(x) x$p.value),
    p_value_in_vivo = sapply(ttest_in_vivo, function(x) x$p.value)
  ) %>%
    mutate(
    sig_level_untreated_vs_exvivo = case_when(
      p_value_ex_vivo < 0.001 ~ "***",
      p_value_ex_vivo < 0.01 ~ "**",
      p_value_ex_vivo < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_level_untreated_vs_invivo = case_when(
      p_value_in_vivo < 0.001 ~ "***",
      p_value_in_vivo < 0.01 ~ "**",
      p_value_in_vivo < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>% select(Patient, sig_level_untreated_vs_exvivo ,sig_level_untreated_vs_invivo)



#Intra clone heterogeneity

intra_clone_ITH <- function(patient){
    rds_dir = file.path(main_dir,patient, "tree30_change_name","rds")
    seurat_patient <- readRDS( file.path(rds_dir,paste0(patient,"three_samples_final_Clones.rds")))
    #change the name of sample
    meta <- seurat_patient@meta.data
    meta <- meta %>% mutate(sample = case_when(
      grepl("-no", orig.ident) ~ sub(paste0(patient,"-no"), "Untreated", orig.ident),
      grepl("-vitro", orig.ident) ~ sub(paste0(patient,"-vitro"), "Ex vivo", orig.ident),
      grepl("-vivo", orig.ident) ~ sub(paste0(patient,"-vivo"), "In vivo", orig.ident),
      TRUE ~ orig.ident))
    seurat_patient@meta.data <- meta
    pca_embeddings <- Embeddings(seurat_patient, reduction = "pca")[, 1:30]
    pca_df <- data.frame(pca_embeddings,
                         Sample = seurat_patient@meta.data$sample,
                         Clone = seurat_patient@meta.data$Clone )

    library(stats)  # For dist function

  # Assess each clone existing in all samples at least 2 cells
  mat = table(pca_df$Clone, pca_df$Sample)
  Clones_kept = rownames(mat)[apply(mat,1,min) >= 2] 
  
#Functions
#calculate distance of random pair of cell 
random_pair_distance = function(iter, sample_data){
  num_cells <- nrow(sample_data)
  selected_indices <- sample(1:num_cells, 2, replace = F)
  dist_value <- dist(sample_data[selected_indices, ], method = "euclidean")
}

#mean of each clone and sample
mean_random_pair_distance_sample <- function(sample, pca_clone_df, n_iter = 1000) {
    sample_data = pca_clone_df[pca_clone_df$Sample == sample,1:30]
    distances = map(1:n_iter, ~random_pair_distance(.x,sample_data= sample_data))
    mean(unlist(distances))
  }

#calculate all samples for each clone
clone_result <- function(clone){
    # Assess each sample
    Samples <- unique(pca_df$Sample)
    pca_clone_df <- pca_df[pca_df$Clone == clone, ]
    samples_distances = map(Samples, ~mean_random_pair_distance_sample(.x, pca_clone_df = pca_clone_df, n_iter = 1000))
    distance_clone_df = data.frame(Sample = Samples, Distance = unlist(samples_distances))
    distance_clone_df$Clone = clone
    distance_clone_df
}

  Distances = map_dfr( Clones_kept, clone_result )
  rownames(Distances) = NULL
  Distances$Patient = patient
  Distances 
}




Patients_intra_clone_ITH = map_dfr(patients, intra_clone_ITH )
colnames(Patients_intra_clone_ITH) = c("Treatment","ITH","Clone", "Patient")
Patients_intra_clone_ITH = Patients_intra_clone_ITH %>% filter(Clone  != "UnI", Treatment !=  "Ex vivo")
Patients_intra_clone_ITH$Treatment = factor(Patients_intra_clone_ITH$Treatment,levels = c("Untreated","In vivo") )


#Plot for each patient 

#OO77
ITH_OO77 <- Patients_intra_clone_ITH %>% filter(Patient == "OO77")

clone_colors_77 <- c(`2` = "lightgreen", `5a` =  "lightblue",`5b` = "dodgerblue4",  `6` = "mediumseagreen")

#Plot
ITH_OO77  %>%  ggplot(aes(x= Treatment, y = ITH, group = Clone))+
  geom_point(aes(colour = Treatment), size = 2, position = position_dodge(width = 0)) +
  geom_line(aes( colour= Clone) ,linewidth = 1, alpha = 0.8, position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", clone_colors_77))+
  scale_y_continuous(limits = c(0,115))+
  ylab("Intra-clone heterogeneity")   +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))

  ggsave(file.path(plot_dir,paste0("OO77_Intra_clone_Transcriptomic_ITH_in_vivo.pdf")),width= 4.5,height= 3) 


#OO99
ITH_OO99 <- Patients_intra_clone_ITH %>% filter(Patient == "OO99")

clone_colors_99 <- clone_colors_99 <- c(  `3` = "moccasin" , `4` = "dodgerblue4", `6` ="yellow2",
            `7` = "lightseagreen", `8` = "lightsalmon4" , `9`= "plum" , 
             `10` = "tomato", `11` = "mediumseagreen", `12` = "purple")

#Plot
ITH_OO99  %>%  ggplot(aes(x= Treatment, y = ITH, group = Clone))+
  geom_point(aes(colour = Treatment), size = 2, position = position_dodge(width = 0)) +
  geom_line(aes( colour= Clone) ,linewidth = 1, alpha = 0.8, position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", clone_colors_99))+
  scale_y_continuous(limits = c(0,115))+
  ylab("Intra-clone heterogeneity") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))

  ggsave(file.path(plot_dir,"OO99_Intra_clone_Transcriptomic_ITH_in_vivo.pdf"),width= 4.5,height= 3) 


#OO100
ITH_OO100 <- Patients_intra_clone_ITH %>% filter(Patient == "OO100")

clone_colors_100 <- c( `2` = "burlywood2",
            `3` = "mediumseagreen",
            `4` = "red",
            `5` = "tan4", 
            `7` = "plum", 
            `8` = "dodgerblue4", 
            `9` = "purple",
            `10` = "yellow2",
            `11` = "lightpink",
            `12` = "darkorange")

#Plot
ITH_OO100  %>%  ggplot(aes(x= Treatment, y = ITH, group = Clone))+
  geom_point(aes(colour = Treatment), size = 2, position = position_dodge(width = 0)) +
  geom_line(aes( colour= Clone) ,linewidth = 1, alpha = 0.8, position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", clone_colors_100))+
  scale_y_continuous(limits = c(0,115))+
  ylab("Intra-clone heterogeneity")  +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))

  ggsave(file.path(plot_dir,"OO100_Intra_clone_Transcriptomic_ITH_in_vivo.pdf"),width= 4.5,height= 3) 