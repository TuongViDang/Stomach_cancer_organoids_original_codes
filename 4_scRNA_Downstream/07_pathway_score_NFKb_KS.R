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
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(ggpubr)



#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
txt_files<-file.path(main_path,"txt_files")
rds_path<-file.path(main_path,'rds')
rds_path_3 <-file.path(main_path,'rds',"three_patients")
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_3_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_3_path,"pathway_score")
dir.create(plot_result_path)

#Figures data path
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)


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


DefaultAssay(seurat_merged)<-'RNA'
options(future.globals.maxSize = 17 * 1024^3)
seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, npcs = 30, )
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:30, verbose = F)

seurat_merged$Patient<- str_extract(seurat_merged$orig.ident,"^[^[_-]]*")
seurat_merged$Treatment<-str_extract(seurat_merged$orig.ident,"[^-]*$")

meta <- seurat_merged@meta.data
meta <- meta %>%
  mutate(Sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ),
   Treatment = sub("^[^ ]* ","",Sample))


colnames(meta)[14:21] = Pathway_names 
seurat_merged@meta.data <- meta

#Plot cummulative curve for NFkb
names(Pathway_names) <- NULL
data = meta %>% select(Sample, Patient, Treatment, all_of(Pathway_names ))
data_long = data %>% pivot_longer(cols = 4:11, names_to = "pathway",values_to = "score" )


data_long = data_long %>% filter(pathway == "TNFA_SIGNALING_VIA_NFKB")
data_long$Treatment = factor(data_long$Treatment, levels = c( "In vivo","Untreated", "Ex vivo"))
data_long$Patient = factor(data_long$Patient, levels = c("OO77", "OO99", "OO100"))

# Manual tatistics

# Prepare data: for each patient and treatment, sort scores and compute empirical CDF (cumulative proportion)
df_cdf <- data_long %>%
  group_by(Patient, Treatment) %>%
  arrange(score) %>%
  mutate(cum_prob = ecdf(score)(score)) %>%
  ungroup()

colors = c("plum","deepskyblue2","tomato")
# Plot cumulative distribution lines
ggplot(df_cdf, aes(x = score, y = cum_prob, color = Treatment, group =  Treatment)) +
  geom_line(alpha = 1, linewidth= 0.7) +
  scale_color_manual(values = colors) +
  labs(title = "Pathway score TNFA_SIGNALING_VIA_NFKB",
       x = "Score",
       y = "Cumulative Proportion") +
  theme_classic () +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, color = "black" ),
        legend.text = element_text(size = 15, color = "black" ),
        plot.title =  element_text(size = 15, color = "black" ),
        axis.title.y.left = element_text(size = 15, color = "black" ),
        axis.text.y.left = element_text(size = 15,  color = "black"),
        axis.text.x.bottom = element_text(size = 15,  color = "black"),
        axis.title.x.bottom = element_text(size = 15,  color = "black"),
        strip.text = element_text(size = 15, color = "black")  )+
  facet_wrap(vars(Patient))
ggsave(file.path(plot_result_path, "CDF_Pathways_score.pdf"), width = 12, height = 3.5, scale = 1)

#save Figure data
write.table(df_cdf , file.path(Figures_data, "Figure_2C_NFKB_pathway_score_Cummulative_Curve.txt"), sep = "\t", quote = F)


#Statistics
ks_results <- data_long %>%
  group_by(patient) %>%
  summarise(
    ks_exvivo = {
      untreated_scores <- score[Treatment == "Untreated"]
      exvivo_scores <- score[Treatment == "Ex vivo"]
      if (length(untreated_scores) > 1 & length(exvivo_scores) > 1) {
        ks.test( exvivo_scores,untreated_scores, alternative = "greater")$p.value
      } else {
        NA
      }
    },
    ks_invivo = {
      untreated_scores <- score[Treatment == "Untreated"]
      invivo_scores <- score[Treatment == "In vivo"]
      if (length(untreated_scores) > 1 & length(invivo_scores) > 1) {
        ks.test(invivo_scores, untreated_scores, alternative = "greater")$p.value
      } else {
        NA
      }
    }
  )

# Flatten p-values into a single vector
pvals_vector <- unlist(ks_results[, c("ks_exvivo", "ks_invivo")])

# Apply multiple testing correction
pvals_adj <- p.adjust(pvals_vector, method = "BH")  # You can choose other methods like "bonferroni"

# Reshape adjusted p-values back into a data frame
ks_results$ks_exvivo_adj <- pvals_adj[1:3]
ks_results$ks_invivo_adj <- pvals_adj[4:6]


#### Get genetic clone labelling

OO77_rds_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA","OO77","tree30_change_name","rds")
OO77_no = readRDS(file.path(OO77_rds_dir, paste0("OO77","-no_Clones_merged_separated_c5_filtered.RDS" )))
OO77_vitro = readRDS(file.path(OO77_rds_dir, paste0("OO77","-vitro_Clones_merged_separated_c5_filtered.RDS")))
OO77_vivo =  readRDS(file.path(OO77_rds_dir, paste0("OO77","-vivo_Clones_merged_separated_c5_filtered.RDS" )))

OO99_rds_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA","OO99","tree30_change_name","rds")
OO99_no = readRDS(file.path(OO99_rds_dir, paste0("OO99","-no_Clones_filtered.RDS" )))
OO99_vitro = readRDS(file.path(OO99_rds_dir, paste0("OO99","-vitro_Clones_filtered.RDS")))
OO99_vivo =  readRDS(file.path(OO99_rds_dir, paste0("OO99","-vivo_Clones_filtered.RDS" )))

OO100_rds_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA","OO100","tree30_change_name","rds")
OO100_no = readRDS(file.path(OO100_rds_dir, paste0("OO100","-no_Clones_filtered.RDS" )))
OO100_vitro = readRDS(file.path(OO100_rds_dir, paste0("OO100","-vitro_Clones_filtered.RDS")))
OO100_vivo =  readRDS(file.path(OO100_rds_dir, paste0("OO100","-vivo_Clones_filtered.RDS" )))

meta_77 = rbind(OO77_no@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")], 
                 OO77_vitro@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")]) %>% 
            rbind(OO77_vivo@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")])

meta_99 = rbind(OO99_no@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")], 
                 OO99_vitro@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")]) %>% 
            rbind(OO99_vivo@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")])

meta_100 = rbind(OO100_no@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")], 
                 OO100_vitro@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")]) %>% 
            rbind(OO100_vivo@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")])


##### Compare  Untreated and In vivo of the same clone
##OO77
seurats_OO77 = subset( seurat_merged, subset = orig.ident %in% c("OO77-no","OO77-vivo"))
meta_77_no_vivo = seurats_OO77@meta.data 
meta_77_no_vivo$CB_seurat = rownames(meta_77_no_vivo )
meta_77_no_vivo = meta_77_no_vivo %>% left_join(meta_77[,c("CB_seurat", "Clone")], by = "CB_seurat")
seurats_OO77@meta.data = meta_77_no_vivo
rownames(seurats_OO77@meta.data) = seurats_OO77@meta.data$CB_seurat
table(seurats_OO77$Sample,seurats_OO77$Clone) #remove 5b,7,8 not at least 3 cells in both


meta_77_no_vivo =  meta_77_no_vivo %>% select(CB_seurat,Patient, Treatment, Clone, TNFA_SIGNALING_VIA_NFKB)
meta_77_no_vivo = meta_77_no_vivo   %>% filter(Clone %in% c("2","5a","6") )

OO77_summary = meta_77_no_vivo %>% group_by(Clone, Treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB))
OO77_summary$Treatment = factor(OO77_summary$Treatment , levels = c("Untreated", "In vivo"))
OO77_summary$Clone = factor(OO77_summary$Clone , levels  = c("2","5a","6"))

clone_colors_77 <- c(`2` = "lightgreen", `5a` =  "lightblue",`5b` = "dodgerblue4", 
  `6` = "mediumseagreen",`7` = "plum", `8` = "lightgoldenrod1", `UnI` = "gray92"  # Color for NA cells
)
OO77_summary %>% 
  ggplot(aes(x = Treatment, y = mean_pathway, group = Clone)) +
  geom_point(aes(colour = Treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = Clone, colour = Clone), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                clone_colors_77)) +
  scale_y_continuous(limits = c(-0.05, 0.1)) +
  ylab("Mean Pathway Score") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    panel.spacing = unit(0.5, "lines")
  )

ggsave(file.path(plot_result_path , "LinePlot_Mean_NFKb_each_Clone_OO77_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)

#save Figure data
write.table(meta_77_no_vivo , file.path(Figures_data, "Figure_5A_NFKB_pathway_score_intra_genetic_clone_OO77.txt"), sep = "\t", quote = F)

##OO99
seurats_OO99 = subset( seurat_merged, subset = orig.ident %in% c("OO99-no","OO99-vivo"))
meta_99_no_vivo = seurats_OO99@meta.data 
meta_99_no_vivo$CB_seurat = rownames(meta_99_no_vivo )
meta_99_no_vivo = meta_99_no_vivo %>% left_join(meta_99[,c("CB_seurat", "Clone")], by = "CB_seurat")
seurats_OO99@meta.data = meta_99_no_vivo
rownames(seurats_OO99@meta.data) = seurats_OO99@meta.data$CB_seurat
table(seurats_OO99$Sample,seurats_OO99$Clone) #remove 2,7,9,12 not at least 3 cells in both


meta_99_no_vivo =  meta_99_no_vivo %>% select(CB_seurat,Patient, Treatment, Clone, TNFA_SIGNALING_VIA_NFKB)

meta_99_no_vivo = meta_99_no_vivo %>% filter(Clone %in% c("3","4","6","8","10","11"))

OO99_summary = meta_99_no_vivo %>% group_by(Clone, Treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB))
OO99_summary$treatment = factor(OO99_summary$Treatment , levels = c("Untreated", "In vivo"))
OO99_summary$Clone = factor(OO99_summary$Clone , levels  = c("3","4","6","8","10","11"))

clone_colors_99 <- c( `2`= "aquamarine2", `3` = "moccasin" , `4` = "dodgerblue4", `6` ="yellow2",
            `7` = "lightseagreen", `8` = "lightsalmon4" , `9`= "plum" , 
             `10` = "tomato", `11` = "mediumseagreen", `12` = "purple", `UnI` = "gray92"  # Color for NA cells
)



OO99_summary  %>% 
  ggplot(aes(x = Treatment, y = mean_pathway, group = Clone)) +
  geom_point(aes(colour = Treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = Clone, colour = Clone), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                clone_colors_99)) +
  scale_y_continuous(limits = c(-0.01, 0.12)) +
  ylab("Mean Pathway Score") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    panel.spacing = unit(0.5, "lines")
  )


ggsave(file.path(plot_result_path, "LinePlot_Mean_NFKb_each_Clone_OO99_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)


#save Figure data
write.table(meta_99_no_vivo , file.path(Figures_data, "Figure_5A_NFKB_pathway_score_intra_genetic_clone_OO99.txt"), sep = "\t", quote = F)


##OO100
seurats_OO100 = subset( seurat_merged, subset = orig.ident %in% c("OO100-no","OO100-vivo"))
meta_100_no_vivo = seurats_OO100@meta.data 
meta_100_no_vivo$CB_seurat = rownames(meta_100_no_vivo )
meta_100_no_vivo = meta_100_no_vivo %>% left_join(meta_100[,c("CB_seurat", "Clone")], by = "CB_seurat")
seurats_OO100@meta.data = meta_100_no_vivo
rownames(seurats_OO100@meta.data) = seurats_OO100@meta.data$CB_seurat
table(seurats_OO100$Sample,seurats_OO100$Clone) #remove 2,5,7,9, 10,11,12 not at least 3 cells in both


meta_100_no_vivo =  meta_100_no_vivo %>% select(CB_seurat,Patient, Treatment, Clone, TNFA_SIGNALING_VIA_NFKB)
meta_100_no_vivo = meta_100_no_vivo  %>% filter(Clone %in% c("3","4","8"))

OO100_summary = meta_100_no_vivo %>%  group_by(Clone, Treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB))
OO100_summary$Treatment = factor(OO100_summary$Treatment , levels = c("Untreated", "In vivo"))
OO100_summary$Clone = factor(OO100_summary$Clone , levels  = c("3","4","8"))

clone_colors_100 <- c( `2` = "burlywood2",
            `3` = "mediumseagreen",
            `4` = "red",
            `5` = "tan4", 
            `7` = "plum", 
            `8` = "dodgerblue4", 
            `9` = "purple",
            `10` = "yellow2",
            `11` = "lightpink",
            `12` = "darkorange",
            `UnI` = "gray92"  # Color for NA cells
)



OO100_summary %>%  
  ggplot(aes(x = Treatment, y = mean_pathway, group = Clone)) +
  geom_point(aes(colour = Treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = Clone, colour = Clone), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                clone_colors_100)) +
  scale_y_continuous(limits = c(-0.01, 0.12)) +
  ylab("Mean Pathway Score") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    panel.spacing = unit(0.5, "lines")
  )

ggsave(file.path(plot_result_path, "LinePlot_Mean_NFKb_each_Clone_OO100_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)

#save Figure data
write.table(meta_100_no_vivo , file.path(Figures_data, "Figure_5A_NFKB_pathway_score_intra_genetic_clone_OO100.txt"), sep = "\t", quote = F)



#### Plot heatmap of intra clone common gene in pseudobulk ###
genes_common = c( "DUSP1"   , "SQSTM1" ,  "EDN1" ,   "PPP1R15A" ,"PTGS2"   ,
                   "EGR1" ,   "HBEGF" ,   "CXCL1"  ,  "SOCS3" ,   "GADD45B" ,
                    "IER5"   ,  "BTG2"  , "EGR2"  ,   "BCL6" )


#OO77
DefaultAssay(seurats_OO77) = "RNA"
options(future.globals.maxSize = 6 * 1024^3)
seurats_OO77 = SCTransform(seurats_OO77,vst.flavor='v1',
                          variable.features.n = nrow(seurats_OO77))

expression_mat_OO77 = seurats_OO77$SCT$scale.data
expression_intersect_OO77 = expression_mat_OO77[genes_common,rownames(seurats_OO77@meta.data )]

expression_intersect_OO77_t = t(expression_intersect_OO77) %>% as.data.frame() 

expression_intersect_OO77_t  = expression_intersect_OO77_t[meta_77_no_vivo$CB_seurat,]
expression_intersect_OO77_t$CB_seurat = rownames(expression_intersect_OO77_t)
expression_intersect_OO77_t = expression_intersect_OO77_t %>% left_join(meta_77_no_vivo, by = "CB_seurat")

df_long_OO77 = expression_intersect_OO77_t %>% pivot_longer(cols = 1:14, names_to = "gene", values_to = "expression" )
Gene_summary_OO77 = Combined_df_long %>% group_by(gene, Clone, Treatment) %>% summarize( mean_expression = mean(expression))

log2fc_OO77 <- Gene_summary_OO77  %>%
  pivot_wider(names_from = Treatment, values_from = mean_expression) %>%
  mutate(logFC_ln = `In vivo` - Untreated,
         log2FC = logFC_ln / log(2))


log2fc_OO77 <- log2fc_OO77 %>% filter( Clone %in% c("2","5a","6"))
log2fc_OO77$Clone <- factor(log2fc_OO77$Clone  , levels  = c("2","5a","6"))
log2fc_OO77$patient = "OO77"


#OO99
DefaultAssay(seurats_OO99) = "RNA"
options(future.globals.maxSize = 6 * 1024^3)
seurats_OO99 = SCTransform(seurats_OO99,vst.flavor='v1',
                          variable.features.n = nrow(seurats_OO99))

expression_mat_OO99 = seurats_OO99$SCT$scale.data
expression_intersect_OO99 = expression_mat_OO99[genes_common,rownames(seurats_OO99@meta.data )]

expression_intersect_OO99_t = t(expression_intersect_OO99) %>% as.data.frame() 
expression_intersect_OO99_t  = expression_intersect_OO99_t[meta_99_no_vivo$CB_seurat,]
expression_intersect_OO99_t$CB_seurat = rownames(expression_intersect_OO99_t)
expression_intersect_OO99_t = expression_intersect_OO99_t %>% left_join(meta_99_no_vivo, by = "CB_seurat")


df_long_OO99 = expression_intersect_OO99_t %>% pivot_longer(cols = 1:14, names_to = "gene", values_to = "expression" )
Gene_summary_OO99 = Combined_df_long %>% group_by(gene, Clone, Treatment) %>% summarize( mean_expression = mean(expression))


log2fc_OO99 <- Gene_summary_OO99  %>%
  pivot_wider(names_from = Treatment, values_from = mean_expression) %>%
  mutate(logFC_ln = `In vivo` - Untreated,
         log2FC = logFC_ln / log(2))

log2fc_OO99 <- log2fc_OO99 %>% filter( Clone %in%  c("3","4","6","8","10","11"))
log2fc_OO99$Clone <- factor(log2fc_OO99$Clone  , levels  =  c("3","4","6","8","10","11") )
log2fc_OO99$patient = "OO99"

#OO100
DefaultAssay(seurats_OO100) = "RNA"
options(future.globals.maxSize = 6 * 1024^3)
seurats_OO100 = SCTransform(seurats_OO100,vst.flavor='v1',
                          variable.features.n = nrow(seurats_OO100))

expression_mat_OO100 = seurats_OO100$SCT$scale.data
expression_intersect_OO100 = expression_mat_OO100[genes_common,rownames(seurats_OO100@meta.data )]

expression_intersect_OO100_t = t(expression_intersect_OO100) %>% as.data.frame() 
expression_intersect_OO100_t  = expression_intersect_OO100_t[meta_100_no_vivo$CB_seurat,]
expression_intersect_OO100_t$CB_seurat = rownames(expression_intersect_OO100_t)
expression_intersect_OO100_t = expression_intersect_OO100_t %>% left_join(meta_100_no_vivo, by = "CB_seurat")

df_long_OO100 = expression_intersect_OO100_t %>% pivot_longer(cols = 1:14, names_to = "gene", values_to = "expression" )
Gene_summary_OO100 = Combined_df_long %>% group_by(gene, Clone, Treatment) %>% summarize( mean_expression = mean(expression))

log2fc_OO100 <- Gene_summary_OO100  %>%
  pivot_wider(names_from = Treatment, values_from = mean_expression) %>%
  mutate(logFC_ln = `In vivo` - Untreated,
         log2FC = logFC_ln / log(2))

log2fc_OO100 <- log2fc_OO100 %>% filter( Clone %in%  c("3","4","8"))
log2fc_OO100$Clone <- factor(log2fc_OO100$Clone  , levels  =  c("3","4","8"))
log2fc_OO100$patient = "OO100"


#Combine 
Combined_df = rbind(log2fc_OO77, log2fc_OO99, log2fc_OO100) %>% select(-Untreated,- `In vivo` ,-logFC_ln)
Combined_df_wide = Combined_df %>% pivot_wider(id_cols= "gene", names_from = c("patient","Clone"), values_from ="log2FC" ,names_sep = "_")
matrix = Combined_df_wide[,2:13] %>% as.data.frame()
rownames(matrix) = Combined_df_wide$gene


Clone_colors <- c(
  `1` = "lightgreen", `2` = "lightblue",  `3` = "mediumseagreen", 
  `4` = "moccasin", `5` = "dodgerblue4",`6` = "yellow2", `7` = "lightsalmon4", `8` = "tomato", `9` = "mediumseagreen", 
  `10` = "mediumseagreen",`11` = "red", `12` = "dodgerblue4"
)

patient_colors <- c("OO77" = "goldenrod3", "OO99" = "darkcyan", "OO100" = "indianred1")

anno_colors <- list(Clone = Clone_colors, 
                   patient = patient_colors )

colnames(matrix) = names(Clone_colors)
anno_col = data.frame(Clone =  names(Clone_colors), 
                      patient = c(rep("OO77",3 ), rep("OO99",
                      6), rep("OO100",3)))
matrix = matrix[genes_common,]

pheatmap::pheatmap(as.matrix(matrix),  
                   color = colorRampPalette(brewer.pal(n = 11, name = "RdBu")[1:7])(100),
                   show_rownames = TRUE, 
                   show_colnames = FALSE,   
                   cluster_row = FALSE, 
                   cluster_col = FALSE, 
                   annotation_col = anno_col,
                  annotation_colors = anno_colors,
                   gaps_col = c( 3 , 9 ),
                   filename = file.path( plot_result_path , "Intra_clone_Heatmap_gene_common_NFkb.pdf"),
                   width = 7, 
                   height = 3.7)