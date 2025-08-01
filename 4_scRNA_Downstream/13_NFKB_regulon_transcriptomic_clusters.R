
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


meta <- seurat_merged@meta.data
meta <- meta %>%
  mutate(sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ),
   treatment = sub("^[^ ]* ","",sample),
   patient = sub(" .*","",sample)
   )


colnames(meta)[14:21] = Pathway_names 
seurat_merged@meta.data <- meta

seurat_merged@meta.data$CB = rownames(seurat_merged@meta.data )

#Seurat object with transcriptomic cluster identity
seurat_integrated <- readRDS(file.path(rds_3_path, "clustered_0.35_integrated_seurat_ref.rds"))

seurat_integrated$CB = rownames(seurat_integrated@meta.data)

seurat_integrated@meta.data  = seurat_integrated@meta.data %>% left_join(seurat_merged@meta.data[, c("CB","patient", "treatment", "TNFA_SIGNALING_VIA_NFKB")], by = "CB")




#load regulon AUC score 
auc <- read.delim("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCENIC/results/auc_mtx_NFKB.tsv")
colnames(auc) = paste0(colnames(auc), "_regulon")
seurat_integrated@meta.data = cbind(seurat_integrated@meta.data, auc)
rownames(seurat_integrated@meta.data ) = seurat_integrated@meta.data$CB
# FeaturePlot(seurat_integrated, features = c( "RELB_regulon", "NFKB1_regulon","NFKB2_regulon", "REL_regulon"), raster=T, pt.size = 1.5) 
# ggsave(file.path(plot_result_path, "FeaturePlot_regulons.pdf"), width = 15, height = 10)

VlnPlot(seurat_integrated,group.by= "seurat_clusters", features = c( "RELB_regulon", "NFKB1_regulon","NFKB2_regulon", "REL_regulon"),ncol=2, alpha = 0.01) 
# ggsave(file.path(plot_result_path, "VlnPlot_regulons.pdf"), width = 10, height = 5)


FeaturePlot(seurat_integrated, features =  "NFKB1_regulon", raster=T, pt.size = 1.5)
ggsave(file.path(plot_result_path, "FeaturePlot_NFKB1_regulons.pdf"), width = 12, height = 7)

df = seurat_integrated@meta.data  

#Plot Violin plot for NFKB1 regulon activities



custom_palette <- c(  "#80B1D3", "limegreen", "#FCCDE5" ,
                  "#FDB462","#D9D9D9", "#CCEBC5" ,
                 "#BC80BD" , "khaki1",  "brown1") 

df$seurat_clusters = factor(df$seurat_clusters, levels = 1:9)

plot <- df %>% ggplot(aes(x= seurat_clusters, y = NFKB1_regulon, color = seurat_clusters  )) +  
      geom_violin(linewidth = 1)+
      geom_boxplot(width = 0.2, linewidth = 1) +
      scale_color_manual(values = custom_palette ) +
      ylab("NFKB1 regulon activity")+
      theme_classic()+
      theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 31, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 31, color = "black"),
    axis.title.y = element_text(size = 31, color = "black"),
    axis.text.y = element_text(size = 31, color = "black"),
  ) + guides(color = guide_legend(nrow = 1))

plot
ggsave(file.path(plot_result_path , "NFKB1_regulon_across_seurat_clusters.pdf"),width= 8,height= 4.8,scale = 1)

#Save Figure data
df <- df %>%
  mutate(
    Patient = word(Sample, 1),
    Treatment = word(Sample, 2)
  ) %>% select(Patient, Treatment, seurat_clusters, NFKB1_regulon)

write.table(df, file.path(Figures_data, "Figure_6C_NFKB1_regulon_activity_across_clusters"), sep = "\t", quote = F)


#Regulon across Sample
df$Treatment = factor (df$Treatment , levels =  c("Untreated", "Ex" ,"In" ))
df$Patient = factor(df$Patient, levels = c("OO77", "OO99", "OO100"))


plot <- df %>% ggplot(aes(x= Treatment, y = NFKB1_regulon, color = Treatment)) +  
      geom_violin(linewidth = 1)+
      geom_boxplot(width = 0.2, linewidth = 1) +
      scale_color_manual(values = c("deepskyblue",  "tomato", "plum"),labels= c("Untreated" = "Untreated","Ex"= "Ex vivo", "In"="In vivo"))+
      ylab("NFKB1 regulon activity")+
      theme_classic()+
      theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
  ) + guides(color = guide_legend(nrow = 1)) + facet_wrap( vars(Patient), nrow= 1)

plot
ggsave(file.path(plot_result_path , "NFKB1_regulon_across_sample.pdf"),width= 12,height= 4.8,scale = 1)






#OO77

df_OO77 = df %>% as.data.frame() %>% filter(patient == "OO77", treatment %in% c("Untreated", "In vivo")) 
OO77_summary = df_OO77 %>% group_by(seurat_clusters, treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB), 
                                                                              mean_NFKB1_regulon = mean(NFKB1_regulon))
OO77_summary$treatment = factor(OO77_summary$treatment , levels = c("Untreated", "In vivo"))
OO77_summary$seurat_clusters = factor(OO77_summary$seurat_clusters, levels = 1:9)

names(custom_palette ) = 1:9


#pathway
OO77_summary %>% 
  ggplot(aes(x = treatment, y = mean_pathway, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_NFKb_each_clusters_OO77_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)

#regulon
OO77_summary %>% 
  ggplot(aes(x = treatment, y = mean_NFKB1_regulon, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
  ylab("Mean NFKB1 regulon AUC") +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_AUC_NFKB1_regulon_each_clusters_OO77_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)



#OO99

df_OO99 = df %>% as.data.frame() %>% filter(patient == "OO99", treatment %in% c("Untreated", "In vivo")) 
OO99_summary = df_OO99 %>% group_by(seurat_clusters, treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB),
                                                                                             mean_NFKB1_regulon = mean(NFKB1_regulon))
OO99_summary$treatment = factor(OO99_summary$treatment , levels = c("Untreated", "In vivo"))
OO99_summary$seurat_clusters = factor(OO99_summary$seurat_clusters, levels = 1:9)

names(custom_palette ) = 1:9
OO99_summary %>% 
  ggplot(aes(x = treatment, y = mean_pathway, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_NFKb_each_clusters_OO99_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)


#save Figure data
data = df %>% select(CB, patient, treatment,seurat_clusters ,TNFA_SIGNALING_VIA_NFKB)
write.table(data, file.path(Figures_data, "Figure_5D_NFKB_pathway_score_intra_transcriptional_cluster.txt"), sep = "\t", quote = F)


#regulon
OO99_summary %>% 
  ggplot(aes(x = treatment, y = mean_NFKB1_regulon, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
  ylab("Mean NFKB1 regulon AUC") +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_AUC_NFKB1_regulon_each_clusters_OO99_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)


#OO100

df_OO100 = df %>% as.data.frame() %>% filter(patient == "OO100", treatment %in% c("Untreated", "In vivo")) 
OO100_summary = df_OO100 %>% group_by(seurat_clusters, treatment) %>% summarize(mean_pathway = mean( TNFA_SIGNALING_VIA_NFKB),
                                                                                  mean_NFKB1_regulon = mean(NFKB1_regulon))
OO100_summary$treatment = factor(OO100_summary$treatment , levels = c("Untreated", "In vivo"))
OO100_summary$seurat_clusters = factor(OO100_summary$seurat_clusters, levels = 1:9)

names(custom_palette ) = 1:9
OO100_summary %>% 
  ggplot(aes(x = treatment, y = mean_pathway, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_NFKb_each_clusters_OO100_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)



#regulon
OO100_summary %>% 
  ggplot(aes(x = treatment, y = mean_NFKB1_regulon, group = seurat_clusters)) +
  geom_point(aes(colour = treatment), 
             size = 2, 
             position = position_dodge(width = 0)) + 
  geom_line(aes(group = seurat_clusters, colour = seurat_clusters), 
            linewidth = 1.5, 
            alpha = 1, 
            position = position_dodge(width = 0)) +
  scale_color_manual(values = c("Untreated" = "deepskyblue2", 
                                "In vivo" = "plum", 
                                custom_palette )) +
 # scale_y_continuous(limits = c(-0.05, 0.15)) +
  ylab("Mean NFKB1 regulon AUC") +
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

ggsave(file.path(plot_result_path , "LinePlot_Mean_AUC_NFKB1_regulon_each_clusters_OO100_Untreated_In_Vivo.pdf"),width= 700,height= 500,dpi=150,units="px",scale = 1)


