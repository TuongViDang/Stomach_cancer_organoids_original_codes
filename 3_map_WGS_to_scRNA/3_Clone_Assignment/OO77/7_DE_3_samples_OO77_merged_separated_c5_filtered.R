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

patient="OO77"


Integrate_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA",patient,"tree30_change_name")
plot_dir = file.path(Integrate_dir,"plots","DE_3_samples_merged_separated_c5_filtered")
dir.create(plot_dir, recursive = T)
txt_files = file.path(Integrate_dir,"txt_files")
dir.create(txt_files,recursive=T)
rds_dir = file.path(Integrate_dir,"rds")
dir.create(rds_dir)

#Figures data path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)

seurat_patient = readRDS(file.path(rds_dir,paste0(patient,"three_samples_final_Clones.rds")))
table(seurat_patient$Clone)    
seurat_patient$Clone = as.character(seurat_patient$Clone)
#Combine 6,7,8
seurat_patient@meta.data = seurat_patient@meta.data %>% mutate(Clone_new = ifelse(Clone %in% c("7","8"), "6",Clone))
                                                        

seurat_patient_sub  = subset(seurat_patient , subset = Clone_new %in% c("2","5a","5b","6") )     
table(seurat_patient_sub$Clone_new)


seurat_patient_sub_2 =  subset(seurat_patient , subset = Clone_new %in% c("2","5a","5b") )     

seurat_patient_sub_6 =  subset(seurat_patient , subset = Clone_new %in% c("6","5a","5b") )     

#Differential analysis
#DE 
Idents(seurat_patient_sub_2) = "Clone_new"
Markers_2 <- FindAllMarkers(seurat_patient_sub_2 , verbose = TRUE, only.pos = FALSE)


Idents(seurat_patient_sub_6) = "Clone_new"
Markers_6 <- FindAllMarkers(seurat_patient_sub_6 , verbose = TRUE, only.pos = FALSE)


#save markers
write.table(Markers_2 , file.path(txt_files,paste0(patient,"_clone_markers_2_5a_5b.txt")), sep = "\t", quote = F)
write.table(Markers_6 , file.path(txt_files,paste0(patient,"_clone_markers_6_5a_5b.txt")), sep = "\t", quote = F)


#reload data
Markers_2 <- read.delim(file.path(txt_files,paste0(patient,"_clone_markers_2_5a_5b.txt")))
Markers_6 <- read.delim(file.path(txt_files,paste0(patient,"_clone_markers_6_5a_5b.txt")))

#Choose avg_log2FC >= 1, pct.1 >= 0.7
Markers_sig_2 = Markers_2 %>% filter(p_val_adj < 0.05, avg_log2FC >= 0.5, pct.1 >= 0.7)
table(Markers_sig_2$cluster)
 

Markers_sig_6 = Markers_6 %>% filter(p_val_adj < 0.05, avg_log2FC >= 0.5, pct.1 >= 0.7)
table(Markers_sig_6$cluster)


#GO ANALYSIS
#Prepare for GO enrichment analysis
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(ggpubr)

columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
All_genes<- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(seurat_patient),keytype='SYMBOL',column=c('ENSEMBL'),muliVals=first)
All_genes<-All_genes[!(is.na(All_genes))]


H<-msigdbr(species = 'Homo sapiens',category='H')
H.ensemble<-H%>%dplyr::select(gs_name,ensembl_gene)


Find_GO_seurat <- function(seurat, Marker_sig){
  Find_GO <-function(seurat, gene_list){
    markers_genes <-AnnotationDbi::mapIds(org.Hs.eg.db,keys= gene_list,keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
    markers_genes<-markers_genes[!(is.na(markers_genes))]
    enrich.cluster <- enricher(markers_genes,universe = All_genes,TERM2GENE = H.ensemble)
     enrich.cluster_df <- data.frame()
     if (!(is.null(enrich.cluster ))){
        enrich.cluster_df <- enrich.cluster@result %>% cbind(Clone = clone,.) %>% filter( p.adjust < 0.05)}
        enrich.cluster_df 
} 
  GO_df = data.frame()
  for (clone in unique(seurat$Clone_new)){
    gene_list_clone = Marker_sig %>% filter(cluster == clone) %>% pull(gene)
    GO_df = rbind(GO_df, Find_GO(seurat, gene_list_clone))
}
GO_df
}
#For seurat_patient_sub_2

Clone_2_5a_5b_GO = Find_GO_seurat(seurat_patient_sub_2,Markers_sig_2)


Clone_2_5a_5b_GO$ID = sub("HALLMARK_", "", Clone_2_5a_5b_GO$ID)

#Dot plot

# Convert GeneRatio to numeric values
Clone_2_5a_5b_GO <- Clone_2_5a_5b_GO %>%
  mutate(GeneRatio_numeric = sapply(GeneRatio, function(x) {
    ratio <- strsplit(x, "/")[[1]]
    as.numeric(ratio[1]) / as.numeric(ratio[2])
  })) %>% mutate(neg_log_p_adj = - log10( p.adjust))



pathways_2 <- rev(c(
  "TNFA_SIGNALING_VIA_NFKB", "P53_PATHWAY",
   "HYPOXIA", "APOPTOSIS","TGF_BETA_SIGNALING", "ANDROGEN_RESPONSE"  ,   "FATTY_ACID_METABOLISM",
    "OXIDATIVE_PHOSPHORYLATION","ADIPOGENESIS", "E2F_TARGETS",
  "G2M_CHECKPOINT", "MITOTIC_SPINDLE", "MYC_TARGETS_V1", "MTORC1_SIGNALING"
))

Clone_2_5a_5b_GO$ID = factor(Clone_2_5a_5b_GO$ID , levels = pathways_2)

plot <- ggplot(Clone_2_5a_5b_GO , aes(x = Clone, y = ID, fill = neg_log_p_adj , size = GeneRatio_numeric)) +
  geom_point(shape = 21, color = "black") + # shape = 21 allows fill and border colors
  scale_fill_viridis_c( limits = c(0,12),breaks = seq(0,10,4)) + 
  scale_size_continuous(limits = c(0, 0.45), breaks = seq(0.05,0.45,0.15)) +  # Set fixed range for size
  labs(
    x = "",
    y = "HALLMARK gene sets",
    fill = "-log10(adjusted p)",
    size = "Gene Ratio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"), # Adjust text size for clarity
    axis.text.x = element_text(size = 16, color = "black"),
  )
plot
ggsave(file.path(plot_dir, "GO_2_5a_5b_DotPlot.pdf"), width = 7.5, height = 5, scale = 1)

#save Figures data
data <- plot$data

write.table(data, file.path(Figures_data,"Figure_4D_GO_Clone_OO77_2_vs_5a_5b.txt" ), sep = "\t", quote = F)

#For seurat_patient_sub_6

Clone_6_5a_5b_GO = Find_GO_seurat(seurat_patient_sub_6,Markers_sig_6)


Clone_6_5a_5b_GO$ID = sub("HALLMARK_", "", Clone_6_5a_5b_GO$ID)
Clone_6_5a_5b_GO$Clone = factor(Clone_6_5a_5b_GO$Clone , levels = c("6","5a","5b"))
#Dot plot

# Convert GeneRatio to numeric values
Clone_6_5a_5b_GO <- Clone_6_5a_5b_GO %>%
  mutate(GeneRatio_numeric = sapply(GeneRatio, function(x) {
    ratio <- strsplit(x, "/")[[1]]
    as.numeric(ratio[1]) / as.numeric(ratio[2])
  })) %>% mutate(neg_log_p_adj = - log10( p.adjust))



pathways_6 <- rev(c(
  "TNFA_SIGNALING_VIA_NFKB", "P53_PATHWAY","XENOBIOTIC_METABOLISM",
   "OXIDATIVE_PHOSPHORYLATION", "MYC_TARGETS_V1",  "FATTY_ACID_METABOLISM",  "ESTROGEN_RESPONSE_LATE"
))

Clone_6_5a_5b_GO$ID = factor(Clone_6_5a_5b_GO$ID , levels = pathways_6)

plot <- ggplot(Clone_6_5a_5b_GO , aes(x = Clone, y = ID, fill = neg_log_p_adj , size = GeneRatio_numeric)) +
  geom_point(shape = 21, color = "black") + # shape = 21 allows fill and border colors
  scale_fill_viridis_c(limits = c(0,12), breaks = seq(0,12,4)) + 
  scale_size_continuous(limits = c(0, 0.45), breaks = seq(0.05,0.45,0.15)) +  # Set fixed range for size
  labs(
    x = "",
    y = "HALLMARK pathways",
    fill = "-log10(adjusted p)",
    size = "Gene Ratio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"), # Adjust text size for clarity
    axis.text.x = element_text(size = 16, color = "black")
  ) + 
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5, nrow = 1))
plot
ggsave(file.path(plot_dir, "GO_6_5a_5b_DotPlot.pdf"), width = 7.5, height = 2.8, scale = 1)

#save Figures data
data <- plot$data

write.table(data, file.path(Figures_data,"Figure_4D_GO_Clone_OO77_6_vs_5a_5b.txt" ), sep = "\t", quote = F)
