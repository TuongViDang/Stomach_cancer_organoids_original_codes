## Load libraries
library(Seurat)
library(tidyverse)
library(HGNChelper)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(RColorBrewer)


#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_path,"cluster_markers")


#Figures data path
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)


#Load file
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))



Markers_all_clusters <-read.csv(file.path(txt_files_3,paste0("markers_0.35_ref.txt")),sep = '\t')


library(dplyr)
Markers_all_clusters  <- Markers_all_clusters  %>%
  rowwise() %>%
  mutate(mean_log2FC = mean(c_across(contains("log2FC")), na.rm = TRUE),
         mean_pct.1 = mean(c_across(contains("pct.1")), na.rm = TRUE)) %>%
  ungroup()

top_markers <- Markers_all_clusters %>% 
                group_by(cluster) %>% filter(mean_pct.1 > 0.5) %>%
                top_n(n = 30, wt= mean_log2FC) %>% select(cluster , gene  ,description   ,mean_log2FC,  mean_pct.1,max_pval )


#Markers of cluster 9
markers_9 <-  Markers_all_clusters %>% filter(cluster=="9") %>% pull(gene) %>% head(20)

cytokines = c("CXCL1", "CXCL2", "CXCL3", "CXCL8")

markers_9_select = rev(c("CXCL1", "CXCL2", "CXCL3", "CXCL8", "TNFAIP3","TNFAIP2", "NFKBIA" ,  "NFKBIZ" ,   "NFKB1", "RELB" ,  "IER3" ))
plot <- DotPlot(seurat_integrated,features= markers_9_select) + coord_flip()
plot
ggsave(filename="Dotplot_markers_9.pdf",path=plot_result_path,width= 5.7,height= 3.6)

data <- plot$data
head(data)
write.table(data, file.path(Figures_data,"Figure_6B_Dotplot_cluster_9_markers.txt" ), sep = "\t", quote = F)



#FeaturePlot(seurat_integrated,features= cytokines ,label=F, ncol = 2, raster = T, pt.size = 3)
#ggsave(filename="FeaturePlot_markers_9.pdf",path=plot_result_path,width= 13000,height= 9000,dpi=1500,units="px",device="pdf")





#GO ANALYSIS
#Prepare for GO enrichment analysis
Markers_all_clusters_select <- Markers_all_clusters %>% 
                              filter(mean_log2FC > 0.5 ,  mean_pct.1 > 0.5, max_pval < 0.05 ) %>% select(cluster , gene  ,description   ,mean_log2FC,max_pval )


Markers_genes_list <- split(Markers_all_clusters_select$gene,Markers_all_clusters_select$cluster)
sapply( Markers_genes_list, length)


columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
All_genes<- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(seurat_integrated),keytype='SYMBOL',column=c('ENSEMBL'),muliVals=first)
All_genes<-All_genes[!(is.na(All_genes))]


H<-msigdbr(species = 'Homo sapiens',category='H')
H.ensemble<-H%>%dplyr::select(gs_name,ensembl_gene)

#All cluster
Find_GO<-function(cluster){
  markers_genes <-AnnotationDbi::mapIds(org.Hs.eg.db,keys=Markers_genes_list[[cluster]],keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
  markers_genes <- markers_genes[!(is.na(markers_genes))]
  enrich.cluster <- enricher(markers_genes,universe = All_genes,TERM2GENE = H.ensemble)
  enrich.cluster_df = data.frame()
  if (!(is.null(enrich.cluster ))){
  enrich.cluster_df <-enrich.cluster@result %>% cbind(cluster_id=cluster,.)}
  enrich.cluster_df
}

Markers_genes_list <- Markers_genes_list[c("1","2","3","4","5","6","7","8", "9")]
enrich_GO_all_clusters<- purrr::map_dfr(names(Markers_genes_list),Find_GO)


range(enrich_GO_all_clusters$p.adjust)


enrich_GO_all_clusters_wide<-enrich_GO_all_clusters%>%pivot_wider(id_cols='cluster_id',
                                                                  names_from = 'Description',
                                                                  values_from = 'p.adjust',values_fill=1)
nb_col = ncol(enrich_GO_all_clusters_wide)
colnames(enrich_GO_all_clusters_wide)[2:nb_col]<-gsub("HALLMARK_","", colnames(enrich_GO_all_clusters_wide)[2:nb_col])
colnames(enrich_GO_all_clusters_wide)[2:nb_col]<-gsub("_"," ",colnames(enrich_GO_all_clusters_wide)[2:nb_col])

pathways_keep = apply(enrich_GO_all_clusters_wide[,2:nb_col],2,min) <= 0.05
pathways_keep = c("cluster_id",names(pathways_keep[pathways_keep==TRUE]))
enrich_GO_all_clusters_wide  = enrich_GO_all_clusters_wide[,pathways_keep] 



annotation_clusters = as.data.frame(enrich_GO_all_clusters_wide[1])
enrich_GO_all_clusters_wide = enrich_GO_all_clusters_wide[,c(2:ncol(enrich_GO_all_clusters_wide))]

library(pheatmap)
library(RColorBrewer)
library(brewer.pal)

matrix<-t(enrich_GO_all_clusters_wide)
matrix <- - log10(matrix)
colnames(matrix)<-annotation_clusters$cluster_id
colnames(matrix)
rownames(matrix)

cluster_order<-c('1','2','3','4','5','6','7', '8' , '9')

pathways = c(  "E2F TARGETS"  , "G2M CHECKPOINT" , "DNA REPAIR",   
              "MYC TARGETS V1","MYC TARGETS V2"  , "MITOTIC SPINDLE"   ,"SPERMATOGENESIS",  
               "HYPOXIA"  , "INTERFERON GAMMA RESPONSE"  ,   "UV RESPONSE UP" , 
              "IL6 JAK STAT3 SIGNALING", "INFLAMMATORY RESPONSE"  , "TNFA SIGNALING VIA NFKB"  )


matrix <- matrix[pathways,cluster_order]
               
pheatmap(matrix,cluster_rows = F,
         color = colorRampPalette(c("white",brewer.pal(n = 9, name ="Reds")[c(1,2,4,7)]))(100),
         breaks = seq(0, 2, 0.02),
         legend_breaks = c(0, 0.5, 1, 1.5, 2 , 2.5 ),
         cluster_cols = F,
         annotation_names_col = T,annotation_legend = F,show_colnames = T,fontsize = 12,
         filename=file.path(plot_result_path,"GO_cluster_heatmap.pdf"),
         width = 5.8 , height = 4)

write.table(matrix, file.path(Figures_data, "Figure_5C_GO_transcriptional_clusters.txt"), sep = "\t", quote = F)



#Plot expression
meta <- seurat_integrated@meta.data
meta$CB_seurat <- rownames(meta)
meta <- meta %>%
  mutate(treatment = case_when(
    grepl("-no", orig.ident) ~ "Untreated",
    grepl("-vitro", orig.ident) ~ "Ex vivo",
    grepl("-vivo", orig.ident) ~ "In vivo",
    TRUE ~ orig.ident)
  )
meta <- meta %>% mutate(patient = sub("-.*","",orig.ident))

seurat_integrated@meta.data <- meta
rownames(seurat_integrated@meta.data) <- rownames(meta)

norm_data = seurat_integrated$SCT$data

#ALDH1 
Genes = c("ALDH1A1" )
norm_data_genes = as.data.frame(t(norm_data[Genes,]))
norm_data_genes$gene = "ALDH1A1" 
norm_data_genes_longer = norm_data_genes %>%  pivot_longer(cols = 1:43366,names_to = "CB_seurat", values_to = "Expression")
norm_data_genes_longer = norm_data_genes_longer %>% left_join(meta[, c("CB_seurat","patient", "treatment")], by = "CB_seurat")
norm_data_genes_longer$treatment = factor(norm_data_genes_longer$treatment, levels = c("Untreated", "Ex vivo", "In vivo"))
norm_data_genes_longer$patient = factor(norm_data_genes_longer$patient, levels = c("OO77", "OO100", "OO99"))

norm_data_genes_longer %>% ggplot(aes(x = treatment,y= Expression, fill= treatment ))+
 geom_boxplot( width = 0.1, outlier.size = 0.5)+
 geom_violin(
    trim = TRUE, 
    width = 0.5, 
    alpha = 0.5, 
    color = "black"
  ) +   
 scale_y_continuous(limits = c(-1,9), breaks = c(0,2,4,6,8))+
ylab("ALDH1A1 gene expression")+ 
 scale_fill_manual(values = c("deepskyblue",  "tomato", "plum"),labels= c("no" = "Untreated","vitro"= "Ex vivo", "vivo"="In vivo"))+
 geom_jitter(size = 0.1, alpha = 0.01, aes(fill = treatment))+
  theme_classic() +                         
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17),
        strip.text = element_text(size = 15),
        panel.spacing = unit(1.5, "lines") ,
         strip.background = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 15))   + facet_wrap(vars(patient), nrow = 3)
ggsave('ALDH1A1_Violin.pdf',path = plot_result_path, width = 4.5, height = 7.5)

#Perform statistic
norm_data_genes_longer %>%
  group_by(patient, gene, treatment) %>%
  summarise(mean = mean(Expression))


results <- norm_data_genes_longer %>%
  group_by(patient, gene) %>%
  summarise(
    ttest_ex_vivo = list(t.test(Expression[treatment == "Untreated"], 
                                Expression[treatment == "Ex vivo"], 
                                alternative = "less",
                                var.equal = TRUE)),
    ttest_in_vivo = list(t.test(Expression[treatment == "Untreated"], 
                                Expression[treatment == "In vivo"], 
                                alternative = "less",
                                var.equal = TRUE))
  ) %>%
  mutate(
    p_value_ex_vivo = sapply(ttest_ex_vivo, function(x) x$p.value),
    p_value_in_vivo = sapply(ttest_in_vivo, function(x) x$p.value)
  ) %>%  mutate(
    sig_level_untreated_vs_exvivo = case_when(
      p_value_ex_vivo < 0.001 ~ "***",
      p_value_ex_vivo < 0.01 ~ "**",
      p_value_ex_vivo < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_level_untreated_vs_invivo = case_when(
      p_value_in_vivo < 0.001 ~ "***",
      p_value_in_vivo< 0.01 ~ "**",
      p_value_in_vivo < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>%
  select(patient, gene, sig_level_untreated_vs_exvivo , sig_level_untreated_vs_invivo)



#Feature plot for ALDH1A1

DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 8 * 1024^3)  
seurat_integrated <- SCTransform(seurat_integrated, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = F)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)




#plot ALDH1A1
FeaturePlot(seurat_integrated,features= "ALDH1A1" ,label=F, raster = TRUE, pt.size = 2)
ggsave(filename="FeaturePlot_ALDH1A1.pdf",path=plot_result_path,width= 5,height= 4)

#Try CD44
VlnPlot(seurat_integrated,features= c("CD44"), group.by= "orig.ident")
ggsave(filename="VlnPlot_CSC.pdf",path=plot_result_path,width= 9000,height= 9000,dpi=1500,units="px",device="pdf")


#Cluster 9
seurat_sub_9 <- subset(seurat_integrated, subset = seurat_clusters == 9)


DefaultAssay(seurat_sub_9)<-'RNA'
seurat_sub_9 <- SCTransform(seurat_sub_9 )

Genes =  c("CXCL1","CXCL2","CXCL3", "CXCL8")
Idents(seurat_sub_9)= "orig.ident"

norm_data = seurat_sub_9$SCT$data
norm_data_genes = as.data.frame(norm_data[Genes,])
norm_data_genes$gene = rownames(norm_data_genes)
norm_data_genes_longer = norm_data_genes %>%  pivot_longer(cols = 1:816,names_to = "CB_seurat", values_to = "Expression")
norm_data_genes_longer = norm_data_genes_longer %>% left_join(meta[, c("CB_seurat","patient", "treatment")], by = "CB_seurat")
norm_data_genes_longer$treatment = factor(norm_data_genes_longer$treatment, levels = c("Untreated", "Ex vivo", "In vivo"))
norm_data_genes_longer$patient = factor(norm_data_genes_longer$patient, levels = c("OO77", "OO100", "OO99"))
norm_data_genes_longer$gene = factor(norm_data_genes_longer$gene, levels = Genes)

plot <- norm_data_genes_longer %>% ggplot(aes(x = treatment,y= Expression,colour = treatment ))+
  geom_violin(
    trim = TRUE, 
    width = 0.7, 
    alpha = 0.5, 
    aes(color = treatment),
  ) +
  geom_boxplot(
    width = 0.15, 
    position = position_dodge(width = 0.75), 
    aes(color = treatment), 
    fill = "white", 
    alpha = 0.7,
    outlier.size = 0.1,
  ) +
 scale_color_manual(values = c("deepskyblue2","tomato","plum"))+
 ylim(-0.5,6)+
 geom_jitter(size = 0.1, alpha = 0.1, aes(fill = treatment))+
 ylab("Normalized expression")+
 theme_classic()+
  theme(legend.position = "bottom", 
       legend.title = element_blank(),
       legend.text = element_text(size = 28, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        strip.text = element_text(size = 28, color = "black"),
        strip.background = element_blank(),
        panel.spacing = unit(1.1, "lines") ) + facet_grid(col = vars(patient), row= vars(gene))

plot 
ggsave('Ligand_9_Violin.pdf',path = plot_result_path, width = 8.5, height = 12)

#Save Figure data
data = plot$data
head(data)

write.table(data, file.path(Figures_data, "Figure_6D_cytokines_expression_across_samples.txt"), sep = "\t", quote = F )


# Perform one-tailed t-tests for each combination of treatment and gene, for each patient
p_values <- norm_data_genes_longer %>% 
  group_by(patient, gene) %>%
  summarise(
    pvalue_untreated_vs_exvivo = t.test(Expression[treatment == "Ex vivo"],Expression[treatment == "Untreated"], alternative = "less")$p.value,
    pvalue_untreated_vs_invivo = t.test(Expression[treatment == "In vivo"],Expression[treatment == "Untreated"], alternative = "less")$p.value
  ) %>%  mutate(
    sig_level_untreated_vs_exvivo = case_when(
      pvalue_untreated_vs_exvivo < 0.001 ~ "***",
      pvalue_untreated_vs_exvivo < 0.01 ~ "**",
      pvalue_untreated_vs_exvivo < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_level_untreated_vs_invivo = case_when(
      pvalue_untreated_vs_invivo < 0.001 ~ "***",
      pvalue_untreated_vs_invivo < 0.01 ~ "**",
      pvalue_untreated_vs_invivo < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>% select(patient, gene, sig_level_untreated_vs_exvivo ,sig_level_untreated_vs_invivo)












#other genes
Genes_bis =  c("NFKBIA","RELB","IER3" )

VlnPlot(object = seurat_sub_9, features = Genes_bis , cols = rep(c("deepskyblue",  "tomato", "plum"),3), ncol = 1) + 
theme(legend.position = "bottom")
ggsave('Markers_9_Violin.pdf',path = plot_result_path, width = 5, height = 8)




#Receptor of ligand cluster 9

Genes =  c("SDC1","SDC3","ADRA2A","DPP4")


norm_data = seurat_integrated$SCT$data
norm_data_genes = as.data.frame(norm_data[Genes,])
norm_data_genes$gene = rownames(norm_data_genes)
norm_data_genes_longer = norm_data_genes %>%  pivot_longer(cols = 1:(ncol(norm_data_genes)-1),names_to = "CB_seurat", values_to = "Expression")
norm_data_genes_longer = norm_data_genes_longer %>% left_join(meta[, c("CB_seurat","patient", "treatment", "seurat_clusters_new")], by = "CB_seurat")
norm_data_genes_longer$treatment = factor(norm_data_genes_longer$treatment, levels = c("Untreated", "Ex vivo", "In vivo"))
norm_data_genes_longer$patient = factor(norm_data_genes_longer$patient, levels = c("OO77", "OO100", "OO99"))
norm_data_genes_longer$gene = factor(norm_data_genes_longer$gene, levels = Genes)

summary_by_cluster = norm_data_genes_longer  %>% group_by(gene,patient,treatment,seurat_clusters_new) %>% 
                                                  summarize(mean = mean(Expression)) %>% 
                                                  rename(seurat_clusters = seurat_clusters_new)

custom_palette <- c(  "#80B1D3", "limegreen", "#FCCDE5" ,
                  "#FDB462","#D9D9D9", "#CCEBC5" ,
                 "#BC80BD" , "khaki1",  "brown1")   




ggplot(summary_by_cluster %>% filter(treatment != "Ex vivo"), aes(x = treatment, y = mean, group = seurat_clusters, color = seurat_clusters)) +
  geom_point( size = 2, position = position_dodge(width = 0.1)) +
  geom_line(aes(group = seurat_clusters, color = seurat_clusters), linewidth = 1.5, alpha = 1, position = position_dodge(width = 0)) +
  ylab("Mean Expression") +
  scale_color_manual(values = custom_palette  ) +
  #geom_text(aes(label = seurat_clusters), vjust = -0.5, position = position_dodge(width = 0.5), size = 3) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text =  element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 22 , color = "black"),
        axis.title.y = element_text(size = 22, color = "black" ),
        axis.text.y = element_text(size = 22, color = "black"),
        strip.background = element_blank(),                 # Remove box
        strip.text = element_text(size = 22, color = "black"),  #
        panel.spacing = unit(1.5 , "lines"))  +
  facet_grid(rows = vars(gene), cols =vars(patient)
             , scales = "free"
             )

ggsave("Receptors_Untreated_vivo.pdf",path=plot_result_path,width=13,height= 10)


write.table(summary_by_cluster , file.path(Figures_data,"Figure_Sup_9B_LinePlot_Receptor_expression.txt" ), sep = "\t", quote = F)

