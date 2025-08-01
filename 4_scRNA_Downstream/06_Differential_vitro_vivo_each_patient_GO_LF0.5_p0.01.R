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
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_3_path <-file.path(main_path,'plots_3_patients')
plot_path_results <- file.path(plot_3_path,"DE_each_patient")
dir.create(plot_path_results)
functions_path <- '/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)



source(file.path(functions_path,'functions.R'))


seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))
seurat_integrated

table(seurat_integrated$orig.ident)
table(seurat_integrated$seurat_clusters)
#Modify the meta data
meta <- seurat_integrated@meta.data
meta <- meta[,c("orig.ident","Phase","seurat_clusters")]
meta <- meta %>% rowwise %>% mutate(patient = str_split_1(orig.ident,"-")[[1]],
                        treatment = str_split_1(orig.ident,"-")[[2]]) %>% as.data.frame
 
rownames(meta) <- rownames(seurat_integrated@meta.data)                
seurat_integrated@meta.data <- meta


table(seurat_integrated$patient)
table(seurat_integrated$seurat_clusters)

######### Pseudobulk using SCE and cluster of 3 patients
library(DESeq2)

DE_and_FA <- function(patient){
    patient_name = patient
    seurat = subset(seurat_integrated , subset = patient == patient_name)

    DefaultAssay(seurat) <- 'RNA'
    counts_matrix <- GetAssayData(object = seurat, slot = "counts")
    coldata <- seurat@meta.data
    #umap_embeddings<-Embeddings(object = seurat_integrated, reduction = "umap")

    sce <- SingleCellExperiment(assays = list(counts=counts_matrix), colData=coldata)
    #reduceDim(sce,"umap")<-umap_embeddings


    #Aggregate count 
    summed <- aggregateAcrossCells(sce, id=colData(sce)[,c("seurat_clusters","treatment")])
    meta <- colData(summed)  %>% as.data.frame %>% mutate(seurat_clusters = as.character(seurat_clusters))%>%
    mutate(sample_id = paste(orig.ident, seurat_clusters,sep="_"))

    dds <- DESeqDataSetFromMatrix(countData = counts(summed), colData = meta, design = ~ seurat_clusters + treatment)

    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)

    png(file=file.path(plot_path_results,paste(patient, "PCA_Plot_pseudobulk_Deseq2.png")))
    plotPCA(rld, intgroup= "treatment")
    dev.off()


    png(file=file.path(plot_path_results,paste(patient, "PCA_Plot_pseudobulk_Deseq2_cluster.png")))
    plotPCA(rld, intgroup= "seurat_clusters")
    dev.off()

    #Heatmap of sample correlation
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)   
    rownames(rld_cor) <- meta$sample_id
    colnames(rld_cor) <- meta$sample_id

    #reorder the sample
    meta_reorder <- meta %>% arrange(treatment)
    rld_cor_reorder <- rld_cor[meta_reorder$sample_id, meta_reorder$sample_id]
    anno = as.data.frame(meta_reorder[,"treatment"])
    colnames(anno) = c("treatment")
    rownames(anno) = rownames(rld_cor_reorder)

    pheatmap(rld_cor_reorder,annotation = anno,cluster_row=F,cluster_col=F, filename=file.path(plot_path_results,paste(patient,"heatmap_pseudobulk_Deseq2.png")))


    # Run DESeq2 differential expression analysis
    dds <- DESeq(dds)

    png(file=file.path(plot_path_results,paste(patient, "Dispersion_Plot_pseudobulk_Deseq2.png")))
    plotDispEsts(dds)
    dev.off()


    ## Define contrasts, extract results table, and shrink the log2 fold changes
    ##vitro
    contrast_vitro <- c("treatment", "vitro", "no")

    res_table_vitro_unshrunken <- results(dds, contrast=contrast_vitro, alpha = 0.05)

    res_table_vitro <- lfcShrink(dds, contrast=contrast_vitro, res=res_table_vitro_unshrunken,type="normal")

    # Turn the results object into a data frame
    res_vitro_df <- data.frame(res_table_vitro)
	  #Save this table
write.table(res_vitro_df,file.path(txt_files_3,paste0(patient,"_res_vitro_df.txt" )), sep = "\t", quote = F)

    # Subset the significant results
    padj.cutoff = 0.05
    lfc.cutoff= 1
    sig_vitro_res <- filter(res_vitro_df, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% arrange(desc(log2FoldChange )) %>% tibble::rownames_to_column(var="gene")


    ## Extract normalized counts from dds object
    normalized_counts <- counts(dds, normalized = TRUE)

    ## Extract top 20 DEG from resLFC (make sure to order by padj)
    top10_up_genes <- sig_vitro_res %>% head(10) %>%
                      dplyr::pull(gene) %>%
    head(n = 10)
    top10_down_genes <- sig_vitro_res %>% tail(10) %>%
                      dplyr::pull(gene) %>%
    head(n = 10)

    top20_vitro_genes <- c(top10_up_genes,top10_down_genes)


    ## Extract matching normalized count values from matrix
    top20_vitro_counts <- normalized_counts[rownames(normalized_counts) %in% top20_vitro_genes, ]
    colnames(top20_vitro_counts) <- colData(dds)$sample_id

    ## Convert wide matrix to long data frame for ggplot2
    top20_vitro_df <- data.frame(top20_vitro_counts)
    top20_vitro_df$gene <- factor(rownames(top20_vitro_counts), levels = c(top10_up_genes,top10_down_genes))

    top20_vitro_df <- pivot_longer(top20_vitro_df, cols = 1:(ncol(top20_vitro_df)-1),
                     values_to = "value",
                     names_to = "sample_id") %>% 
    data.frame()

    ## add treatment column
    top20_vitro_df <- top20_vitro_df %>% mutate(treatment = str_extract(sample_id, "(?<=\\.)[^_]+"))

    ## Generate plot
  ggplot(top20_vitro_df, aes(y = value, x = treatment, col = treatment)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("treatment") +
  ggtitle("Top 10 genes up and 10 genes downs DE vitro vs no") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)


ggsave(file.path(plot_path_results,paste(patient,"top20_vitro_genes.png")))



##vivo
contrast_vivo <- c("treatment", "vivo", "no")

res_table_vivo_unshrunken <- results(dds, contrast=contrast_vivo, alpha = 0.05)

res_table_vivo <- lfcShrink(dds, contrast=contrast_vivo, res=res_table_vivo_unshrunken,type="normal")

# Turn the results object into a data frame
res_vivo_df <- data.frame(res_table_vivo)
#Save this table
write.table(res_vivo_df,file.path(txt_files_3,paste0(patient,"_res_vivo_df.txt" )), sep = "\t", quote = F)
# Subset the significant results
padj.cutoff = 0.05
lfc.cutoff= 1
sig_vivo_res <- filter(res_vivo_df, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)  %>% arrange(desc(log2FoldChange ))  %>% tibble::rownames_to_column(var="gene")

length(intersect(sig_vivo_res$gene, sig_vitro_res$gene))

## Extract top 20 DEG from resLFC (make sure to order by padj)
top10_up_genes <- sig_vivo_res %>% head(10) %>%
  dplyr::pull(gene) %>%
  head(n = 10)
top10_down_genes <- sig_vivo_res %>% tail(10) %>%
  dplyr::pull(gene) %>%
  head(n = 10)

top20_vivo_genes <- c(top10_up_genes,top10_down_genes)


## Extract matching normalized count values from matrix
top20_vivo_counts <- normalized_counts[rownames(normalized_counts) %in% top20_vivo_genes, ]
colnames(top20_vivo_counts) <- colData(dds)$sample_id

## Convert wide matrix to long data frame for ggplot2
top20_vivo_df <- data.frame(top20_vivo_counts)
top20_vivo_df$gene <- factor(rownames(top20_vivo_counts), levels = c(top10_up_genes,top10_down_genes))

top20_vivo_df <- pivot_longer(top20_vivo_df, cols = 1:(ncol(top20_vivo_df)-1),
                     values_to = "value",
                     names_to = "sample_id") %>% 
  data.frame()

## add treatment column
top20_vivo_df <- top20_vivo_df %>% mutate(treatment = str_extract(sample_id, "(?<=\\.)[^_]+"))

## Generate plot
ggplot(top20_vivo_df, aes(y = value, x = treatment, col = treatment)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("treatment") +
  ggtitle("Top 10 genes up and 10 genes downs DE vivo vs no") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)


ggsave(file.path(plot_path_results,paste(patient,"top20_vivo_genes.png")))



### Extract normalized counts for significant genes only
only_vivo_genes = setdiff(sig_vivo_res$gene,sig_vitro_res$gene)
all_sig_genes = c(sig_vitro_res$gene,only_vivo_genes )


sig_counts <- normalized_counts[rownames(normalized_counts) %in% all_sig_genes, ]

## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

#set colnames of matrix
colnames(sig_counts) <- meta$sample_id

#reorder the sample
sig_counts_reorder <-sig_counts[all_sig_genes , meta_reorder$sample_id]
#use the same anno as the correlation heatmap

## Run pheatmap using the metadata data frame for the annotation
anno_colors = list(treatment = c("no"="deepskyblue2", 
                                "vitro" = "tomato", 
                                "vivo" = "plum"))

pheatmap(sig_counts_reorder, 
         color = heat_colors,
         cluster_col = F, cluster_rows = F, 
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = anno, 
         annotation_colors = anno_colors ,
         border_color = NA, 
         fontsize = 14, 
         scale = "row", 
         fontsize_row = 14, 
         width = 6, height = 7,
         filename=file.path(plot_path_results,paste(patient,"heatmap_all_sig_genes_Deseq2.png")))

write.table(sig_counts_reorder, file.path(Figures_data, paste0( "Figure_Sup_2B_Heatmap_significant_genes_", patient, ".txt")), sep = "\t", quote = F)



#Scatter plot
res_vitro_df$gene<-rownames(res_vitro_df)
res_vivo_df$gene<-rownames(res_vivo_df)
both<-dplyr::full_join(res_vitro_df,res_vivo_df,by="gene",suffix=c(".vitro",".vivo"))

sum(is.na(both$padj.vitro))
both$padj.vitro[is.na(both$padj.vitro)]=1
sum(is.na(both$padj.vitro))

sum(is.na(both$padj.vivo))
both$padj.vivo[is.na(both$padj.vivo)]=1
sum(is.na(both$padj.vivo))

sum(is.na(both$log2FoldChange.vitro))
both$log2FoldChange.vitro[is.na(both$log2FoldChange.vitro)]=0
sum(is.na(both$log2FoldChange.vitro))

sum(is.na(both$log2FoldChange.vivo))
both$log2FoldChange.vivo[is.na(both$log2FoldChange.vivo)]=0
sum(is.na(both$log2FoldChange.vivo))

both<-both%>%
          mutate(vitro_sig=ifelse(padj.vitro <= 0.05 & abs(log2FoldChange.vitro) >= 1 ,"yes","no"), 
                vivo_sig=ifelse(padj.vivo <= 0.05 & abs(log2FoldChange.vivo) >= 1,"yes","no"))%>%
                mutate(sig=case_when(
                vitro_sig=="yes" & vivo_sig=="yes" ~ 'both',
                vitro_sig=='yes' & vivo_sig=='no' ~ 'vitro',
                vitro_sig=='no' & vivo_sig=='yes' ~ 'vivo',
                TRUE~'none'))
correlation_value <- cor(both$log2FoldChange.vitro, both$log2FoldChange.vivo)

both%>%ggplot(aes(x=log2FoldChange.vitro,,y=log2FoldChange.vivo))+
    geom_point(aes(color=sig),size= 1.5 )+
    scale_colour_manual(labels=c("both"="both","none"="none","vitro"="ex vivo","vivo"="in vivo"),
                        values=c("both"="darkseagreen1","none"="grey","vitro"="coral1","vivo"="plum"),name="DE")+
    scale_x_continuous(limits=c(-5,5))+
    scale_y_continuous(limits=c(-12,10))+
    geom_text(x = 5,y = 9.5,label = paste("r=", round(correlation_value, 4)), hjust = 1,vjust = 0, size= 7.5,  fontface ="plain")+
    theme_classic()+  
    ylab("LogFC In vivo vs Untreated")+
    xlab("LogFC Ex vivo vs Untreated")+
    theme(axis.title.x.bottom=element_text(size=22, color = "black"),
          axis.text.x.bottom=element_text(size=22, color = "black"),
          axis.title.y.left=element_text(size=22, color = "black"),
          axis.text.y.left = element_text(size=22, color = "black"),
          legend.position = "None")
  

ggsave(paste(patient,"_pseudobulk_Deseq2_scatter.pdf"),path=plot_path_results,width=5,height=5)

#save both 
write.table(both, file.path(txt_files_3 , paste0(patient,"_DE_significant_both.txt")), sep = "\t", quote = F)

write.table(both, file.path(Figures_data, paste0( "Figure_Sup_2A_ScatterPlot_Log2FC_", patient, ".txt")), sep = "\t", quote = F)


#######FA analysis###########
#vitro
vitro_genes <- rownames(res_vivo_df)
vitro_enriched_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = vitro_genes, keytype = 'SYMBOL', column = 'ENTREZID', multiVals = 'first')
vitro_enriched_genes <- vitro_enriched_genes[!is.na(vitro_enriched_genes)]
vitro_enriched_genes <- vitro_enriched_genes[which(duplicated(vitro_enriched_genes) == F) ]

res_vitro_df_ENT <-res_vitro_df[names(vitro_enriched_genes),]
res_vitro_df_ENT$entrezID <- vitro_enriched_genes

## Extract the foldchanges
vitro_foldchanges <- res_vitro_df_ENT$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(vitro_foldchanges) <- res_vitro_df_ENT$entrezID

## Sort fold changes in decreasing order
vitro_foldchanges <- sort(vitro_foldchanges, decreasing = TRUE)




#Use HALLMARK
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
msig_GSEA_vitro <- GSEA(vitro_foldchanges, TERM2GENE = m_t2g, verbose = FALSE , pvalueCutoff = 1)
msig_GSEA_vitro_result <- msig_GSEA_vitro@result %>% arrange(desc(NES))


# vitro_top_NES <- msig_GSEA_vitro$Description[order(msig_GSEA_vitro$NES,decreasing=T)][1:2]
# vitro_bottom_NES <- msig_GSEA_vitro$Description[order(msig_GSEA_vitro$NES)][1:2]

# set=vitro_top_NES
# p1=gseaplot(msig_GSEA_vitro, geneSetID = set[1],title = gsub("HALLMARK_","",set[1]))
# p2=gseaplot(msig_GSEA_vitro, geneSetID = set[2],title = gsub("HALLMARK_","",set[2]))
# ggarrange(p1,p2,nrow=1,ncol=2)
# ggsave(filename=file.path(plot_path_results,paste(patient, "GSE_vitro_top_plot.png")),width=11,height=8)

# set=vitro_bottom_NES
# p1=gseaplot(msig_GSEA_vitro, geneSetID = set[1],title = gsub("HALLMARK_","",set[1]))
# p2=gseaplot(msig_GSEA_vitro, geneSetID = set[2],title = gsub("HALLMARK_","",set[2]))
# ggarrange(p1,p2,nrow=1,ncol=2)
# ggsave(filename=file.path(plot_path_results,paste(patient,"GSE_vitro_bottom_plot.png")),width=11,height=8)


# gseaplot(msig_GSEA_vitro, geneSetID = "HALLMARK_E2F_TARGETS",title = gsub("HALLMARK_","","HALLMARK_E2F_TARGETS"))
# ggsave(filename=file.path(plot_path_results,"GSE_E2F_vitrpo_plot.png"),width=5.5,height=8)
# gseaplot(msig_GSEA_vitro, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",title = gsub("HALLMARK_","","HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
# ggsave(filename=file.path(plot_path_results,"GSE_TNFA_vitro_plot.png"),width=5.5,height=8)



msig_GSEA_vitro_result$ID <- gsub("HALLMARK_","",msig_GSEA_vitro_result$ID)
msig_GSEA_vitro_result$ID <- factor(msig_GSEA_vitro_result$ID,levels=msig_GSEA_vitro_result$ID )


msig_GSEA_vitro_result%>%ggplot(aes(x=ID,y=NES,fill=ID))+geom_col(width=0.5)+theme_bw()+
theme(axis.text.y.left=element_text(size=12),axis.text.x.bottom=element_blank(),legend.position="bottom",legend.title = element_blank())

ggsave(filename=file.path(plot_path_results,paste(patient,"_NES_vitro_plot.png")),width=11,height=8)

#add gene count and patient, treatment
msig_GSEA_vitro_result$counts<-sapply(msig_GSEA_vitro_result$core_enrichment,function(x){str_split(x,"/")[[1]]%>%length()})
msig_GSEA_vitro_result <- msig_GSEA_vitro_result%>%select("ID","NES","counts","p.adjust", "core_enrichment")%>%mutate(patient= patient,treatment="vitro")


#vivo
vivo_genes <- rownames(res_vivo_df)
vivo_enriched_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = vivo_genes, keytype = 'SYMBOL', column = 'ENTREZID', multiVals = 'first')
vivo_enriched_genes <-vivo_enriched_genes[!is.na(vivo_enriched_genes)]
vivo_enriched_genes <- vivo_enriched_genes[which(duplicated(vivo_enriched_genes) == F) ]

res_vivo_df_ENT <-res_vivo_df[names(vivo_enriched_genes),]
res_vivo_df_ENT$entrezID <- vivo_enriched_genes

## Extract the foldchanges
vivo_foldchanges <- res_vivo_df_ENT$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(vivo_foldchanges) <- res_vivo_df_ENT$entrezID

## Sort fold changes in decreasing order
vivo_foldchanges <- sort(vivo_foldchanges, decreasing = TRUE)




#Use HALLMARK
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
msig_GSEA_vivo <- GSEA(vivo_foldchanges, TERM2GENE = m_t2g, verbose = FALSE, pvalueCutoff = 1)
msig_GSEA_vivo_result <- msig_GSEA_vivo@result  %>% arrange(desc(NES)) 

msig_GSEA_vivo_result$ID <- gsub("HALLMARK_","",msig_GSEA_vivo_result$ID)
msig_GSEA_vivo_result$ID <- factor(msig_GSEA_vivo_result$ID,levels=msig_GSEA_vivo_result$ID )
msig_GSEA_vivo_result %>%
 ggplot(aes(x=ID,y=NES,fill=ID)) + 
 geom_col(width=0.5)+
 theme_bw()+
 theme(axis.text.y.left=element_text(size=12),
      axis.text.x.bottom=element_blank(),
      legend.position="bottom",
      legend.title = element_blank())

ggsave(filename=file.path(plot_path_results,paste(patient,"_NES_vivo_plot.png")),width=12,height=8)

#add gene count and patient, treatment
msig_GSEA_vivo_result$counts<-sapply(msig_GSEA_vivo_result$core_enrichment,function(x){str_split(x,"/")[[1]]%>%length()})
msig_GSEA_vivo_result <- msig_GSEA_vivo_result%>%select("ID","NES","counts", "p.adjust", "core_enrichment")%>%mutate(patient=patient,treatment="vivo")

#summary for each patient
msig_GSEA_both_pseudobulk <- rbind(msig_GSEA_vitro_result,msig_GSEA_vivo_result)

msig_GSEA_both_pseudobulk

}

patients = c("OO100","OO99","OO77")
all_results =  purrr::map_dfr(patients, DE_and_FA)
write.table(all_results, file.path(txt_files_3,"DE_FA_each_patient.txt"), sep = "\t", quote = F)


#Plot FA result
all_results = read.delim(file.path(txt_files_3,"DE_FA_each_patient.txt"))
all_results <- all_results %>% mutate(comparision_ID  = paste(patient,treatment,sep="_"))

#Plot bubble plot for vivo 
all_vivo_results <- all_results %>% filter(treatment == "vivo")
all_vivo_results$comparision_ID <- factor(all_vivo_results$comparision_ID, levels = c("OO77_vivo","OO99_vivo","OO100_vivo"))


all_vivo_results <- all_vivo_results %>% filter( ID %in% ordered_pathways)
all_vivo_results$ID <- factor(all_vivo_results$ID, levels = rev(ordered_pathways ))
rownames(all_vivo_results) <- NULL

#only significant pathways
all_vivo_results_sig = all_vivo_results %>% filter(p.adjust < 0.05)
ggplot(all_vivo_results_sig, aes(x = comparision_ID, y = ID, color = NES, size = counts)) +
  geom_point(alpha = 0.7) + scale_x_discrete(labels = c( "OO77_vivo" = "OO77 in vivo", "OO100_vivo" = "OO100 in vivo","OO99_vivo" = "OO99 in vivo"))+
  scale_color_gradientn( breaks = c(-2,-1,0,1,2),colours = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(7))+
  theme_classic()+
  theme(axis.title.x.bottom=element_blank(), axis.text.x.bottom=element_text(size = 10),
        axis.title.y.left=element_blank(), legend.position = "top")

ggsave(filename=file.path(plot_path_results,paste("Bubble_vivo_plot_sig_FA.png")),width=8,height=7)

#Bubble plot for all
# Add a new column to indicate significant points
all_results <- all_results %>% 
  mutate(significant = ifelse(p.adjust <= 0.05, "Yes", "No"),
         neg_log_adj_p = -log10(p.adjust))  # Calculate -log10(p.adjust)

all_results$ID <- as.character(all_results$ID)
remaining_pathways = setdiff(unique(all_results$ID) , ordered_pathways )
all_ordered_pathways = c(ordered_pathways,remaining_pathways )

all_results$ID <- factor(all_results$ID, levels = rev(all_ordered_pathways ))
rownames(all_results) <- NULL
all_results$comparision_ID <- factor(all_results$comparision_ID, levels = c("OO100_vivo","OO99_vivo","OO77_vivo","OO100_vitro","OO99_vitro","OO77_vitro"))

#Only  ordered_pathways
#all_results <- all_results %>% filter(ID %in% ordered_pathways)
#Plot the buble plot
ggplot(all_results, aes(x = comparision_ID, y = ID, color = NES, size = neg_log_adj_p, shape = significant)) +
  geom_point(alpha = 1, stroke = 1.3) + scale_x_discrete(labels = c( "OO77_vivo" = "OO77 In","OO77_vitro" = "OO77 Ex",
                                                         "OO99_vivo" = "OO99 In ","OO99_vitro" = "OO99 Ex",
                                                         "OO100_vivo" = "OO100 In","OO100_vitro" = "OO100 Ex"
                                                        ))+
  scale_color_gradientn( breaks = c(-2,-1,0,1,2),colours = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(7))+
    scale_shape_manual(values = c("Yes" = 16, "No" = 1)) +  # Use different shapes for significant points
      labs(size = "-log10(p.adjust)") + 
  theme_classic()+
  theme(axis.title.x.bottom=element_blank(), axis.text.x.bottom=element_text(size = 12),
        axis.title.y.left=element_blank(),axis.text.y.left =element_text(size = 12),
         legend.position = "top")

ggsave(filename=file.path(plot_path_results,paste("Bubble_vitro_vivo_all_pathways_plot_FA.png")),width= 12,height= 14)


########## GO analysis ###########

#GO analysis 

Perform_GO_patient <- function(patient){
 DefaultAssay(seurat_integrated)<-"RNA"
 All_genes<- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(seurat_integrated),keytype='SYMBOL',column=c('ENSEMBL'),muliVals=first)
 All_genes<-All_genes[!(is.na(All_genes))]


 H<-msigdbr(species = 'Homo sapiens',category='H')
 H.ensemble<-H%>%dplyr::select(gs_name,ensembl_gene)

 #Function to perform GO with Hallmark and return name of gene set
 perform_GO <- function(genes, universe_genes) {
  enriched_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes, keytype = 'SYMBOL', column = 'ENSEMBL', multiVals = 'first')
  enriched_genes <- enriched_genes[!is.na(enriched_genes)]
  enrich_H <- enricher(enriched_genes, universe = All_genes,TERM2GENE = H.ensemble )
  enrich_H_df <- enrich_H@result %>%  filter(p.adjust <= 0.01) 
  enrich_H_df<-enrich_H_df %>% mutate(Description = gsub("HALLMARK_", "", Description))
  return(enrich_H_df)
}

padj.cutoff = 0.05
lfc.cutoff= 0.5  #0.75
#vitro
res_vitro_df = read.delim(file.path(txt_files_3,paste0(patient,"_res_vitro_df.txt" )))
# Subset the significant results

sig_vitro_res <- filter(res_vitro_df, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)  %>% arrange(desc(log2FoldChange ))  %>% tibble::rownames_to_column(var="gene")


vitro_up_genes = sig_vitro_res %>% filter(log2FoldChange > 0) %>% pull(gene)
GO_vitro_up = perform_GO(vitro_up_genes , universe_genes = All_genes)
if(nrow(GO_vitro_up) > 0) {
  GO_vitro_up$Direction = "Up"} else {GO_vitro_up$Direction <- character() }

vitro_down_genes = sig_vitro_res %>% filter(log2FoldChange < 0) %>% pull(gene)
GO_vitro_down = perform_GO(vitro_down_genes , universe_genes = All_genes)

if(nrow(GO_vitro_down) > 0) {
  GO_vitro_down$Direction = "Down"} else {GO_vitro_down$Direction <- character() }
GO_vitro = rbind(GO_vitro_up, GO_vitro_down)
GO_vitro$Comparison = "Ex"

#vivo
res_vivo_df = read.delim(file.path(txt_files_3,paste0(patient,"_res_vivo_df.txt" )))
sig_vivo_res <- filter(res_vivo_df, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)  %>% arrange(desc(log2FoldChange ))  %>% tibble::rownames_to_column(var="gene")

vivo_up_genes = sig_vivo_res %>% filter(log2FoldChange > 0) %>% pull(gene)
GO_vivo_up = perform_GO (vivo_up_genes , universe_genes = All_genes)

if(nrow(GO_vivo_up) > 0) {
  GO_vivo_up$Direction = "Up"} else {GO_vivo_up$Direction <- character() }
vivo_down_genes = sig_vivo_res %>% filter(log2FoldChange < 0) %>% pull(gene)

vivo_down_genes = sig_vivo_res %>% filter(log2FoldChange < 0) %>% pull(gene)
GO_vivo_down = perform_GO(vivo_down_genes , universe_genes = All_genes)

if(nrow(GO_vivo_down) > 0) {
  GO_vivo_down$Direction = "Down"} else {GO_vivo_down$Direction <- character() }
GO_vivo = rbind(GO_vivo_up, GO_vivo_down)
GO_vivo$Comparison = "In"

GO_patient = rbind(GO_vitro, GO_vivo)
GO_patient$patient = patient

GO_patient
}


patients = c("OO100","OO99","OO77")
GO_results <- purrr::map_dfr(patients, Perform_GO_patient)
GO_results$ID <- sub("HALLMARK_","", GO_results$ID )
GO_results$neg_log_adj_p <- -log10(GO_results$p.adjust)
GO_results$comparision_ID <- paste(GO_results$patient, GO_results$Comparison)

ordered_com = c(  "OO77 Ex", "OO99 Ex","OO100 Ex",  "OO77 In", "OO99 In","OO100 In"  )
GO_ordered_pathway = c(  "TNFA_SIGNALING_VIA_NFKB",  "HYPOXIA"      ,"INTERFERON_ALPHA_RESPONSE" , 
                          "UV_RESPONSE_UP" , 
                          "APOPTOSIS"        ,  "P53_PATHWAY"  ,   "INTERFERON_GAMMA_RESPONSE"  , "ESTROGEN_RESPONSE_EARLY" , "CHOLESTEROL_HOMEOSTASIS",  "APICAL_JUNCTION"   , 
                          "KRAS_SIGNALING_UP"  , "EPITHELIAL_MESENCHYMAL_TRANSITION" ,             
                          "ESTROGEN_RESPONSE_LATE",  "GLYCOLYSIS","XENOBIOTIC_METABOLISM", 
                     "MYOGENESIS",    "COMPLEMENT"  ,  "MTORC1_SIGNALING"  , "MITOTIC_SPINDLE","E2F_TARGETS"  , "G2M_CHECKPOINT"   )

GO_results$ID <- factor(GO_results$ID , levels = rev(GO_ordered_pathway))

GO_results$comparision_ID <- factor(GO_results$comparision_ID , levels = ordered_com)

plot  <-  ggplot(GO_results, aes(x = comparision_ID, y = ID, color = Direction , size = neg_log_adj_p)) +
  geom_point(alpha = 1, stroke = 1.3) + 
  scale_color_manual(values = c( "#440154FF", "#FDE725FF"))+
      labs(size = "-log10(p.adjust)") + 
  theme_bw()+
  theme(axis.title.x.bottom=element_blank(),
        axis.text.x.bottom=element_text(size = 18, color = "black"),
        axis.title.y.left=element_blank(),
        axis.text.y.left = element_text(size = 18, color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), 
        legend.box = "vertical",            # Stack legends vertically
        legend.box.just = "left" ,
        legend.spacing.y = unit(2, "pt")         
  ) +
  guides(
    color = guide_legend(nrow = 1),
    size = guide_legend(nrow = 1)
  )

plot
ggsave(filename=file.path(plot_path_results,paste("Bubble_vitro_vivo_GO_plot_LF0.5_p0.01.pdf")),width= 11,height= 8)

#Save data
data <- plot$data

write.table(data, file.path(Figures_data, "Figure_2B_GO_DE_bubleplot.txt"), sep = "\t", quote = F)


## Genes 
Hallmark_df <- read.delim("/group/poetsch_projects/poetsch_sc/Driver_predict/cancer_gene_lists/h.all.v2023.2.Hs.symbols.gmt")
Hallmark_df[,1] <- sub("HALLMARK_","",Hallmark_df[,1] )


#####Characterise CANCER pathways  ###
#Combined vivo res table of 3 patient
OO77 = read.delim(file.path(txt_files_3,"OO77_res_vivo_df.txt" )) 
OO77$patient = "OO77"
OO77$gene = rownames(OO77)
OO77_sig = OO77 %>% filter(padj< 0.05, abs(log2FoldChange)> 0.5)
OO77_up = OO77_sig %>% filter(log2FoldChange > 0) %>% rownames()
OO77_down = OO77_sig %>% filter(log2FoldChange < 0) %>% rownames()

OO99 = read.delim(file.path(txt_files_3,"OO99_res_vivo_df.txt" )) 
OO99$patient = "OO99"
OO99$gene = rownames(OO99)
OO99_sig = OO99 %>% filter(padj< 0.05, abs(log2FoldChange)> 0.5)
OO99_up = OO99_sig %>% filter(log2FoldChange > 0) %>% rownames()
OO99_down = OO99_sig %>% filter(log2FoldChange < 0) %>% rownames()

OO100 = read.delim(file.path(txt_files_3,"OO100_res_vivo_df.txt" )) 
OO100$patient = "OO100"
OO100$gene = rownames(OO100)
OO100_sig = OO100 %>% filter(padj< 0.05, abs(log2FoldChange)> 0.5)
OO100_up = OO100_sig %>% filter(log2FoldChange > 0) %>% rownames()
OO100_down = OO100_sig %>% filter(log2FoldChange < 0) %>% rownames()

All = rbind(OO77,OO99,OO100)


Plot_HM_dir = file.path(plot_path_results, "HALLMARK_vivo")
dir.create(Plot_HM_dir)


fig_list = list()
intersect_df_list = list()

for (i in c(1:1)){
pathway = GO_ordered_pathway[i]
genes_pathway = Hallmark_df[Hallmark_df[1]== pathway,3:ncol(Hallmark_df)] %>% as.vector()  %>% unlist() %>% as.character()
OO77_up_pathways = intersect(OO77_up ,genes_pathway )
OO99_up_pathways = intersect(OO99_up , genes_pathway )
OO100_up_pathways = intersect(OO100_up , genes_pathway )

OO77_down_pathways = intersect(OO77_down ,genes_pathway )
OO99_down_pathways = intersect(OO99_down , genes_pathway )
OO100_down_pathways = intersect(OO100_down , genes_pathway )

up_intersect = intersect(OO77_up_pathways ,intersect(OO99_up_pathways , OO100_up_pathways))
down_intersect = intersect(OO77_down_pathways ,intersect(OO99_down_pathways , OO100_down_pathways))

intersect = union(up_intersect , down_intersect ) 
# Create a new column to highlight genes in the union and label genes in the intersect
All$highlight <- ifelse(All$gene %in% genes_pathway, TRUE, FALSE)
All$label <- ifelse(All$gene %in% intersect, All$gene, NA)
All$shape <- ifelse(All$gene %in% intersect, TRUE, FALSE)

# Calculate -log10(p-value)
All$neg_log10_pvalue <- -log10(All$padj)
All$neg_log10_pvalue <- ifelse(All$neg_log10_pvalue > 300, 300, All$neg_log10_pvalue )
All$patient = factor(All$patient, levels = c("OO77","OO99","OO100"))
library(ggrepel)
# Volcano plot with facet_wrap
p =ggplot(All, aes(x=log2FoldChange, y= neg_log10_pvalue, color = highlight, shape = shape)) +
  geom_point(data = subset(All, highlight == FALSE),color = "gray85", size= 0.8) +
  geom_point(data = subset( All, highlight == TRUE), color = "hotpink", size = 1.5) +
  geom_point(data = subset( All, shape == TRUE), shape = 21, size = 2) +
  scale_y_continuous(breaks=seq(0,300,100), limits = c(0,300))+
  scale_x_continuous(breaks=seq(-8,8,2), limits = c(-10,10))+
  theme_classic() +
  labs(title= pathway, x= "Log2(Fold Change)", y= "-Log10(p-adj)") +
  theme( plot.title = element_text(size = 13) ,
       legend.position = "none",
        axis.text.y.left = element_text(size = 14), 
        axis.title.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_text(size = 14),
        #strip.text = element_text(size = 16, face = "bold"),
        strip.text = element_blank(),
        strip.background = element_blank(),  
        panel.spacing = unit(1, "lines") ) +
  geom_text_repel(aes(label = label), vjust = 1.5, size = 3)+
  facet_wrap(vars(patient), nrow = 3) 

p
fig_list[[i]]= p
ggsave(file.path(Plot_HM_dir ,paste0(pathway,"_vivo_vs_untreated.pdf")), width = 5, height = 10)

#Prepare expression df
intersect_up = data.frame(Gene = up_intersect ,Direction = rep("Up", length(up_intersect)), Pathway= rep(pathway, length(up_intersect)))
intersect_down = data.frame(Gene = down_intersect ,Direction = rep("Down", length(down_intersect)), Pathway= rep(pathway, length(down_intersect)))
intersect_df_list[[i]] = rbind(intersect_up,intersect_down )
}

# Combine the plots using patchwork
combined_plot = fig_list[[1]] + fig_list[[2]] + fig_list[[3]] +fig_list[[4]]  +  plot_layout(ncol = 4, widths = c(1, 1, 1, 1, 1))  # Adjust ncol or nrow as needed for layout
ggsave(file.path(Plot_HM_dir ,"Cancer_vivo_vs_untreated.pdf"), width = 18, height = 8)


#Heatmap of single cell of similar trend genes
intersect_df = do.call(rbind, intersect_df_list) %>% distinct(Gene, .keep_all = T) %>% arrange(desc(Direction))

seurat_sub = subset(seurat_integrated, subset = orig.ident %in% c("OO100-no","OO100-vivo","OO99-no","OO99-vivo","OO77-no","OO77-vivo"))
DefaultAssay(seurat_sub) = "RNA"
options(future.globals.maxSize = 6 * 1024^3)
seurat_sub = SCTransform(seurat_sub,vst.flavor='v1',
                          variable.features.n = nrow(seurat_sub))


seurat_sub$patient = factor(seurat_sub$patient, levels = c("OO77", "OO99","OO100"))
seurat_sub@meta.data = seurat_sub@meta.data %>% arrange(patient)

expression_mat = seurat_sub$SCT$scale.data
expression_intersect = expression_mat[intersect_df$Gene,rownames(seurat_sub@meta.data )]

anno_gene = intersect_df %>% tibble::column_to_rownames(var ="Gene")
anno_gene = anno_gene[,c("Pathway","Direction")]

meta = seurat_sub@meta.data
meta = meta %>% mutate(Treatment = case_when(treatment == "no" ~ "Untreated",
                                            treatment == "vitro" ~ "Ex vivo",
                                            treatment == "vivo" ~ 'In vivo'))
meta = meta %>% mutate(sample = paste(patient, treatment))                                           
seurat_sub@meta.data = meta
rownames(seurat_sub@meta.data) = rownames(meta)

anno_cell = seurat_sub@meta.data[,"Treatment",drop= F]

pathways_color = brewer.pal(4, "Set3")
names(pathways_color) = unique(anno_gene$Pathway)
annotation_colors <- list( Treatment = c("Untreated" = "deepskyblue2", "In vivo" = "violet"),
                           Pathway = pathways_color ,
                           Direction =  c("Up" = "cornflowerblue", "Down" = "tomato"))



pheatmap::pheatmap(expression_intersect, 
                   color = viridis::viridis(100),
                   show_rownames = TRUE, 
                   show_colnames = FALSE,   
                   breaks = seq(- 1, 1, length.out = 101),
                   cluster_row = FALSE, 
                   cluster_col = FALSE, 
                   gaps_col = c( 8558 , 17879 ),
                   gaps_row = c(2),
                   annotation_row = anno_gene, 
                   annotation_col = anno_cell, 
                   annotation_colors = annotation_colors,
                   filename = file.path(Plot_HM_dir, "Heatmap_gene_common_pathway_viridis.png"),
                   width = 11, 
                   height = 9)
#Only NFKb

#arrange gene by gene expression gap In vivo vs Untreated
expression_intersect_bis <- as.data.frame(t(expression_intersect) )
expression_intersect_bis$CB <- rownames(expression_intersect_bis)
meta$CB = rownames(meta)
expression_intersect_bis <- expression_intersect_bis %>% left_join(meta[,c("CB","patient","Treatment")], by = "CB")
expression_intersect_bis <- expression_intersect_bis %>% pivot_longer(cols = 1:15, values_to = "expression", names_to = "gene")
expression_intersect_bis <- expression_intersect_bis %>% group_by(gene, Treatment) %>% 
                             summarize(mean =  mean(expression)) %>%
                             pivot_wider(names_from = "Treatment", values_from = "mean") %>%
                             mutate(gap = Untreated - `In vivo`) %>%
                             arrange(desc(gap))

gene_order = expression_intersect_bis %>% pull(gene)
expression_intersect = expression_intersect[gene_order,]


pheatmap::pheatmap(expression_intersect, 
                   color = viridis::viridis(100),
                   show_rownames = TRUE, 
                   show_colnames = FALSE,   
                   breaks = seq(- 1, 1, length.out = 101),
                   cluster_row = FALSE, 
                   cluster_col = FALSE, 
                   gaps_col = c( 8558 , 17879 ),
                   annotation_col = anno_cell, 
                   annotation_colors = annotation_colors,
                   filename = file.path(Plot_HM_dir, "Heatmap_gene_common_pathway_viridis_NFkb.png"),
                   width = 7.5, 
                   height = 3.7)

write.table(expression_intersect, file.path(Figures_data ,"Figure_2C_Common_genes_NFKB.txt"), sep = "\'t", quote = F)



#GO for DE genees both ex vivo and in vivo
patient = "OO77"

both_df <- read.delim(file.path(txt_files_3 , paste0(patient,"_DE_significant_both.txt"))) %>% filter(sig == "both")
both_genes <- both_df$gene
All_genes<- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(seurat_integrated),keytype='SYMBOL',column=c('ENSEMBL'),muliVals=first)
All_genes<-All_genes[!(is.na(All_genes))]


H <- msigdbr(species = 'Homo sapiens',category='H')
H.ensemble <- H%>%dplyr::select(gs_name,ensembl_gene)


 #Function to perform GO with Hallmark and return name of gene set
 perform_HM <- function(genes, universe_genes) {
  enriched_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes, keytype = 'SYMBOL', column = 'ENSEMBL', multiVals = 'first')
  enriched_genes <- enriched_genes[!is.na(enriched_genes)]
  enrich_H <- enricher(enriched_genes, universe = All_genes,TERM2GENE = H.ensemble )
  enrich_H_df <- enrich_H@result %>%  filter(p.adjust <= 0.05) 
  enrich_H_df<-enrich_H_df %>% mutate(Description = gsub("HALLMARK_", "", Description))
  return(enrich_H_df)
}


 # Perform GO enrichment
perform_GO <- function(genes, universe_genes) {
  # Convert gene symbols to Entrez IDs (required for clusterProfiler GO functions)
  entrez_genes <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  entrez_genes <- entrez_genes[!is.na(entrez_genes)]
  # Convert universe gene symbols to Entrez IDs
  universe_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  universe_entrez <- universe_entrez[!is.na(universe_entrez)]
  enrich_GO <- enrichGO(
    gene =   entrez_genes  ,
    universe = universe_entrez ,
    OrgDb = org.Hs.eg.db,
    ont = "All",           # "BP" for Biological Process, "MF", or "CC"
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  enrich_GO_df <- enrich_GO@result %>%  filter(p.adjust <= 0.05) 
  return(enrich_GO_df)
}

perform_GO(both_genes,rownames(seurat_integrated) )