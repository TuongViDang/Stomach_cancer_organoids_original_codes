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
library(GSEABase)
library(fgsea)
library(clusterProfiler)
library(SingleCellExperiment)
library(scater)
library(ggpubr)
library(ggtern)
library(variancePartition)
library(edgeR)

#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_path_3 <-file.path(main_path,'rds',"three_patients")
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <- file.path(main_path,'txt_files', "three_patients")
dir.create(txt_files_3)
plot_3_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_3_path,"by_patient_separate_in_ex")
dir.create(plot_result_path)
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)

#Load data
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))
seurat_integrated


#Before integration
DefaultAssay(seurat_integrated)<-'RNA'
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
seurat_integrated <- SCTransform(seurat_integrated, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)
seurat_integrated  <- RunPCA(seurat_integrated, npcs = 30, )
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

#Edit the meta data
seurat_integrated$patient<- str_extract(seurat_integrated$orig.ident,"^[^[_-]]*")
seurat_integrated$treatment<-str_extract(seurat_integrated$orig.ident,"[^-]*$")
meta <- seurat_integrated@meta.data
meta <- meta %>%
  mutate(Sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ), 
   Treatment = case_when(
    treatment == "no" ~ "Untreated",
    treatment == "vitro" ~ "Ex vivo",
    treatment == "vivo" ~ "In vivo"
   ))  %>% rename(Patient = patient)

seurat_integrated@meta.data <- meta



seurat_integrated$Sample = factor(seurat_integrated$Sample, levels = c("OO77 Untreated",  "OO77 Ex vivo", "OO77 In vivo",
                                  "OO99 Untreated", "OO99 Ex vivo", "OO99 In vivo",
                                  "OO100 Untreated", "OO100 Ex vivo", "OO100 In vivo" ))
Idents(seurat_integrated) <- "Sample"


color_palette <- c(brewer.pal(12, "Set3")[1:8],  brewer.pal(8, "Set2")[4])

plot <- DimPlot(seurat_integrated, reduction = "umap",label=TRUE,label.size = 8, cols= color_palette, repel=T, raster=T, pt.size = 1.5)  + theme_void() +
        theme(legend.text = element_text( size = 15))
plot
ggsave(file.path(plot_result_path, 'bf_integration_by_patients.pdf'),width = 10,height = 8)


#save data for this Figure
data = plot$data
colnames(data)[3] = "Sample"

head(data)
write.table(data , file.path(Figures_data, "Figure_1B_UMAP_before_integration_all_samples.txt"), sep = "\t", quote = F)

#Plot heatmap for top 200 most variable genes

expr_mat <- seurat_integrated$RNA$counts
metadata <- seurat_integrated@meta.data[,c("Patient","Treatment")]
# Normalize counts using voom
dge <- DGEList(counts = expr_mat)
dge <- calcNormFactors(dge)
expr_voomed <- voom(dge)$E  # Log2-normalized expression

# Assuming expr_voomed is your log2-transformed expression data
gene_variances <- apply(expr_voomed, 1, var)  # Variance for each gene (rows = genes, columns = samples)

# Order the genes by variance (highest variance first)
ordered_genes_var <- sort(gene_variances, decreasing = TRUE)

#Heatmap of top 200 variable genes
top_genes <- names(ordered_genes_var)[1:200]
expression_matrix <- seurat_integrated@assays$SCT$data
expression_matrix <- expression_matrix[rownames(expression_matrix)  %in% top_genes,]
expression_matrix <- expression_matrix[match(top_genes ,rownames(expression_matrix)),]


annotation <- data.frame( 
                    "Treatment"  = seurat_integrated@meta.data['Treatment'],
                    "Patient" = seurat_integrated@meta.data['Patient'])


# Define colors for annotations
annotation_colors <- list(
  Patient = c("OO77" = "goldenrod3", "OO99" =  "darkcyan", "OO100" = "indianred1"),
  Treatment = c("Untreated" = "deepskyblue", "Ex vivo" = "tomato", "In vivo" = "plum"))

# Generate heatmap
#Use viridis color 

plot <- pheatmap(expression_matrix,
         color = viridis::plasma(100),  
         breaks = seq(-2, 2, length.out = 101),
         annotation_col = annotation,  
         annotation_colors = annotation_colors, 
         scale = "row",
         show_colnames = FALSE, show_rownames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         filename = file.path(plot_result_path, "top_200_variable_genes_patient_all_sample_heatmap_viridis_plasma.png"),
         width = 10, height = 6.5)

#save data for this Figure
write.table(expression_matrix , file.path(Figures_data, "Figure_1C_Gene_expression_top_200_genes.txt"), sep = "\t", quote = F)
write.table(annotation , file.path(Figures_data, "Figure_1C_Cell_annotation.txt"), sep = "\t", quote = F)

### Fit variance partition separately for in vivo and ex vivo 

#Ex vivo
seurat_ex <- subset(seurat_integrated, subset = Treatment %in% c("Untreated", "Ex vivo") )
expr_mat_ex <- seurat_ex$RNA$counts
metadata_ex <- seurat_ex@meta.data[,c("Patient","Treatment")]
# Normalize counts using voom
dge_ex <- DGEList(counts = expr_mat_ex)
dge_ex <- calcNormFactors(dge_ex)
expr_voomed_ex <- voom(dge_ex)$E  # Log2-normalized expression

# Assuming expr_voomed is your log2-transformed expression data
gene_variances_ex <- apply(expr_voomed_ex, 1, var)  # Variance for each gene (rows = genes, columns = samples)

# Order the genes by variance (highest variance first)
ordered_genes_var_ex <- sort(gene_variances_ex, decreasing = TRUE)

# Select the top N most variated genes (e.g., top 2000 genes)
top_genes_ex <- names(ordered_genes_var_ex)[1:2000]
# Subset the expression data to include only the top variated genes
expr_voomed_top_ex <- expr_voomed_ex[top_genes_ex, ]


#Formular: both patient and treatment as fixed effect
formula <- ~ as.factor(Treatment) + as.factor( Patient)  +  as.factor (Patient): as.factor(Treatment)


varPart_ex <- fitExtractVarPartModel(expr_voomed_top_ex, formula, metadata_ex)
colnames(varPart_ex) = c( "Treatment","Patient", "Interaction", "Residuals")
varPart_ex$gene = rownames(varPart_ex)
varPart_ex_long = varPart_ex %>% pivot_longer( cols = 1:4,names_to= "Factor", values_to = "Variance" ) %>% select(gene, Factor, Variance)
varPart_ex_long$Variance = varPart_ex_long$Variance * 100
varPart_ex_long$Factor = factor(varPart_ex_long$Factor , levels = c("Patient", "Treatment", "Interaction", "Residuals"))

cols = c("mediumturquoise", "hotpink", "gold", "gray")
plot <- varPart_ex_long %>% ggplot(aes( x= Factor, y = Variance, group = Factor, fill = Factor )) +
            geom_boxplot( width = 0.3, outlier.size = 0.5)+
            geom_violin( width = 1.1, 
                        alpha = 0.5) +    
            geom_jitter(size = 0.1, alpha = 0.05, aes(fill = Factor))+
            scale_fill_manual(values = cols)+
            ylab("% Variance explained")+
            theme_classic() +
            theme(legend.position = "None",
                  axis.text.x.bottom = element_text(size = 20, color = "black"),
                  axis.title.x.bottom = element_text(size = 20, color = "black"),
                  axis.text.y.left = element_text(size = 20, color = "black"),
                  axis.title.y.left = element_text(size = 20, color = "black"))
plot
ggsave(file.path(plot_result_path , "variancePartition_ex.pdf"), width = 7, height = 3.7)

#Save data for this Figure
data = plot$data
head(data)
write.table(data , file.path(Figures_data, "Figure_1D_VariancePartition_Ex_vivo.txt"), sep = "\t", quote = F)


# Ternary plot
# Modify so that treatment + patient + interaction = 1
varPart_ex <- varPart_ex %>% as.data.frame %>%
  mutate(Total = Treatment + Patient + Interaction) %>%
  mutate(Treatment = Treatment / Total,
         Patient = Patient / Total,
         Interaction = Interaction / Total) %>%
  select(-Total,-Residuals)  # Drop the total and residuals column

#save this file
write.table(varPart_ex, file.path(txt_files_3, "varPart_ex_vivo.txt"), sep = "\t", quote = F)

#reload
varPart_ex = read.delim(file.path(txt_files_3, "varPart_ex_vivo.txt"))

#color dot
patients_ex = varPart_ex %>% filter(Patient > 0.5, Treatment < 0.5, Interaction < 0.5) 
treatments_ex = varPart_ex %>% filter(Patient < 0.5, Treatment > 0.5, Interaction < 0.5) 
interactions_ex = varPart_ex %>% filter(Patient < 0.5, Treatment < 0.5, Interaction > 0.5) 

#with all genes are labelled 
ggtern(data = varPart_ex, aes(x = Treatment, y = Patient, z = Interaction 
  )) +
  geom_point(size = 0.5, alpha = 1) +  # Plot points
  geom_point( data = patients_ex , size = 1, alpha = 1, color = "mediumturquoise",)+
  geom_point( data = treatments_ex, size = 1, alpha = 1, color = "hotpink")+
  geom_point( data = interactions_ex, size = 1, alpha = 1, color = "gold")+
  geom_text(
    vjust = -1,hjust = 1, size = 1.5, aes(label = gene), color = "red", check_overlap = TRUE) + 
  labs(T = "Patient Effect",
       L = "Treatment Effect",
       R = "Interaction Effect",
       color = "Patient Effect") +
  theme_bw()+  
  theme(aspect.ratio = 1)  # Makes the triangle taller

ggsave(file.path(plot_result_path , "ternary_plot_ex_allgene.pdf"), width = 15, height = 12)



#text label  
patient_label_ex = c()
treatment_label_ex = c("FABP5",  "ME1", "PPP1R14B")
interaction_label_ex = c( "APBB2", "ID1" ,"HNRNPAB" )


#plot
plot <- ggtern(data = varPart_ex, aes(x = Treatment, y = Patient, z = Interaction  )) +
  geom_point(size = 1, alpha = 1, color = "gray") +  # Plot points
  geom_point( data = patients_ex , size = 1, alpha = 1, color = "mediumturquoise",)+
  geom_point( data = treatments_ex, size = 1, alpha = 1, color = "hotpink")+
  geom_point( data = interactions_ex, size = 1, alpha = 1, color = "gold")+
  geom_text(
    data = varPart_ex %>% filter(gene %in% patient_label_ex),
    vjust = 1.1, hjust = 0.57, size = 4, aes(label = gene), color = "red", check_overlap = TRUE) + 
  geom_text(
    data = varPart_ex %>% filter(gene %in% treatment_label_ex),
    vjust = 1.5, hjust = 0.5, size = 4, aes(label = gene), color = "hotpink", check_overlap = TRUE) + 
  geom_text(
    data = varPart_ex %>% filter(gene %in% interaction_label_ex),
    vjust = -0.6, hjust = 1, size = 4, aes(label = gene), color = "gold", check_overlap = TRUE) + 
  # Circle intersect genes with in vivo (below)
  geom_text(
    data = varPart_ex %>% filter(gene %in% setdiff(intersect_genes,c("EPB41L2", "POF1B","ICA1","LGR4"))),
    vjust = -1, hjust = 1.1, size = 4, aes(label = gene), color = "red", check_overlap = TRUE) + 
  geom_text(
    data = varPart_ex %>% filter(gene %in% c("EPB41L2", "POF1B" )),
    vjust = 1.9, hjust = 0.5, size = 4, aes(label = gene), color = "red", check_overlap = TRUE)+
  geom_text(
    data = varPart_ex %>% filter(gene %in% c("ICA1","LGR4")),
    vjust = 1.8, hjust = 0.5, size = 4, aes(label = gene), color = "red", check_overlap = TRUE)+
  geom_point(data = varPart_ex %>% filter(gene %in% intersect_genes ),
             size = 3, shape = 21, stroke = 0.5, fill = NA, color = "red") +
  labs(T = "Patient Effect",
       L = "Treatment Effect",
       R = "Interaction Effect",
       color = "Patient Effect") +
  theme_bw()+  
    theme(
  aspect.ratio = 0.8,
  axis.text = element_text(color = "black"),  # general axis text
  tern.axis.text.T = element_text(size = 16, color = "black"),  # Top axis
  tern.axis.text.L = element_text(size = 16, color = "black"),  # Left axis
  tern.axis.text.R = element_text(size = 16, color = "black")   # Right axis
)

  theme(aspect.ratio = 0.8)  # Makes the triangle taller


plot
ggsave(file.path(plot_result_path , "ternary_plot_ex.pdf"), width = 9, height = 8.2)


data <- plot$data
write.table(data , file.path(Figures_data, "Figure_1D_TernaryPlot_Ex_vivo.txt"), sep = "\t", quote = F)





#gmt_file <- "/group/poetsch_projects/poetsch_sc/Driver_predict/cancer_gene_lists/h.all.v2023.2.Hs.symbols.gmt" 
#gene_sets <- lapply(hallmark_gs, geneIds)
#names(gene_sets) <- names(hallmark_gs)
#gene_sets <- gene_sets[sapply(gene_sets, function(x) is.character(x) && length(x) > 0)]
#term2gene <- do.call(rbind, lapply(names(gene_sets), function(term) {
#  data.frame(term = term, gene = gene_sets[[term]], stringsAsFactors = FALSE)
#}))

perform_HALLMARK <- function(gene_list){
  gene_list_EN <- AnnotationDbi::mapIds(org.Hs.eg.db,keys= gene_list,keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
  gene_list_EN  <- gene_list_EN[!is.na(gene_list_EN )]
  All_genes <- AnnotationDbi::mapIds(org.Hs.eg.db,keys= rownames(seurat_integrated),keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
  All_genes <- All_genes[!is.na(All_genes)]
  H <- msigdbr(species = 'Homo sapiens',category='H')
  H.ensemble <- H%>%dplyr::select(gs_name,ensembl_gene)
  enrich_HALLMARK <- enricher(gene = gene_list_EN,universe = All_genes,TERM2GENE = H.ensemble)
  enrich_HALLMARK_df <- enrich_HALLMARK@result %>% filter(p.adjust <= 0.05)
  enrich_HALLMARK_df
}

perform_GO <- function(gene_list){
  gene_list_EN <- AnnotationDbi::mapIds(org.Hs.eg.db,keys= gene_list,keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
  gene_list_EN  <- gene_list_EN[!is.na(gene_list_EN )]
  All_genes <- AnnotationDbi::mapIds(org.Hs.eg.db,keys= rownames(seurat_integrated),keytype='SYMBOL',column=c('ENSEMBL'),multiVals='first')
  All_genes <- All_genes[!is.na(All_genes)]
  enrich_BP <- enrichGO(gene_list_EN ,universe = All_genes,keyType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db',ont = 'BP')
  enrich_BP_df <- enrich_BP@result %>% filter(p.adjust <= 0.05, Count >= 5)
  enrich_BP_df
}

#patient genes
patient_ex_genes = patients_ex %>% pull(gene)
perform_HALLMARK(patient_ex_genes )
perform_GO(patient_ex_genes)
#treatment genes
treatment_ex_genes = treatments_ex %>% pull(gene)
perform_HALLMARK(treatment_ex_genes )
perform_GO(treatment_ex_genes )
#interaction genes
interaction_ex_genes = interactions_ex %>% pull(gene)
perform_HALLMARK(interaction_ex_genes )
perform_GO(interaction_ex_genes )

####   In vivo  #####
seurat_in <- subset(seurat_integrated, subset = Treatment %in% c("Untreated", "In vivo") )
expr_mat_in <- seurat_in$RNA$counts
metadata_in <- seurat_in@meta.data[,c("Patient","Treatment")]
# Normalize counts using voom
dge_in <- DGEList(counts = expr_mat_in)
dge_in <- calcNormFactors(dge_in)
expr_voomed_in <- voom(dge_in)$E  # Log2-normalized expression

# Assuming expr_voomed is your log2-transformed expression data
gene_variances_in <- apply(expr_voomed_in, 1, var)  # Variance for each gene (rows = genes, columns = samples)

# Order the genes by variance (highest variance first)
ordered_genes_var_in <- sort(gene_variances_in, decreasing = TRUE)

# Select the top N most variated genes (e.g., top 2000 genes)
top_genes_in <- names(ordered_genes_var_in)[1:2000]
# Subset the expression data to include only the top variated genes
expr_voomed_top_in <- expr_voomed_in[top_genes_in, ]


#Formular: both patient and treatment as fixed effect
formula <- ~ as.factor(Treatment) + as.factor( Patient)  +  as.factor (Patient): as.factor(Treatment)


varPart_in <- fitExtractVarPartModel(expr_voomed_top_in, formula, metadata_in)
colnames(varPart_in) = c( "Treatment","Patient", "Interaction", "Residuals")
varPart_in$gene = rownames(varPart_in)
varPart_in_long = varPart_in %>% pivot_longer( cols = 1:4,names_to= "Factor", values_to = "Variance" ) %>% select(gene, Factor, Variance)
varPart_in_long$Variance = varPart_in_long$Variance * 100
varPart_in_long$Factor = factor(varPart_in_long$Factor , levels = c("Patient", "Treatment", "Interaction", "Residuals"))

cols = c("mediumturquoise", "hotpink", "gold", "gray")
plot <- varPart_in_long %>% ggplot(aes( x= Factor, y = Variance, group = Factor, fill = Factor )) +
            geom_boxplot( width = 0.3, outlier.size = 0.5)+
            geom_violin( width = 1.1, 
                        alpha = 0.5) +    
            geom_jitter(size = 0.1, alpha = 0.05, aes(fill = Factor))+
            scale_fill_manual(values = cols)+
            ylab("% Variance explained")+
            theme_classic() +
            theme(legend.position = "None",
                  axis.text.x.bottom = element_text(size = 20, color = "black"),
                  axis.title.x.bottom = element_text(size = 20, color = "black"),
                  axis.text.y.left = element_text(size = 20, color = "black"),
                  axis.title.y.left = element_text(size = 20, color = "black"))

plot
ggsave(file.path(plot_result_path , "variancePartition_in.pdf"), width = 7, height = 3.7)

#Save data for this Figure
data = plot$data
head(data)
write.table(data , file.path(Figures_data, "Figure_1D_VariancePartition_In_vivo.txt"), sep = "\t", quote = F)



# Ternary plot
# Modify so that treatment + patient + interaction = 1
varPart_in <- varPart_in %>% as.data.frame %>%
  mutate(Total = Treatment + Patient + Interaction) %>%
  mutate(Treatment = Treatment / Total,
         Patient = Patient / Total,
         Interaction = Interaction / Total) %>%
  select(-Total, -Residuals)  # Drop the total and residuals column


#save this file
write.table(varPart_in, file.path(txt_files_3, "varPart_in_vivo.txt"), sep = "\t", quote = F)

#reload
varPart_in = read.delim(file.path(txt_files_3, "varPart_in_vivo.txt"))

#color dot
patients_in = varPart_in %>% filter(Patient > 0.5, Treatment < 0.5, Interaction < 0.5) 
treatments_in = varPart_in %>% filter(Patient < 0.5, Treatment > 0.5, Interaction < 0.5) 
interactions_in = varPart_in %>% filter(Patient < 0.5, Treatment < 0.5, Interaction > 0.5) 


#patient genes
patient_in_genes = patients_in %>% pull(gene)
perform_HALLMARK(patient_in_genes )
perform_GO(patient_in_genes)
#treatment genes
treatment_in_genes = treatments_in %>% pull(gene)
perform_HALLMARK(treatment_in_genes )
perform_GO(treatment_in_genes )
#interaction genes
interaction_in_genes = interactions_in %>% pull(gene)
perform_HALLMARK(interaction_in_genes )
perform_GO(interaction_in_genes )


#Intersect genes driven by treatment/interaction in ex vivo and in vivo
intersect_genes = unique(intersect(c(interaction_in_genes, treatment_in_genes), c(interaction_ex_genes, treatment_ex_genes)))

perform_HALLMARK(intersect_genes)
perform_GO(intersect_genes)
#with all genes are labelled 
ggtern(data = varPart_in, aes(x = Treatment, y = Patient, z = Interaction 
  )) +
  geom_point(size = 0.5, alpha = 1) +  # Plot points
  geom_point( data = patients_in , size = 1, alpha = 1, color = "mediumturquoise",)+
  geom_point( data = treatments_in, size = 1, alpha = 1, color = "hotpink")+
  geom_point( data = interactions_in, size = 1, alpha = 1, color = "gold")+
  geom_text(
    vjust = -1,hjust = 1, size = 1.5, aes(label = gene), color = "red", check_overlap = TRUE) + 
  labs(T = "Patient Effect",
       L = "Treatment Effect",
       R = "Interaction Effect",
       color = "Patient Effect") +
  theme_bw()+  
  theme(aspect.ratio = 1)  # Makes the triangle taller

ggsave(file.path(plot_result_path , "ternary_plot_in_allgene.pdf"), width = 15, height = 12)



#text label  
patient_label_in = c()
treatment_label_in = c("GDA", "MBNL2",  "AC012501.2")
interaction_label_in = c()


#plot
plot <- ggtern(data = varPart_in, aes(x = Treatment, y = Patient, z = Interaction  )) +
  geom_point(size = 1, alpha = 1, color = "gray") +  # Plot points
  geom_point( data = patients_in , size = 1, alpha = 1, color = "mediumturquoise",)+
  geom_point( data = treatments_in, size = 1, alpha = 1, color = "hotpink")+
  geom_point( data = interactions_in, size = 1, alpha = 1, color = "gold")+
  geom_text(
    data = varPart_in %>% filter(gene %in% treatment_label_in),
    vjust = -1,hjust = 0.57, size = 4, aes(label = gene), color = "hotpink", check_overlap = TRUE) + 
  geom_text(
    data = varPart_in %>% filter(gene %in% interaction_label_in),
    vjust = -1, hjust = 0.6, size = 4, aes(label = gene), color = "gold", check_overlap = TRUE) + 
  # Label and circle genes intersect with ex vivo
   geom_text(
    data = varPart_in %>% filter(gene %in% setdiff(intersect_genes,c("ICA1","TNFRSF12A", "EPB41L2" ))),
    vjust = -0.7, hjust = 1, size = 4, aes(label = gene), color = "red", check_overlap = TRUE) + 
  geom_text(
    data = varPart_in %>% filter(gene == "TNFRSF12A"),
    vjust = 0, hjust = 1.2, size = 4, aes(label = gene), color = "red", check_overlap = FALSE) + 
  geom_text(
    data = varPart_in %>% filter(gene == "EPB41L2"),
    vjust = 1.7, hjust = 0.5, size = 4, aes(label = gene), color = "red", check_overlap = FALSE) + 
  geom_text(
    data = varPart_in %>% filter(gene == "ICA1"),
    vjust = 0.7, hjust = 1.2, size = 4, aes(label = gene), color = "red", check_overlap = TRUE) + 
  geom_point(data = varPart_in %>% filter(gene %in% intersect_genes),
             size = 3, shape = 21, stroke = 0.5, fill = NA, color = "red") +
  labs(T = "Patient Effect",
       L = "Treatment Effect",
       R = "Interaction Effect",
       color = "Patient Effect") +
  theme_bw()+  
  theme(
  aspect.ratio = 0.8,
  axis.text = element_text(color = "black"),  # general axis text
  tern.axis.text.T = element_text(size = 16, color = "black"),  # Top axis
  tern.axis.text.L = element_text(size = 16, color = "black"),  # Left axis
  tern.axis.text.R = element_text(size = 16, color = "black")   # Right axis
)

  theme(aspect.ratio = 0.8)  # Makes the triangle taller


plot
ggsave(file.path(plot_result_path , "ternary_plot_in.pdf"), width = 9, height = 8.2)


data <- plot$data
write.table(data , file.path(Figures_data, "Figure_1D_TernaryPlot_In_vivo.txt"), sep = "\t", quote = F)



