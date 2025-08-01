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
library(diptest)

patient="OO77"


Integrate_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA",patient,"tree30_change_name")
plot_dir = file.path(Integrate_dir,"plots","Map_validate_3_samples")
dir.create(plot_dir, recursive = T)
txt_files = file.path(Integrate_dir,"txt_files")
dir.create(txt_files,recursive=T)
rds_dir = file.path(Integrate_dir,"rds")
dir.create(rds_dir)




#Path to save Figures data
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)

#Add Annotation (Annovar + VEP) : see Mutagen script
#Annovar 
Annovar <- read.csv(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2",patient,"/Somatic/MuTect2/Annovar/",paste0(patient,".filtered.PASS.variants.annovar.hg38_multianno.csv")))
Annovar <- Annovar %>% mutate(mutation_id= paste0(":",sub("chr","",Annovar$Chr),":", Annovar$Start,":", Annovar$Ref, ":", Annovar$Alt),
                             Func_Annovar = paste0(Annovar$Func.refGene,"_",Annovar$ExonicFunc.refGene), 
                             Gene_Annovar = Gene.refGene) %>%                 
                        select(mutation_id,Func_Annovar,Gene_Annovar) 
#VEP 
VEP <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/",patient,"Somatic/MuTect2",paste0(patient,".filtered.PASS.vep_clean.vcf")), header = F )

colnames(VEP) <- c("Chr","Start","End","Ref","Alt","Func_VEP","Gene_VEP")
VEP <- VEP %>% mutate(mutation_id = paste0(":", sub("chr","",VEP$Chr),":", VEP$Start,":", VEP$Ref,":", VEP$Alt) ) %>% select(mutation_id,Func_VEP,Gene_VEP)

annotation <- left_join(Annovar,VEP, by = "mutation_id")

#check gene name lack and add VEP if necessary
sum(annotation$Gene_Annovar == "."  )
sum(annotation$Gene_VEP == "" )

annotation <- annotation %>%mutate(Gene_name=ifelse(Gene_Annovar == "." & Gene_VEP != "", Gene_VEP, Gene_Annovar ))

#check
sum(annotation$Gene_name == "." )    #fewer than Gene_Annovar 

#load the mutation detected in single cell
Combined_all = read.table(file.path(txt_files,paste0(patient,"_Combined_all.txt")))

#Specify the path


all_sample_rds_path <- "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/rds"
argu<-"clustered_0.8integrated_doublet_rm_filtered_all_patients_qc_cca"
seurat_integrated <-readRDS(file.path(all_sample_rds_path,paste0(argu,".rds")))

seurat_integrated$patient<- str_extract(seurat_integrated$orig.ident,"^[^[_-]]*")
seurat_integrated$treatment<-str_extract(seurat_integrated$orig.ident,"[^-]*$")

seurat_patient = subset(seurat_integrated, subset = patient == "OO77")

DefaultAssay(seurat_patient)<-'RNA'
options(future.globals.maxSize = 4294967296)  # 4 GB
seurat_patient <- SCTransform(seurat_patient,  vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_patient <- RunPCA(seurat_patient , npcs = 30, verbose = F)
seurat_patient <- RunUMAP(seurat_patient, reduction = "pca", dims = 1:30, verbose = F)
  

#Load the subclone information
seurat_no = readRDS(file.path(rds_dir, paste0(patient,"-no_Clones.RDS" )))


seurat_vitro = readRDS(file.path(rds_dir, paste0(patient,"-vitro_Clones.RDS")))


seurat_vivo =  readRDS(file.path(rds_dir, paste0(patient,"-vivo_Clones.RDS" )))

meta_all = rbind(seurat_no@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")], 
                seurat_vitro@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")]) %>% 
            rbind(seurat_vivo@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")])

meta_all = meta_all %>% mutate(Clone = ifelse( is.na(Clone), "UnI",as.character(Clone)))
table(meta_all$Clone)
#Assigne Clone to seurat_patient
seurat_patient@meta.data$CB_seurat = colnames(seurat_patient)
meta = seurat_patient@meta.data
meta = meta %>% left_join(meta_all[,c("CB_seurat","Clone")], by = "CB_seurat")
seurat_patient@meta.data = meta[,c("orig.ident", "treatment", "CB_seurat", "Clone")]
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat
table(seurat_patient@meta.data$Clone)


clones_list_all =   c( "2", "4","5","6", "7" ,"8", "UnI")
seurat_patient@meta.data = seurat_patient@meta.data %>% 
                               mutate( treatment = case_when(orig.ident == paste0(patient,"-no") ~ "Untreated",
                                                             orig.ident == paste0(patient,"-vitro") ~ "Ex vivo",
                                                             orig.ident  == paste0(patient,"-vivo") ~ "In vivo"))

seurat_patient@meta.data$treatment = factor(seurat_patient@meta.data$treatment , levels = c("Untreated", "Ex vivo", "In vivo"))
seurat_patient@meta.data$Clone = factor(seurat_patient@meta.data$Clone, levels = clones_list_all)


#save seurat
saveRDS(seurat_patient, file.path(rds_dir, paste0(patient,"_three_samples_raw_Clones.rds")))


#by Treatment
colors = c("cadetblue2", "tomato","plum")
plot <- DimPlot(seurat_patient, group.by= "treatment", cols = colors) + labs(title = patient)
plot 
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_DimPlot_all_cell.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')



#Save data of Plot 4A
data = plot$data 
all(rownames(data) == rownames(seurat_patient@meta.data))

#Combine with Clone information
data = cbind(data, Clone = seurat_patient@meta.data$Clone)

write.table(data, file.path(Figures_data, paste0("Figure_Sup_4A_initial_mapping_Clones_to_scRNA_OO77.txt")), sep = "\t", quote = F)





#Include the non-genotyped cells
# Define colors for clones and NA cells
clone_colors <- c(`2` = "lightgreen", `4` =  "lightpink",`5` = "dodgerblue4", 
  `6` = "mediumseagreen",`7` = "plum", `8` = "lightgoldenrod1", `UnI` = "gray92"  # Color for NA cells
)

# Extract UMAP embeddings and clone information
p <-  DimPlot(seurat_patient, reduction = "umap", group.by = "Clone") + theme_minimal()
plot_data <- p$data
# Convert `Clone` to factor and ensure `NA` is handled properly
plot_data$Clone <- factor(plot_data$Clone, levels =  clones_list_all)

# Plotting
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = Clone)) +
  # First layer: NA cells in background (gray)
  geom_point(data = subset(plot_data, Clone == "UnI"), size = 0.5) +
  geom_point(data = subset(plot_data, Clone == "2"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "4"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "5"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "6"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "7"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "8"), size = 0.75) +
  # Manual color scale for specific clone colors
  scale_color_manual(values = clone_colors)+
  guides(color = guide_legend(reverse = TRUE))+
   theme(
    legend.position = "right",
    legend.text = element_text(size = 14),        # Increase legend text size
    legend.title = element_text(size = 16),       # Increase legend title size (if used)
    legend.key.size = unit(10, "cm")             # Increase size of legend keys
  ) +
  theme_void()
  
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_by_Clone_all_cell.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')

#Only genotyped cells
#Extract cells with known clones and ncell > 3 
clones_list =  c( "2", "4","5","6", "7" ,"8")
seurat_patient_sub = subset(seurat_patient , subset = Clone %in% clones_list)
seurat_patient_sub@meta.data$Clone = factor(seurat_patient_sub@meta.data$Clone, levels = clones_list)

#by Treatment
colors = c("cadetblue2", "tomato","plum")
DimPlot(seurat_patient_sub, group.by= "treatment", cols = colors) + labs(title = patient)
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_DimPlot.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')
#by Clone
colors <- c(`2` = "lightgreen", `4` = "lightpink",`5` =  "dodgerblue4" , 
  `6` = "mediumseagreen",`7` = "plum", `8` = "lightgoldenrod1"
)

DimPlot(seurat_patient_sub, group.by= "Clone", cols = colors) + labs(title = patient) 
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_by_Clone.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')



#Mutation_matrix
combined_filtered = Combined_all %>% filter(mutated == TRUE, treeCLUSTER %in% clones_list)  %>% distinct(mutation_id, CB_seurat, .keep_all= T)
mutation_matrix <- as.matrix(table(combined_filtered$mutation_id, combined_filtered$CB_seurat))
dim(mutation_matrix)
#prepare the mutation meta and cell meta
#Rearrange mutation_matrix by mutation Clones and cell clones to plot heatmap

#mutation
Mutation_id_meta = data.frame(mutation_id = rownames(mutation_matrix)) %>% left_join(combined_filtered[,c("mutation_id","treeCLUSTER")] %>% distinct(mutation_id, .keep_all=TRUE), by ="mutation_id")
Mutation_id_meta$Clone= factor(Mutation_id_meta$treeCLUSTER , levels = clones_list ) 
Mutation_id_meta = Mutation_id_meta %>% arrange(Clone)
anno_mut = data.frame(Clone = Mutation_id_meta$Clone)
rownames(anno_mut) = Mutation_id_meta$mutation_id
dim(anno_mut)
#cell
Cell_meta =  data.frame( CB_seurat = colnames(mutation_matrix)) %>% left_join(seurat_patient_sub@meta.data[,c( "CB_seurat" , "Clone")], by = "CB_seurat")
Cell_meta = Cell_meta %>% filter(!(is.na(Clone)))
Cell_meta$Clone = factor(Cell_meta$Clone , levels = clones_list )
Cell_meta = Cell_meta %>% arrange(Clone)
anno_cell = data.frame(Clone = Cell_meta$Clone)
rownames(anno_cell) =  Cell_meta$CB_seurat
dim(anno_cell)
#rearrange mutation_matrix based on mutation meta and cell meta
mutation_matrix = mutation_matrix[Mutation_id_meta$mutation_id,Cell_meta$CB_seurat]


Clone_colors = list(Clone = colors)


#Mutations likely to be detected
row_sum  = rowSums(mutation_matrix)
top_mut = head(order(row_sum,decreasing = T),20)
top_mut_id = Mutation_id_meta[top_mut, ] 
top_mut_info = top_mut_id  %>% left_join(annotation, by = "mutation_id") 


#gap for row and col
table(anno_mut$Clone)
rows_gap = c(489,489+116, 489+ 116 + 535)
table(anno_cell$Clone)
cols_gap = c(2318 ,2318+ 521, 2318+  521 + 2914 )

pheatmap(mutation_matrix, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        legend_breaks = c(0, 1),
        gaps_row = rows_gap , gaps_col = cols_gap ,
        filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap.png"), width = 6,height=7,  units = "in", res = 600)

library(magick)

# Load the heatmap image
heatmap_img <- image_read(file.path(plot_dir, "mutation_matrix_3_samples_heatmap.png"), density = 600)
# Apply Gaussian Blur
blurred_img <- image_blur(heatmap_img, radius = 10, sigma = 1.5)

# Save the blurred image
image_write(blurred_img, path = file.path(plot_dir, "mutation_matrix_3_samples_heatmap_blurred.pdf"),  density = 600,format = "pdf")

#Save Figures data
write.table(mutation_matrix, file.path(Figures_data, "Figure_Sup_4B_mutation_matrix_initial_clone_assignment_OO77.txt"), sep = "\t", quote= F)
write.table(anno_mut,  file.path(Figures_data, "Figure_Sup_4B_mutation_matrix_initial_clone_assignment_OO77_annotation_mutation.txt"), sep = "\t", quote= F)
write.table(anno_cell,  file.path(Figures_data, "Figure_Sup_4B_mutation_matrix_initial_clone_assignment_OO77_annotation_cell.txt"), sep = "\t", quote= F)







#Cell mutation count
mut_count =  colSums(mutation_matrix)
table(mut_count )

# Highly detected mutations + PIK3CA
#mut_id = combined_filtered %>% filter(Gene_name == "PIK3CA")%>% pull(mutation_id) %>% unique()
muts = c(":12:71132813:G:C",":3:143016296:C:T" ,":3:179234297:A:G",":12:22065183:G:A")
names(muts) = c("TSPAN","U2SURP","PIK3CA","CMAS")
rownames(mutation_matrix) = ifelse(rownames(mutation_matrix) %in% muts, rownames(mutation_matrix) ,""  )

#heatmap with rownames
pheatmap(mutation_matrix, color = c("whitesmoke","darkred"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = T, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        legend_breaks = c(0, 1),
        legend_labels = c("no","yes"),
        gaps_row = rows_gap , gaps_col = cols_gap ,
        filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap_rownames.png"), width = 6,height=7)
 

#heatmap of gene expression 

DefaultAssay(seurat_patient_sub)<-'RNA'
options(future.globals.maxSize = 4294967296)  # 4 GB
seurat_patient_sub <- SCTransform(seurat_patient_sub,  vst.flavor = "v1",
                                 variable.features.n = nrow(seurat_patient_sub),
                                 verbose = FALSE)
scaled_data = as.data.frame(seurat_patient_sub$SCT$scale.data)
data = scaled_data %>% mutate(Gene_name = rownames(scaled_data)) 
mut_df = data.frame(mutation_id = rownames(mutation_matrix))

#Modify gene_name in annotation so that integenic mutations take the gene_name for the gene right before
Modify_Gene_name <- function(Gene_name) {
  # Check if Gene_name contains a semicolon
  if (grepl(";", Gene_name)) {
    # Split Gene_name by ";"
    gene_list <- str_split(Gene_name, ";")[[1]]
    # Find the first gene in gene_list that is in rownames(seurat_patient_sub)
    matching_gene <- gene_list[gene_list %in% rownames(seurat_patient_sub)][1]
    # Return the matching gene if found, otherwise return the 1st gene
    return(ifelse(!is.na(matching_gene), matching_gene, gene_list[1]))
  } else {
    # If there's no semicolon, return Gene_name as it is
    return(Gene_name)
  }
}


annotation = annotation %>% rowwise() %>% mutate(Gene_name_new = Modify_Gene_name(Gene_name) )

mut_df = mut_df %>% left_join(annotation[,c("mutation_id","Gene_name_new")], by = "mutation_id") %>% rename("Gene_name" = "Gene_name_new")

mut_scale_df = mut_df %>% left_join(data, by = "Gene_name")

mut_scale_mat = as.matrix(mut_scale_df[,3: ncol(mut_scale_df )])
#reorder cell CB according to mutation_matrix 
order = match(colnames(mutation_matrix),colnames(mut_scale_mat))
mut_scale_mat = mut_scale_mat[,order]
all(colnames(mut_scale_mat) == colnames(mutation_matrix))
rownames(mut_scale_mat) = mut_scale_df$mutation_id

pheatmap(mut_scale_mat, colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[3:9])(7),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut , breaks = c(-6,-4,-2,0,2,4,6,8),
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        gaps_row = rows_gap , gaps_col = cols_gap ,
        filename = file.path(plot_dir,"mutation_gene_expression_3_samples_heatmap.pdf"), width = 6,height=7)



#cluster row
#drop NA value
na_rows = which(apply(mut_scale_df,1, function(row){sum(is.na(row))>0}))
length(na_rows)
non_na_rows = setdiff(1:nrow(mut_scale_df),na_rows )
non_na_rows_id = rownames(mut_scale_mat)[non_na_rows]
mut_scale_mat_nonNA = mut_scale_mat[non_na_rows_id,]
anno_mut_nonNA = anno_mut[non_na_rows_id,"Clone", drop = FALSE]

pheatmap(mut_scale_mat_nonNA, colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[3:9])(7),
        cluster_rows = T, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut_nonNA , breaks = c(-6,-4,-2,0,2,4,6,8),
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        filename = file.path(plot_dir,"mutation_gene_expression_3_samples_heatmap_cluster_row.png"), width = 6,height=7)
 
#Clone 2
combined_filtered_distinct = combined_filtered %>% distinct(CB_seurat,.keep_all = T)
Cell_meta_c = Cell_meta %>% filter(Clone == 2) %>% left_join(combined_filtered_distinct[,c("CB_seurat","sample")], by ="CB_seurat")
mutation_matrix_c = mutation_matrix[,Cell_meta_c$CB_seurat]

anno_cell_sample = data.frame(sample = Cell_meta_c$sample) %>% 
                    mutate(sample = case_when(sample == "OO77-no" ~ "Untreated",
                                              sample == "OO77-vitro" ~ 'Ex vivo',
                                              sample == "OO77-vivo" ~ 'In vivo'))
anno_cell_sample$sample = factor(anno_cell_sample$sample , levels = c("Untreated", "Ex vivo", "In vivo"))                                    
rownames(anno_cell_sample) = Cell_meta_c$CB_seurat

table(anno_cell_sample$sample)
anno_colors = list(
  # Column annotation (for 'sample' in anno_cell_sample)
  sample = c("Untreated" = "tomato", 
             "Ex vivo" = "deepskyblue2",  
             "In vivo" = "khaki3"),
  Clone = Clone_colors$Clone)

pheatmap(mutation_matrix_c[,1:337], color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = head(anno_cell_sample, 337), annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_2_337cells.png"), width = 8 ,height = 8)



muts = c(":12:71132813:G:C",":3:143016296:C:T" ,":3:179234297:A:G",":12:22065183:G:A", ":1:1218852:T:G", ":22:18092823:T:G")

rownames(mutation_matrix_c) = ifelse(rownames(mutation_matrix_c) %in% muts, rownames(mutation_matrix_c) ,""  )
pheatmap(mutation_matrix_c[,1:337], color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = T, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = head(anno_cell_sample, 337), annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_2_337cells_rownames.png"), width = 10 ,height= 8)



#HOW depleted THE CELL CLONE of mutations of clones from separated branches? ( how pure the clone defined is)
cell_anno = anno_cell
cell_anno = cell_anno %>% mutate(diff_branch_clones = case_when(Clone == 2 ~ paste(4,5,6,7,8),
                                                                Clone == 4 ~ paste(2,5,6,7,8),
                                                                Clone == 5 ~ paste(2,4,6,7,8),
                                                                Clone %in% c(6,7,8) ~ paste(2,4,5))
                                                                )
cell_anno$Clone = as.character(cell_anno$Clone)
rownames(mutation_matrix) <- Mutation_id_meta$mutation_id



#Function to calculate Cell clone purity
Purity_calculate <- function(Clone_current){
diff_branch_clones = cell_anno %>% filter(Clone == Clone_current) %>% pull(diff_branch_clones)%>% unique() %>% str_split(,pattern = " ") %>% unlist() 
CB_clone = cell_anno %>% filter(Clone  == Clone_current)%>% rownames()
#mutations in other  branches
diff_branch_muts = anno_mut %>% filter(Clone %in% diff_branch_clones ) %>% rownames()

mutation_matrix_clone_diff_branch = mutation_matrix[diff_branch_muts ,CB_clone]
Nb_contaminated_cells = sum(colSums(mutation_matrix_clone_diff_branch) > 0 )
Purity = 1 - Nb_contaminated_cells/length(CB_clone)
df = data.frame(Clone = Clone_current, Purity = Purity , n_cell = length(CB_clone) )
df
}

Purity = purrr::map_dfr( clones_list,Purity_calculate )
Purity$Clone = factor(Purity$Clone , levels = names(clone_colors))


# Calculate the minimum value of Purity
min_purity <- min(Purity$Purity, na.rm = TRUE)

# Create the bar plot with a horizontal line and display the minimum value
Purity %>% ggplot(aes(x = Clone, y = Purity, fill = Clone)) + 
  scale_fill_manual(values = clone_colors) +
  scale_y_continuous(limits = c(0,1.15), breaks = seq(0,1,0.2))+
  geom_bar(stat = "identity", width = 0.5) + 
  xlab("Cell clone")+
  geom_text(aes(label = n_cell), vjust = -0.5, size = 11, color = "blue" ) +  # Add n_cell values as labels above bars
  geom_hline(yintercept = min_purity, linetype = "dashed", color = "red", size = 1) +  # Add hline
  annotate("text", x = 5, y = min_purity, label = paste("Min:", round(min_purity, 2)),
           vjust = 2, hjust = 1.1, color = "red", size = 11) +  # Add the text annotation
  theme_classic() +
  theme(legend.title =  element_text(size = 32, color = "black"),
        legend.text =  element_text(size = 32, color = "black"),
        axis.text.y.left = element_text(size = 32, color = "black"),
        axis.title.y.left = element_text(size = 32, color = "black"),
        axis.title.x.bottom = element_text(size = 32, color = "black"),
        axis.text.x.bottom = element_text(size = 32, color = "black"))

# Save the plot
ggsave(file.path(plot_dir, "Purity.pdf"), width = 11, height = 5, scale = 1)

#save Figure data
write.table(Purity, file.path(Figures_data, "Figure_Sup_4C_Purity_cell_clone_initial_clone_assignment_OO77.txt"), sep = "\t", quote = F)

#Proportion of each clone
Clone_sample = as.data.frame(table(seurat_patient_sub$orig.ident, seurat_patient_sub$Clone))
colnames(Clone_sample) = c("sample","Clone","count")


Clone_sample <- Clone_sample %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100) # Convert counts to percentages

# Plot proportions as percentages
# Extract the unique order of clones based on the new order
ordered_clones <- Clone_sample %>%
  filter(sample == "OO77-no") %>%  # Filter for the "Untreated" sample
  arrange(desc(proportion)) %>%      # Arrange by proportion in descending order
  pull(Clone)


# Plot proportions as percentages
Clone_sample %>%
  ggplot(aes(x = factor(sample, levels = rev(levels(Clone_sample$sample))), y = proportion, fill = factor(Clone, levels = rev(ordered_clones)))) +
  geom_bar(stat = "identity", width = 0.8) + 
  scale_fill_manual(values = clone_colors ) + 
  coord_flip() +
  ylab("Clone percentage (%)") +
  scale_x_discrete(labels= c("OO77-no" = "Untreated", "OO77-vitro" = "Ex vivo", "OO77-vivo" = "In vivo")) +
  theme_minimal() +
  theme(legend.position = "top",
     legend.direction = "horizontal",   
    legend.title = element_blank(), 
    legend.text = element_text(size = 30),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_text(size= 30),
    axis.text.x.bottom = element_text(size= 30),
    axis.title.x.bottom = element_text(size= 30)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE) )

# Save the plot
ggsave(file.path(plot_dir, "OO77_clone_proportion.pdf"), width = 10, height = 6)





## CONTAMINATE

#Contaminated cells of each clone (fraction of cells contain mutation from other clone)
Contaminate_cell_fraction = function(clone){
   CB_clone = anno_cell %>% filter(Clone == clone) %>% rownames()
   mut_mat_CB_clone = mutation_matrix[,CB_clone]
   clone_others = clones_list
   #calculate contaminate clone
   conta_clone_other = function(clone_other){
        mut_clone_other = anno_mut %>% filter(Clone == clone_other) %>% rownames()
        mut_mat_CB_clone_mut_other = mut_mat_CB_clone[mut_clone_other,]
        Nb_contaminated_cells = sum(colSums(mut_mat_CB_clone_mut_other) > 0 )
        fraction = 100*Nb_contaminated_cells/length(CB_clone)
        df = data.frame(Clone = clone, Clone_other = clone_other, Fraction = fraction)
        df
}
  conta_df = purrr::map_dfr(clone_others, conta_clone_other)
  conta_df
}


Contaminate_cell_each_clone = purrr::map_dfr(clones_list,Contaminate_cell_fraction)
Contaminate_cell_each_clone$Clone = factor(Contaminate_cell_each_clone$Clone , levels = clones_list)
Contaminate_cell_each_clone$Clone_other = factor(Contaminate_cell_each_clone$Clone_other , levels = clones_list)

ggplot(Contaminate_cell_each_clone , aes(x = Clone, y = Fraction, fill = as.factor(Clone_other))) +
  geom_bar(stat = "identity", position = "dodge") +   # Stacked bar for each clone
  scale_y_continuous(breaks = seq(0,100,50))+
  labs(x = "Cell Clone", y = "Fraction of cell", fill = "Clone") +
  theme_classic() +
  scale_fill_manual(values= clone_colors)    +      
  theme(plot.title =   element_text(size = 32, color = "black", hjust = 0.5),
        legend.title =  element_text(size = 32, color = "black"),
        legend.text =  element_text(size = 32, color = "black"),
        panel.grid.major.y = element_line(color = "grey80", size = 0.3),
        panel.grid.minor.y = element_line(color = "grey80", size = 0.3),
        panel.grid.major.x = element_line(color = "grey80", size = 0.3),
        axis.text.y.left = element_text(size = 32, color = "black"),
        axis.title.y.left = element_text(size = 32, color = "black"),
        axis.title.x.bottom = element_text(size = 32, color = "black"),
        axis.text.x.bottom = element_text(size = 32, color = "black"),
        strip.text = element_text(size = 3, color = "black"),
        strip.background = element_blank())+
  facet_wrap(vars(Clone), ncol = 1)

ggsave(file.path(plot_dir, paste0(patient,"_fraction_cell_in_mutation_clone.pdf")), width = 11, height = 10, scale = 1)

#save Figure data
colnames(Contaminate_cell_each_clone ) = c("Cell_clone", "Mutation_clone", "Fraction")
write.table(Purity, file.path(Figures_data, "Figure_Sup_4C_Fraction_mutation_clone_in_cell_clone_initial_clone_assignment_OO77.txt"), sep = "\t", quote = F)


## CONTAMINATE MUT
#do for mutation clone 
clones_list_al2 =  names(table(anno_mut$Clone))[table(anno_mut$Clone) >= 2]

Contaminate_mut_fraction = function(clone,clones_cell_list, mutation_matrix , anno_mut, anno_cell){
   mut_clone = anno_mut %>% filter(Clone == clone) %>% rownames()
   mut_mat_mut_clone = mutation_matrix[mut_clone,]
   #calculate contaminate clone
   detected_clone_other = function(clone_other){
        CB_clone_other = anno_cell %>% filter(Clone == clone_other) %>% rownames()
        mut_mat_mut_clone_CB_other = mut_mat_mut_clone[,CB_clone_other]
        Nb_contaminated_mut = sum(rowSums(mut_mat_mut_clone_CB_other) > 0 )
        fraction = 100*  Nb_contaminated_mut /length( mut_clone)
        df = data.frame(Clone = clone, Clone_other = clone_other, Fraction = fraction)
        df
}
  conta_df = purrr::map_dfr(clones_cell_list,  detected_clone_other )
  conta_df
}


Contaminate_mut_each_clone = purrr::map_dfr(clones_list_al3 ,~Contaminate_mut_fraction(.x, clones_cell_list = clones_list_al3 , mutation_matrix = mutation_matrix, anno_mut = anno_mut, anno_cell=anno_cell))
Contaminate_mut_each_clone$Clone = factor(Contaminate_mut_each_clone$Clone , levels = clones_list)
Contaminate_mut_each_clone$Clone_other = factor(Contaminate_mut_each_clone$Clone_other , levels = clones_list)

ggplot(Contaminate_mut_each_clone , aes(x = Clone, y = Fraction, fill = Clone_other)) +
  geom_bar(stat = "identity", position = "dodge") +   # Stacked bar for each clone
  labs(x = "Mutation Clone", y = "Fraction", fill = "Cell Clone") +
  theme_minimal() +
  scale_fill_manual(values= colors)    +      
  ggtitle("Fraction of mutations detected in scRNA-seq cell clones")

ggsave(file.path(plot_dir, paste0(patient,"_clone_contaminated_mut.pdf")), width = 1500, height = 500, dpi = 150, units = 'px')


###### BOOSTRAP ####
# Function to shuffle cell identities and recalculate contamination
plot_dir_boostrap = file.path(plot_dir,"boostrap_mut_clone")
dir.create(plot_dir_boostrap,recursive = T)


Bootstrap_contaminate_mut_fraction <- function(clone_mut,clones_cell_list,mutation_matrix,anno_mut, anno_cell,plot_dir_boostrap, n_bootstrap = 1000) {
  mut_clone = anno_mut %>% filter(Clone == clone_mut) %>% rownames()
  mut_mat_mut_clone = mutation_matrix[mut_clone, ]
  
  # Create an empty data frame to store p-values
  p_value_results <- data.frame(
    clone_mut = character(),
    clone_cell = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  Each_clone_cell_function <- function(clone_cell, n_bootstrap = n_bootstrap) {
    # Observed contamination
    observed_value <- Contaminate_mut_fraction(clone_mut,clones_cell_list, mutation_matrix = mutation_matrix, anno_mut = anno_mut, anno_cell = anno_cell) %>% 
      filter(Clone_other == clone_cell) %>% 
      pull(Fraction)
    
    # Initialize vector to store bootstrap contamination fractions
    bootstrap_results <- vector("numeric", length = n_bootstrap)
    
    for (i in 1:n_bootstrap) {
      shuffled_anno_cell <- anno_cell
      cells_to_shuffle <- shuffled_anno_cell %>% filter(Clone != clone_mut)
      n_clone_cell <- anno_cell %>% filter(Clone == clone_cell) %>% nrow()
       #shuffle the cell clone label
      set.seed(i)
      shuffled_CB <- sample(rownames(cells_to_shuffle), size = min(n_clone_cell,nrow(cells_to_shuffle)))
      shuffled_subset <- cells_to_shuffle[shuffled_CB, , drop = FALSE]
      CB_clone_cell <- rownames(shuffled_subset)
      
      mut_mat_clone_other <- mut_mat_mut_clone[, CB_clone_cell]
      num_contaminated_mut <- sum(rowSums(mut_mat_clone_other) > 0)
      contamination_fraction <- 100 * num_contaminated_mut / length(mut_clone)
      
      bootstrap_results[i] <- contamination_fraction
    }
    
    p_value <- mean(bootstrap_results >= observed_value)
    
    # Append the results to the data frame
    p_value_results <<- rbind(
      p_value_results, 
      data.frame(clone_mut = clone_mut, clone_cell = clone_cell, p_value = p_value)
    )
    
    # Generate plot (optional)
    p <- ggplot(data.frame(bootstrap_results), aes(x = bootstrap_results)) +
      geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
      annotate("text", x = 75, y = 1000, label = paste("p =", format.pval(p_value, digits = 2)), size = 4) +
      xlim(-10, 100) +
      ylim(0, 1000) +
      geom_vline(aes(xintercept = observed_value), color = "red", linetype = "dashed", size = 1) +
      labs(title = paste("Muts of", clone_mut, "in cell of", clone_cell),
           x = "Fraction", y = "Frequency") +
      theme_minimal()
    
    return(p)
  }
  
  # Create list of plots for each clone in clone_cell_list
  p_list <- purrr::map(clones_cell_list, ~ Each_clone_cell_function(.x, n_bootstrap = n_bootstrap))
  
  # Combine plots using ggarrange
  do.call(ggarrange, c(p_list, nrow = 1))
  
  # Save the final arranged plot
  ggsave(file.path(plot_dir_boostrap, paste0(patient, "_Bootstrap_clone_contaminated_m", clone_mut, ".pdf")), width = 22, height = 2.5, dpi = 300)
  
  # Return the p-value data frame for further use
  return(p_value_results)
}


# Apply the function for each clone_mut in clones_list
pvalue_results = purrr::map_dfr(clones_list_al2,~ Bootstrap_contaminate_mut_fraction(.x, clones_cell_list = clones_list_al2, mutation_matrix = mutation_matrix,anno_mut= anno_mut, anno_cell = anno_cell, plot_dir_boostrap= plot_dir_boostrap))
pvalue_results$patient = "OO77"
#save pvalue_results
write.table(pvalue_results, file.path(txt_files,paste0(patient,"_bootstrap_pvalues.txt")), sep = "\t", quote = F)

#Load all pvalues of 3 patients
OO77 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO77/tree30_change_name/txt_files/OO77_bootstrap_pvalues.txt")
OO99 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO99/tree30_change_name/txt_files/OO99_bootstrap_pvalues.txt")
OO100 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO100/tree30_change_name/txt_files/OO100_bootstrap_pvalues_remove_HLA.txt")


All_pvalue = rbind(OO77, OO99, OO100)
All_pvalue$padjust =  p.adjust(All_pvalue $p_value, method = "BH")

pvalue_results = All_pvalue  %>% filter(patient == "OO77")
pvalue_result_wider = pvalue_results %>% pivot_wider(id_cols = c("clone_mut"),names_from = "clone_cell", values_from = "padjust") %>% as.data.frame()

rownames(pvalue_result_wider) = pvalue_result_wider$clone_mut 
pvalue_result_wider = pvalue_result_wider[,2:ncol(pvalue_result_wider)]

pheatmap(pvalue_result_wider , colorRampPalette(brewer.pal(n = 11, name ="RdBu")[c(2,4,5,8)])(4), 
         breaks = c(0,0.005,0.01,0.05,0.1), legend_breaks = c(0,0.005,0.01,0.05,0.1),
         show_colnames = T, show_rownames = T,fontsize = 8,
         cluster_cols = F,cluster_rows = F,
         gaps_col = c(1,2,3), gaps_row = c(1,2,3),
         filename = file.path(plot_dir_boostrap,"pvalue_all.pdf"),
          width = 4, height = 4)

#Merge clone 4 to 5
Combined_all = Combined_all %>% mutate( treeCLUSTER_merged = ifelse(treeCLUSTER == 4,5, treeCLUSTER)) 



# CALCULATE INTRA CLONE DISTANCE

plot_Intra_Clone_score_dir = file.path(plot_dir,"Intra_clone_Score")
dir.create(plot_Intra_Clone_score_dir)

#Use 30 PC
pca_embeddings <- Embeddings(seurat_patient, reduction = "pca")
pca_30 <- t(as.matrix(pca_embeddings[,1:30]))

#Function to calculate the intr-clone pairwise distance
Distance_2_muts <- function(set_2_muts){
mut1 = set_2_muts[1]
mut2 = set_2_muts[2]
CB_mut1 = Combined_all %>% filter((mutation_id == mut1) & (mutated == T)) %>% pull(CB_seurat)
CB_mut2 = Combined_all %>% filter((mutation_id == mut2) & (mutated == T)) %>% pull(CB_seurat)
meta = seurat_patient@meta.data 
meta$CB = rownames(meta)
meta = meta %>% mutate(mut1 = ifelse(CB %in% CB_mut1,TRUE, FALSE),
                       mut2 = ifelse(CB %in% CB_mut2, TRUE, FALSE))

ALT_mut1 <- meta %>% filter(mut1 == TRUE) %>% pull(CB) 
ALT_mut2 <-  meta %>% filter(mut2 == TRUE) %>% pull(CB) 


# Centroid method
pca_ALT_mut1 <- pca_30[ , ALT_mut1]
pca_ALT_mut2 <- pca_30[, ALT_mut2]
centroid_ALT_mut1 <- rowMeans(pca_ALT_mut1)
centroid_ALT_mut2 <- rowMeans(pca_ALT_mut2)
centroids_distance <- 1- cor(centroid_ALT_mut1, centroid_ALT_mut2, method="spearman")
df = data.frame(mut1 = mut1, mut2 = mut2, dist = centroids_distance)
#Complete distance method
# ALT_both <- c(ALT_mut1,ALT_mut2)
# spearman_dist_mut <- 1-cor(pca_30[,ALT_both], method = "spearman")
# dist = mean(spearman_dist_mut[ALT_mut1,ALT_mut2])
# df = data.frame(mut1 = mut1, mut2 = mut2, dist = dist)
df
}

#Function to filter mut with ALT >= threshold 
Filter_mut <- function(mut, threshold){
CB_pos = Combined_all %>% filter((mutation_id == mut) & (mutated == T)) %>% pull(CB_seurat)
CB_neg = Combined_all %>% filter((mutation_id == mut) & (mutated == F)) %>% pull(CB_seurat)
meta = seurat_patient@meta.data 
meta$CB = rownames(meta)
meta = meta %>% mutate(mut = case_when( CB %in% CB_pos ~ "ALT",
                                        CB %in% CB_neg ~ "REF",
                                        TRUE ~ "UNK"))

n_ALT <- meta %>% filter(mut == "ALT") %>% pull(CB) %>% length()
n_REF <- meta %>% filter(mut == "REF") %>% pull(CB) %>% length()
return = FALSE
if (n_ALT >= threshold) {return= TRUE} 
return
}

#Calculate Score
Intra_Clone_Dist_score <- function(clone, threshold){
  mutations = Combined_all %>% filter( treeCLUSTER_merged == clone) %>% pull(mutation_id) %>% unique()

  keep = purrr::map(mutations,~Filter_mut(.x, threshold = threshold)) %>% unlist()

  mutations = mutations[keep]
  Score = data.frame(mut1 = "", mut2 = "", dist = 0)
  if (length(mutations) > 1){
  sets_2_muts_mat <- combn(mutations, 2)
  sets_2_muts  <- split(sets_2_muts_mat , col(sets_2_muts_mat ))
  Score = purrr::map_dfr(sets_2_muts , Distance_2_muts)
  } 
  Score$Clone = clone
  Score$n_muts = length(mutations)
  Score
}

clones = unique(Combined_all$treeCLUSTER_merged)

Threshold_test_intra_clone <- function(threshold){
df = purrr::map_dfr(clones,~Intra_Clone_Dist_score(.x, threshold = threshold ))
df$Clone = factor(as.character(df$Clone), levels = clones)
df$threshold = threshold
df
}

thresholds = 6
All_thresholds_intra_clone = purrr::map_dfr(thresholds,Threshold_test_intra_clone)

#Plot histogram of these distances
for (clone in clones){
Distance_Clone_df = All_thresholds_intra_clone %>% filter(Clone == clone) 
nmuts =  Distance_Clone_df %>% pull(n_muts) %>% unique()

#test unimodality
Dists  = Distance_Clone_df %>% pull(dist)
dip_test <- dip.test(Dists)
p_value = dip_test$p.value

ggplot(Distance_Clone_df , aes(x = dist)) +
  geom_histogram(bins = 25, fill = "gray", color = "black") +
  xlim(0, 1.7) +
 ylim(0, 1700) +
  labs(
    title = paste0("Clone ", clone, " (", nmuts, " mutations)"),
    x = "PIC distance",
    y = "Frequency"
  ) +
annotate("text", x = 1.4, y = 1600, label = paste("p =", signif(p_value, 4)), size = 8.5)+
  theme_classic() +
  theme(
    plot.title = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 25, color = "black"),
    axis.title.y = element_text(size = 25, color = "black"),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
ggsave(file.path(plot_Intra_Clone_score_dir,paste0("Clone_", clone,".pdf")), width= 5.5, height = 4.5)
}

#Save Figure data
data = All_thresholds_intra_clone %>% select(Clone, n_muts,mut1, mut2, dist)
write.table(data, file.path(Figures_data, "Figure_Sup_5A_Pairwise_Intra_Clone_Distance_OO77.txt"),sep = "\t", quote= F) 


#Separate clone 5
Centroid_mut <- function(mut){
  CB_mut = Combined_all %>% filter((mutation_id == mut) & (mutated == T)) %>% pull(CB_seurat)
  meta = seurat_patient@meta.data 
  meta$CB = rownames(meta)
  meta = meta %>% mutate(mut = ifelse(CB %in% CB_mut,TRUE, FALSE))
  
  #muatated cells
  ALT_mut <- meta %>% filter(mut == TRUE) %>% pull(CB) 

  # Centroid 
  pca_ALT_mut <- as.data.frame(pca_30[ , ALT_mut])
  if (nrow(pca_ALT_mut) == 1){
    df = as.data.frame( pca_ALT_mut)
    }else{
  centroid_ALT_mut <- rowMeans(pca_ALT_mut)
  df = as.data.frame(t(centroid_ALT_mut))
    }
  df$mutation_id = mut
  df
}

mutations = Combined_all %>% filter(mutated == TRUE, treeCLUSTER_merged == 5) %>% pull(mutation_id) %>% unique()

clone_centroids = purrr::map_dfr(mutations, Centroid_mut)

#Cluster mutations
set.seed(123) # Setting seed for reproducibility
kmeans_result <- kmeans(clone_centroids[,1:30], centers = 2)

# Add the cluster assignments to the data frame
clone_centroids$cluster <- kmeans_result$cluster
clone_centroids <- clone_centroids %>% mutate(cluster = ifelse(cluster == 1, "5b","5a"))

#Save the result of separation
write.table(clone_centroids[,31:32],file.path(txt_files,"OO77_Cluster_5_merged_separation.txt"),quote=F, sep = "\t")




#Plot each clone
plot_dir_sub = file.path(plot_dir, "each_clone")
dir.create(plot_dir_sub)

Combined_all_CB = Combined_all %>% filter(mutated == TRUE) %>% group_by(CB_seurat) %>% summarize(clone = treeCLUSTER_merged)
mutated_CB = as.data.frame(table(Combined_all_CB$CB_seurat,Combined_all_CB$clone)) 
colnames(mutated_CB) = c("CB_seurat","Clone","Freq")
mutated_CB = mutated_CB %>% pivot_wider(names_from = "Clone", values_from = "Freq")
head(mutated_CB)
#Clone_1
mutated_CB  <- mutated_CB %>% 
  mutate(Clone_1 = case_when(
    `1` > 0 ~ "Detected",
    TRUE ~ "Undetected"
  ))

meta = seurat_patient@meta.data
meta$CB_seurat = rownames(meta) 
meta = meta %>% left_join(mutated_CB[,c("CB_seurat","Clone_1")], by = "CB_seurat")
meta$Clone_1[(meta$Clone_1 != "Detected") | (is.na(meta$Clone_1 ))] =    "Undetected"
meta$Clone_1 = factor(meta$Clone_1 , levels = c("Detected",  "Undetected"  ))
seurat_patient@meta.data = meta
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat

#seurat_patient_Clone1 <- subset(seurat_patient, Clone_1 == c("Clone 1"))

colors = c("cadetblue2","gray92")
DimPlot(seurat_patient, reduction = "umap", group.by = "Clone_1" , cols = colors, raster = T, pt.size = 2.5) + 
       labs(title = "Clone 1") +
       theme_void()+
       theme(legend.position = "bottom",
             plot.title = element_text(size = 25),
             legend.text = element_text(size = 25) )  
ggsave(paste0("DimPlot_Clone_1.pdf"),path = plot_dir_sub,width = 4,height = 4.5)


#Clone_2
mutated_CB  <- mutated_CB %>% 
  mutate(Clone_2 = case_when(
    `2` > 0 ~ "Detected",
    TRUE ~ "Undetected"
  ))

meta = seurat_patient@meta.data
meta$CB_seurat = rownames(meta) 
meta = meta %>% left_join(mutated_CB[,c("CB_seurat","Clone_2")], by = "CB_seurat")
meta$Clone_2[(meta$Clone_2 != "Detected") | (is.na(meta$Clone_2 ))] =    "Undetected"
meta$Clone_2 = factor(meta$Clone_2 , levels = c("Detected",  "Undetected"))
seurat_patient@meta.data = meta
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat
#seurat_patient_Clone2 <- subset(seurat_patient, Clone_2 == c("Clone 2"))

colors = c("lightgreen", "gray92")
DimPlot(seurat_patient, reduction = "umap", group.by = "Clone_2" , cols = colors,   raster = T, pt.size = 2.5) +
  labs(title = "Clone 2")+
   theme_void()+
    theme(legend.position = "bottom",
         plot.title = element_text(size = 25),
         legend.text = element_text(size = 25) )  
ggsave(paste0("DimPlot_Clone_2.pdf"),path = plot_dir_sub,width = 4,height = 4.5)


#Clone_5
mutated_CB  <- mutated_CB %>% 
  mutate(Clone_5 = case_when(
    `5` > 0 ~ "Detected",
    TRUE ~ "Undetected"
  ))

meta = seurat_patient@meta.data
meta$CB_seurat = rownames(meta) 
meta = meta %>% left_join(mutated_CB[,c("CB_seurat","Clone_5")], by = "CB_seurat")
meta$Clone_5[(meta$Clone_5 != "Detected") | (is.na(meta$Clone_5 ))] =    "Undetected"
meta$Clone_5 = factor(meta$Clone_5 , levels = c("Detected",  "Undetected"))
seurat_patient@meta.data = meta
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat


colors = c("dodgerblue2", "gray92")
DimPlot(seurat_patient, reduction = "umap", group.by = "Clone_5" , cols = colors,   raster = T, pt.size = 2.5) + 
     labs(title = "Clone 5") +
      theme_void()+
      theme(legend.position = "bottom",
            plot.title = element_text(size = 25),
            legend.text = element_text(size = 25) )  
ggsave(paste0("DimPlot_Clone_5.pdf"),path = plot_dir_sub,width = 4,height = 4.5)


#Clone_6
mutated_CB  <- mutated_CB %>% 
  mutate(Clone_6 = case_when(
    `6` > 0 ~ "Detected",
    TRUE ~ "Undetected"
  ))

meta = seurat_patient@meta.data
meta$CB_seurat = rownames(meta) 
meta = meta %>% left_join(mutated_CB[,c("CB_seurat","Clone_6")], by = "CB_seurat")
meta$Clone_6[(meta$Clone_6 != "Detected") | (is.na(meta$Clone_6 ))] =    "Undetected"
meta$Clone_6 = factor(meta$Clone_6 , levels = c("Detected",  "Undetected"))
seurat_patient@meta.data = meta
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat


colors = c("mediumseagreen", "gray92")
DimPlot(seurat_patient, reduction = "umap", group.by = "Clone_6" , cols = colors,   raster = T, pt.size = 2.5) +
        labs(title = "Clone 6") +
        theme_void()+
         theme(legend.position = "bottom",
            plot.title = element_text(size = 25),
            legend.text = element_text(size = 25) )  
ggsave(paste0("DimPlot_Clone_6.pdf"),path = plot_dir_sub,width = 4,height = 4.5)


#save Figures data
umap_coords <- Embeddings(seurat_patient, reduction = "umap")

all(rownames(umap_coords) == rownames(seurat_patient@meta.data))
data = cbind(umap_coords, seurat_patient@meta.data)

head(data)

write.table(data, file.path(Figures_data, "Figure_Sup_5A_UMAP_clone_detected_OO77.txt"), sep = "\t", quote= F)
