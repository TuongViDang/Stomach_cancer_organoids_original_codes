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
patient="OO100"


main_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA",patient,"tree30_change_name")
plot_dir = file.path(main_dir,"plots","Map_validate_3_samples_filtered")
dir.create(plot_dir, recursive = T)
txt_files = file.path(main_dir,"txt_files")
dir.create(txt_files,recursive=T)
rds_dir = file.path(main_dir,"rds")
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
Combined_all = read.table(file.path(txt_files,paste0(patient,"_Combined_all_filtered.txt")))

#Specify the path


all_sample_rds_path <- "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/rds"
argu<-"clustered_0.8integrated_doublet_rm_filtered_all_patients_qc_cca"
seurat_integrated <-readRDS(file.path(all_sample_rds_path,paste0(argu,".rds")))

seurat_integrated$patient<- str_extract(seurat_integrated$orig.ident,"^[^[_-]]*")
seurat_integrated$treatment<-str_extract(seurat_integrated$orig.ident,"[^-]*$")


seurat_patient = subset(seurat_integrated, subset = patient == "OO100")

DefaultAssay(seurat_patient)<-'RNA'
options(future.globals.maxSize = 4294967296)  # 4 GB
seurat_patient <- SCTransform(seurat_patient,  vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_patient <- RunPCA(seurat_patient , npcs = 30, verbose = F)
seurat_patient <- RunUMAP(seurat_patient, reduction = "pca", dims = 1:30, verbose = F)

#Load the subclone information
seurat_no = readRDS(file.path(rds_dir, paste0(patient,"-no_Clones_filtered.RDS" )))


seurat_vitro = readRDS(file.path(rds_dir, paste0(patient,"-vitro_Clones_filtered.RDS")))


seurat_vivo =  readRDS(file.path(rds_dir, paste0(patient,"-vivo_Clones_filtered.RDS" )))

meta_all = rbind(seurat_no@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")], 
                seurat_vitro@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")]) %>% 
            rbind(seurat_vivo@meta.data[,c("CB_seurat","Clone","orig.ident", "Phase", "seurat_clusters")])

meta_all = meta_all %>% mutate(Clone = ifelse( is.na(Clone), "UnI",as.character(Clone)))
table(meta_all$Clone)

#Assigne Clone to seurat_patient
seurat_patient@meta.data$CB_seurat = colnames(seurat_patient)
meta = seurat_patient@meta.data
meta = meta %>% left_join(meta_all[,c("CB_seurat","Clone")], by = "CB_seurat")
seurat_patient@meta.data = meta[,c("orig.ident", "treatment", "CB_seurat", "Clone", "seurat_clusters")]
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat
table(seurat_patient@meta.data$Clone)


clones_list_all =   c( "2", "3","4","5", "7" ,"8","9","10","11","12", "UnI")
seurat_patient@meta.data = seurat_patient@meta.data %>% 
                               mutate( treatment = case_when(orig.ident == paste0(patient,"-no") ~ "Untreated",
                                                             orig.ident == paste0(patient,"-vitro") ~ "Ex vivo",
                                                             orig.ident  == paste0(patient,"-vivo") ~ "In vivo"))

seurat_patient@meta.data$treatment = factor(seurat_patient@meta.data$treatment , levels = c("Untreated", "Ex vivo", "In vivo"))
seurat_patient@meta.data$Clone = factor(seurat_patient@meta.data$Clone, levels = clones_list_all)



#save seurat
saveRDS(seurat_patient, file.path(rds_dir,paste0(patient,"three_samples_final_Clones.rds")))



#by Treatment
colors = c("cadetblue2", "tomato","plum")
plot <- DimPlot(seurat_patient, group.by= "treatment", cols = colors ) + labs(title = patient) + theme_void()
plot
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_DimPlot_all_cell.pdf")),width = 10, height= 8)

#Save data of Plot 4A
data = plot$data 
all(rownames(data) == rownames(seurat_patient@meta.data))

#Combine with Clone information
data = cbind(data, Clone = seurat_patient@meta.data$Clone)

write.table(data, file.path(Figures_data, paste0("Figure_4A_mapping_Clones_to_scRNA_",patient,".txt")), sep = "\t", quote = F)



#Include the non-genotyped cells
# Define colors for clones and NA cells
clone_colors <- c( `2` = "burlywood2",
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

# Extract UMAP embeddings and clone information
p <-  DimPlot(seurat_patient, reduction = "umap", group.by = "Clone") + theme_minimal()
plot_data <- p$data
# Convert `Clone` to factor and ensure `NA` is handled properly
plot_data$Clone <- factor(plot_data$Clone, levels =  clones_list_all)

# Plotting
ggplot(plot_data, aes(x = umap_1, y = umap_2, color = Clone)) +
  # First layer: NA cells in background (gray)
  geom_point(data = subset(plot_data, Clone == "UnI"), size = 0.5) +
  geom_point(data = subset(plot_data, Clone == "7"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "8"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "9"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "10"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "11"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "12"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "2"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "3"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "4"), size = 0.75) +
   geom_point(data = subset(plot_data, Clone == "5"), size = 0.75) +
  # Manual color scale for specific clone colors
  scale_color_manual(values = clone_colors)+
  #guides(color = guide_legend(reverse = TRUE))+
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
clones_list =   c( "2", "3","4","5", "7" ,"8","9","10","11","12")
seurat_patient_sub = subset(seurat_patient , subset = Clone %in% clones_list)
seurat_patient_sub@meta.data$Clone = factor(seurat_patient_sub@meta.data$Clone, levels = clones_list)

clone_colors <- c( `2` = "burlywood2",
            `3` = "mediumseagreen",
            `4` = "red",
            `5` = "tan4", 
            `7` = "plum", 
            `8` = "dodgerblue4", 
            `9` = "purple",
            `10` = "yellow2",
            `11` = "lightpink",
            `12` = "darkorange")






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


Clone_colors = list(Clone = clone_colors)


#Mutations likely to be detected
row_sum  = rowSums(mutation_matrix)
top_mut = head(order(row_sum,decreasing = T),30)
top_mut_id = Mutation_id_meta[top_mut, ] 
top_mut_info = top_mut_id  %>% left_join(annotation, by = "mutation_id") 

#gap for row and col
table(anno_mut$Clone)
rows_gap = c(  6 + 100 + 33 + 2,  6 + 100 + 33 + 2 + 1 + 14 +  4 )
table(anno_cell$Clone)
cols_gap = c( 23 + 1161 + 106+ 7 , 23 + 1161 + 106+ 7 +  1  +  18 +4)
#heatmap
pheatmap(mutation_matrix, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        gaps_row =  rows_gap ,gaps_col = cols_gap,
        #cellwidth = 0.15, cellheight = 1.5,  # Adjust cell size
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap.png"), width = 6,height=7,density = 600)

library(magick)
heatmap_img <- image_read(file.path(plot_dir, "mutation_matrix_3_samples_heatmap.png"), density = 600)
# Apply Gaussian Blur
blurred_img <- image_blur(heatmap_img, radius = 10, sigma = 1.5)

# Save the blurred image
image_write(blurred_img, path = file.path(plot_dir, "mutation_matrix_3_samples_heatmap_blurred.pdf"),  density = 600,format = "pdf")



#Save Figures data
write.table(mutation_matrix, file.path(Figures_data, paste0( "Figure_Sup_6C_mutation_matrix_", patient,".txt")), sep = "\t", quote= F)
write.table(anno_mut,  file.path(Figures_data, paste0( "Figure_Sup_6C_mutation_matrix_annotation_mutation_", patient,".txt")), sep = "\t", quote= F)
write.table(anno_cell,  file.path(Figures_data, paste0("Figure_Sup_6C_mutation_matrix_annotation_cell_", patient, ".txt")), sep = "\t", quote= F)






#Cell mutation count
mut_count =  colSums(mutation_matrix)
table(mut_count )

#PIK3CA
#mut_id = combined_filtered %>% filter(Gene_name == "PIK3CA")%>% pull(mutation_id) %>% unique()
#muts = c(":12:71132813:G:C",":3:143016296:C:T" ,":3:179234297:A:G",":12:22065183:G:A")

#rownames(mutation_matrix) = ifelse(rownames(mutation_matrix) %in% muts, rownames(mutation_matrix) ,""  )

#heatmap with rownames
# pheatmap(mutation_matrix, color = c("whitesmoke","darkred"),
        # cluster_rows = F, cluster_cols = F,
        # show_rownames = T, show_colnames = F,
        # annotation_row = anno_mut ,
        # annotation_col = anno_cell, annotation_colors = Clone_colors ,
        # legend_breaks = c(0, 1),
        # legend_labels = c("no","yes"),
        # filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap_rownames.png"), width = 6,height=7)
 
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
        gaps_row =  rows_gap,gaps_col = cols_gap,
        filename = file.path(plot_dir,"mutation_gene_expression_3_samples_heatmap.png"), width = 6,height=7)


#Save data
write.table(mut_scale_mat, file.path(Figures_data, paste0( "Figure_Sup_6C_expression_matrix_", patient,".txt")), sep = "\t", quote= F)



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
 
#Clone 3
combined_filtered_distinct = combined_filtered %>% distinct(CB_seurat,.keep_all = T)
Cell_meta_c = Cell_meta %>% filter(Clone == 3) %>% left_join(combined_filtered_distinct[,c("CB_seurat","sample")], by ="CB_seurat")
mutation_matrix_c = mutation_matrix[,Cell_meta_c$CB_seurat]

anno_cell_sample = data.frame(sample = Cell_meta_c$sample) %>% 
                    mutate(sample = case_when(sample == "OO100-no" ~ "Untreated",
                                              sample == "OO100-vitro" ~ 'Ex vivo',
                                              sample == "OO100-vivo" ~ 'In vivo'))
anno_cell_sample$sample = factor(anno_cell_sample$sample , levels = c("Untreated", "Ex vivo", "In vivo"))                                    
rownames(anno_cell_sample) = Cell_meta_c$CB_seurat

table(anno_cell_sample$sample)
anno_colors = list(
  # Column annotation (for 'sample' in anno_cell_sample)
  sample = c("Untreated" = "deepskyblue2", 
             "Ex vivo" = "tomato",  
             "In vivo" = "violet"),
  Clone = Clone_colors$Clone)

pheatmap(mutation_matrix_c, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell_sample, annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_3.png"), width = 6.5 ,height = 7,density = 600)


heatmap_img <- image_read(file.path(plot_dir, "mutation_clone_3.png"), density = 600)
# Apply Gaussian Blur
blurred_img <- image_blur(heatmap_img, radius = 10, sigma = 1.5)

# Save the blurred image
image_write(blurred_img, path = file.path(plot_dir, "mutation_clone_3_blurred.pdf"),  density = 600,format = "pdf")


muts = c(":12:66138011:G:C")


rownames(mutation_matrix_c) = ifelse(rownames(mutation_matrix_c) %in% muts, rownames(mutation_matrix_c) ,""  )
pheatmap(mutation_matrix_c, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = T, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell_sample, annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_3_rownames.png"), width = 10 ,height= 8)

#save data
write.table(mutation_matrix_c , file.path(Figures_data, paste0("Figure_4C_mutation_matrix_dominant_invivo_",patient,".txt")), sep = "\t", quote = F)



#HOW depleted THE CELL CLONE of mutations of clones from separated branches? ( how pure the clone defined is)
cell_anno = anno_cell
cell_anno = cell_anno %>% mutate(diff_branch_clones = case_when(Clone %in% c(2,3,4,5) ~ paste(7,8,9,10,11,12),
                                                                Clone %in% c(7,8,9) ~ paste(2,3,4,5,10,11,12),
                                                                 Clone %in% c(10,11,12) ~ paste(2,3,4,5,7,8,9),
                                                                ))
cell_anno$Clone = as.character(cell_anno$Clone)
rownames(mutation_matrix) <- Mutation_id_meta$mutation_id



#Function to calculate Cell clone purity
Purity_calculate <- function(Clone_current){
diff_branch_clones = cell_anno %>% filter(Clone == Clone_current) %>% pull(diff_branch_clones)%>% unique() %>% str_split(,pattern = " ") %>% unlist() 
CB_clone = cell_anno %>% filter(Clone  == Clone_current)%>% rownames()
#mutations in other  branches
diff_branch_muts = anno_mut %>% filter(Clone %in% diff_branch_clones ) %>% rownames()

mutation_matrix_clone_diff_branch = as.matrix(mutation_matrix[diff_branch_muts ,CB_clone])
Nb_contaminated_cells = sum(colSums(mutation_matrix_clone_diff_branch) > 0 )
Purity = 1 - Nb_contaminated_cells/length(CB_clone)
df = data.frame(Clone = Clone_current, Purity = Purity, n_cell = length(CB_clone) )
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
  annotate("text", x = 6, y = min_purity, label = paste("Min:", round(min_purity, 2)),
           vjust = 2.5, hjust = 1.1, color = "red", size = 11) +  # Add the text annotation
  theme_classic() +
  theme(legend.title =  element_text(size = 32, color = "black"),
        legend.text =  element_text(size = 32, color = "black"),
        axis.text.y.left = element_text(size = 32, color = "black"),
        axis.title.y.left = element_text(size = 32, color = "black"),
        axis.title.x.bottom = element_text(size = 32, color = "black"),
        axis.text.x.bottom = element_text(size = 32, color = "black"))

# Save the plot
ggsave(file.path(plot_dir, "Purity.pdf"), width = 11, height = 5, scale =1)

#save Figure data
write.table(Purity, file.path(Figures_data, "Figure_Sup_6C_Purity_cell_clone_OO100.txt"), sep = "\t", quote = F)





#Proportion of each clone
Clone_sample = as.data.frame(table(seurat_patient_sub$orig.ident, seurat_patient_sub$Clone))
colnames(Clone_sample) = c("sample","Clone","count")


Clone_sample <- Clone_sample %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100) # Convert counts to percentages

#Save the file
write.table( Clone_sample, file.path(txt_files, "Clone_proportion.txt"), sep = "\t", quote = F)


# Plot proportions as percentages
# Extract the unique order of clones based on the new order
ordered_clones <- Clone_sample %>%
  filter(sample == "OO100-no") %>%  # Filter for the "Untreated" sample
  arrange(desc(proportion)) %>%      # Arrange by proportion in descending order
  pull(Clone)

ordered_samples = c("OO100-no", "OO100-vitro", "OO100-vivo")
# Plot proportions as percentages
Clone_sample %>%
  ggplot(aes(x = factor(sample, levels = rev(levels(Clone_sample$sample))), y = proportion, fill = factor(Clone, levels = rev(ordered_clones)))) +
  geom_bar(stat = "identity", width = 0.8) + 
  scale_fill_manual(values = colors ) + 
  coord_flip() +
  ylab("Clone percentage (%)") +
  scale_x_discrete(labels= c("OO100-no" = "Untreated", "OO100-vitro" = "Ex vivo", "OO100-vivo" = "In vivo")) +
  theme_minimal() +
  theme(legend.position = "top",
     legend.direction = "horizontal",   
    legend.title = element_blank(), 
    legend.text = element_text(size = 20),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_text(size=20),
    axis.text.x.bottom = element_text(size=20),
    axis.title.x.bottom = element_text(size=20)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE) )

# Save the plot
ggsave(file.path(plot_dir, "OO100_clone_proportion.pdf"), width = 13000, height = 6000, dpi = 1500, units = 'px')

#vertical barplot
Clone_sample$Clone = factor(Clone_sample$Clone, levels = rev(ordered_clones))
plot <- Clone_sample %>%
  ggplot(aes(x = factor(sample, levels = ordered_samples), y = proportion, fill = Clone)) +
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = clone_colors ) + 
  ylab("Clone percentage (%)") +
  scale_x_discrete(labels= c("OO100-no" = "Untreated", "OO100-vitro" = "Ex vivo", "OO100-vivo" = "In vivo")) +
  theme_minimal() +
  theme(legend.position = "right",
    legend.title = element_text(size = 30, color = "black"), 
    legend.text = element_text(size = 30, color = "black"),
    axis.title.y.left = element_text(size = 30, color = "black" ),
    axis.text.y.left = element_text(size = 30,  color = "black"),
    axis.text.x.bottom = element_text(size = 30,  color = "black"),
    axis.title.x.bottom = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = FALSE) )

# Save the plot
plot
ggsave(file.path(plot_dir, "OO100_clone_proportion_vertical.pdf"), width = 8 , height = 8, scale = 1)

# Save Figure data
data <- plot$data %>% mutate( treatment = case_when(sample == paste0(patient,"-no") ~ "Untreated",
                                                    sample == paste0(patient,"-vitro") ~ "Ex vivo",
                                                    sample  == paste0(patient,"-vivo") ~ "In vivo"))

write.table(data, file.path(Figures_data,paste0( "Figure_4B_clone_proportion_", patient,".txt")), sep = "\t", quote = F)



## CONTAMINATE
 
#Contaminated cells of each clone (fraction of cells contain mutation from other clone)
Contaminate_cell_fraction = function(clone){
   CB_clone = anno_cell %>% filter(Clone == clone) %>% rownames()
   mut_mat_CB_clone = as.matrix(mutation_matrix[,CB_clone])
   clone_others = clones_list
   #calculate contaminate clone
   conta_clone_other = function(clone_other){
        mut_clone_other = anno_mut %>% filter(Clone == clone_other) %>% rownames()
        mut_mat_CB_clone_mut_other = as.matrix(mut_mat_CB_clone[mut_clone_other,, drop = FALSE])
        Nb_contaminated_cells = sum(colSums(mut_mat_CB_clone_mut_other) > 0 )
        fraction = 100 * Nb_contaminated_cells/length(CB_clone)
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

ggsave(file.path(plot_dir, paste0(patient,"_fraction_cell_in_mutation_clone.pdf")), width = 11, height = 12.5, scale=1)

#save Figure data
colnames(Contaminate_cell_each_clone ) = c("Cell_clone", "Mutation_clone", "Fraction")
write.table(Contaminate_cell_each_clone, file.path(Figures_data,paste0( "Figure_Sup_6C_Fraction_mutation_clone_in_cell_clone_",patient,".txt")), sep = "\t", quote = F)



## CONTAMINATE MUT
#do for mutation clone 

clones_list_al3 =  names(table(anno_mut$Clone))[table(anno_mut$Clone) >= 3]


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
