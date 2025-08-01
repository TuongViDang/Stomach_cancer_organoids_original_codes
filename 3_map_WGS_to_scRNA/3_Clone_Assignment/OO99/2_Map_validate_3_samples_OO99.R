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
patient="OO99"


main_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA",patient,"tree30_change_name")
plot_dir = file.path(main_dir,"plots","Map_validate_3_samples")
dir.create(plot_dir, recursive = T)
txt_files = file.path(main_dir,"txt_files")
dir.create(txt_files,recursive=T)
rds_dir = file.path(main_dir,"rds")
dir.create(rds_dir)

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


seurat_patient = subset(seurat_integrated, subset = patient == "OO99")

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
#Assigne Clone to seurat_patient
seurat_patient@meta.data$CB_seurat = colnames(seurat_patient)
meta = seurat_patient@meta.data
meta = meta %>% left_join(meta_all[,c("CB_seurat","Clone")], by = "CB_seurat")
seurat_patient@meta.data = meta[,c("orig.ident", "treatment", "CB_seurat", "Clone")]
rownames(seurat_patient@meta.data) = seurat_patient@meta.data$CB_seurat
table(seurat_patient@meta.data$Clone)


clones_list_all =   c( "2", "3","4","6", "7" ,"8","9","10","11","12", "UnI")
seurat_patient@meta.data = seurat_patient@meta.data %>% 
                               mutate( treatment = case_when(orig.ident == paste0(patient,"-no") ~ "Untreated",
                                                             orig.ident == paste0(patient,"-vitro") ~ "Ex vivo",
                                                             orig.ident  == paste0(patient,"-vivo") ~ "In vivo"))

seurat_patient@meta.data$treatment = factor(seurat_patient@meta.data$treatment , levels = c("Untreated", "Ex vivo", "In vivo"))
seurat_patient@meta.data$Clone = factor(seurat_patient@meta.data$Clone, levels = clones_list_all)

#by Treatment
colors = c("cadetblue2", "tomato","plum")
DimPlot(seurat_patient, group.by= "treatment", cols = colors) + labs(title = patient)
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_DimPlot_all_cell.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')

#Include the non-genotyped cells
# Define colors for clones and NA cells
clone_colors <- c( `2`= "aquamarine2", `3` = "moccasin" , `4` = "dodgerblue4", `6` ="yellow2",
            `7` = "lightseagreen", `8` = "lightsalmon4" , `9`= "plum" , 
             `10` = "tomato", `11` = "mediumseagreen", `12` = "purple", `UnI` = "gray92"  # Color for NA cells
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
  geom_point(data = subset(plot_data, Clone == "6"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "7"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "8"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "9"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "10"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "11"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "12"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "2"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "3"), size = 0.75) +
  geom_point(data = subset(plot_data, Clone == "4"), size = 0.75) +
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
clones_list =   c( "2", "3","4","6", "7" ,"8","9","10","11","12")
seurat_patient_sub = subset(seurat_patient , subset = Clone %in% clones_list)
seurat_patient_sub@meta.data$Clone = factor(seurat_patient_sub@meta.data$Clone, levels = clones_list)

#by Treatment
colors = c("cadetblue2", "tomato","plum")
DimPlot(seurat_patient_sub, group.by= "treatment", cols = colors) + labs(title = patient)
ggsave(file.path(plot_dir,paste0(patient,"_3_samples_sub_DimPlot.pdf")),width = 10000, height= 8000 ,dpi = 1500,units = 'px')
#by Clone
colors <- c( `2`= "aquamarine2", `3` = "moccasin" , `4` = "dodgerblue4", `6` ="yellow2",
            `7` = "lightseagreen", `8` = "lightsalmon4" , `9`= "plum" , 
             `10` = "tomato", `11` = "mediumseagreen", `12` = "purple")

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
#remove mutation in no cells


Clone_colors = list(Clone = colors)


#Mutations likely to be detected
row_sum  = rowSums(mutation_matrix)
top_mut = head(order(row_sum,decreasing = T),20)
top_mut_id = Mutation_id_meta[top_mut, ] 
top_mut_info = top_mut_id  %>% left_join(annotation, by = "mutation_id") 

#gap for row and col
table(anno_mut$Clone)
rows_gap = c( 2 + 33 +104 ,2 + 33 +104 +  8 +  6 + 16 +  3 +  3, 2 + 33+ 104 +  8 +  6 + 16 +  3 +  3 + 29 )
table(anno_cell$Clone)
cols_gap = c(  4 + 117 + 329 , 4 + 117 + 329 +56 +  8+  43 +  6 + 35, 4 + 117 + 329 +56 +  8+  43 +  6 + 35 + 109  )


pheatmap(mutation_matrix, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        gaps_row = rows_gap , gaps_col = cols_gap ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap.png"), width = 6,height=7)


#Cell mutation count
mut_count =  colSums(mutation_matrix)
table(mut_count )

#PIK3CA
#mut_id = combined_filtered %>% filter(Gene_name == "PIK3CA")%>% pull(mutation_id) %>% unique()
#muts = c(":12:71132813:G:C",":3:143016296:C:T" ,":3:179234297:A:G",":12:22065183:G:A")
muts = c( ":12:54011146:A:G" ,":2:45078495:T:C" , ":20:49710422:G:T" , ":16:57239788:C:A")
rownames(mutation_matrix) = ifelse(rownames(mutation_matrix) %in% muts, rownames(mutation_matrix) ,""  )

#heatmap with rownames
pheatmap(mutation_matrix, color = c("whitesmoke","darkred"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = T, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell, annotation_colors = Clone_colors ,
        legend_breaks = c(0, 1),
        legend_labels = c("no","yes"),
        filename = file.path(plot_dir,"mutation_matrix_3_samples_heatmap_rownames_mut_remove.png"), width = 6,height=7)
 
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
        filename = file.path(plot_dir,"mutation_gene_expression_3_samples_heatmap.png"), width = 6,height=7)

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
 
#Clone 4
combined_filtered_distinct = combined_filtered %>% distinct(CB_seurat,.keep_all = T)
Cell_meta_c = Cell_meta %>% filter(Clone == 4) %>% left_join(combined_filtered_distinct[,c("CB_seurat","sample")], by ="CB_seurat")
mutation_matrix_c = mutation_matrix[,Cell_meta_c$CB_seurat]

anno_cell_sample = data.frame(sample = Cell_meta_c$sample) %>% 
                    mutate(sample = case_when(sample == "OO99-no" ~ "Untreated",
                                              sample == "OO99-vitro" ~ 'Ex vivo',
                                              sample == "OO99-vivo" ~ 'In vivo'))
anno_cell_sample$sample = factor(anno_cell_sample$sample , levels = c("Untreated", "Ex vivo", "In vivo"))                                    
rownames(anno_cell_sample) = Cell_meta_c$CB_seurat

table(anno_cell_sample$sample)
anno_colors = list(
  # Column annotation (for 'sample' in anno_cell_sample)
  sample = c("Untreated" = "tomato", 
             "Ex vivo" = "deepskyblue2",  
             "In vivo" = "khaki3"),
  Clone = Clone_colors$Clone)

pheatmap(mutation_matrix_c, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = F, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell_sample, annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_4.png"), width = 8 ,height = 7)


muts = c(":11:63976355:C:T",":20:62308749:G:T",":12:54011146:A:G")


rownames(mutation_matrix_c) = ifelse(rownames(mutation_matrix_c) %in% muts, rownames(mutation_matrix_c) ,""  )
pheatmap(mutation_matrix_c, color = c("whitesmoke","red"),
        cluster_rows = F, cluster_cols = F,
        show_rownames = T, show_colnames = F,
        annotation_row = anno_mut ,
        annotation_col = anno_cell_sample, annotation_colors = anno_colors ,
        legend_breaks = c(0, 1),
        filename = file.path(plot_dir,"mutation_clone_4_rownames.png"), width = 10 ,height= 8)



#HOW depleted THE CELL CLONE of mutations of clones from separated branches? ( how pure the clone defined is)
cell_anno = anno_cell
cell_anno = cell_anno %>% mutate(diff_branch_clones = case_when(Clone %in% c(2,3,4) ~ paste(6,7,8,9,10,11,12),
                                                                Clone %in% c(6,7) ~ paste(2,3,4,11,12),
                                                                Clone == 8 ~ paste(2,3,4,9,10,11,12),
                                                                Clone == 9 ~ paste(2,3,4,8,10,11,12),
                                                                Clone == 10 ~ paste(2,3,4,8,9,11,12),
                                                                Clone == 11 ~ paste(2,3,4,6,7,8,9,10,12),
                                                                Clone == 12 ~ paste(2,3,4,6,7,8,9,10,11)
                                                                ))
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
df = data.frame(Clone = Clone_current, Purity = Purity, n_cell = length(CB_clone) )
df
}

Purity = purrr::map_dfr( clones_list,Purity_calculate )
colors = colors$Clone
Purity$Clone = factor(Purity$Clone , levels = names(colors))


# Calculate the minimum value of Purity
min_purity <- min(Purity$Purity, na.rm = TRUE)

# Create the bar plot with a horizontal line and display the minimum value
Purity %>% ggplot(aes(x = Clone, y = Purity, fill = Clone)) + 
  scale_fill_manual(values = colors) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n_cell), vjust = -0.5, size = 3.5, color= "blue" ) +  # Add n_cell values as labels above bars
  geom_hline(yintercept = min_purity, linetype = "dashed", color = "red", size = 1) +  # Add hline
  annotate("text", x = 5, y = min_purity, label = paste("Min:", round(min_purity, 2)),
           vjust = -0.5, hjust = 1.1, color = "red", size = 4) +  # Add the text annotation
  theme_classic()

# Save the plot
ggsave(file.path(plot_dir, "Purity.png"), width = 7, height = 4)




#Proportion of each clone
Clone_sample = as.data.frame(table(seurat_patient_sub$orig.ident, seurat_patient_sub$Clone))
colnames(Clone_sample) = c("sample","Clone","count")


Clone_sample <- Clone_sample %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100) # Convert counts to percentages

# Plot proportions as percentages
# Extract the unique order of clones based on the new order
ordered_clones <- Clone_sample %>%
  filter(sample == "OO99-no") %>%  # Filter for the "Untreated" sample
  arrange(desc(proportion)) %>%      # Arrange by proportion in descending order
  pull(Clone)


# Plot proportions as percentages
Clone_sample %>%
  ggplot(aes(x = factor(sample, levels = rev(levels(Clone_sample$sample))), y = proportion, fill = factor(Clone, levels = rev(ordered_clones)))) +
  geom_bar(stat = "identity", width = 0.8) + 
  scale_fill_manual(values = colors ) + 
  coord_flip() +
  ylab("Clone percentage (%)") +
  scale_x_discrete(labels= c("OO99-no" = "Untreated", "OO99-vitro" = "Ex vivo", "OO99-vivo" = "In vivo")) +
  theme_minimal() +
  theme(legend.position = "top",
     legend.direction = "horizontal",   
    legend.title = element_blank(), 
    legend.text = element_text(size = 16),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_text(size=18),
    axis.text.x.bottom = element_text(size=17),
    axis.title.x.bottom = element_text(size=17)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE) )

# Save the plot
ggsave(file.path(plot_dir, "OO99_clone_proportion.pdf"), width = 13000, height = 6000, dpi = 1500, units = 'px')


## CONTAMINATE

#Contaminated cells of each clone (fraction of cells contain mutation from other clone)
Contaminate_cell_fraction = function(clone){
   CB_clone = anno_cell %>% filter(Clone == clone) %>% rownames()
   mut_mat_CB_clone = mutation_matrix[,CB_clone]
   clone_others = setdiff(clones_list,clone)
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
ggplot(Contaminate_cell_each_clone , aes(x = as.factor(Clone), y = Fraction, fill = as.factor(Clone_other))) +
  geom_bar(stat = "identity", position = "dodge") +   # Stacked bar for each clone
  labs(x = "Cell Clone", y = "Fraction of cell", fill = "Mutation Clone") +
  theme_minimal() +
  scale_fill_manual(values= colors)    +      
  ggtitle("Contamination Fraction in cell clone")

ggsave(file.path(plot_dir, paste0(patient,"_clone_contaminated.pdf")), width = 14000, height = 5000, dpi = 1500, units = 'px')

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

pvalue_results = purrr::map_dfr(clones_list_al3,~ Bootstrap_contaminate_mut_fraction(.x, clones_cell_list = clones_list_al3, mutation_matrix = mutation_matrix,anno_mut= anno_mut, anno_cell = anno_cell, plot_dir_boostrap= plot_dir_boostrap))
pvalue_results$patient = "OO99"
#save pvalue_results
write.table(pvalue_results, file.path(txt_files,paste0(patient,"_bootstrap_pvalues.txt")), sep = "\t", quote = F)


#Load all pvalues of 3 patients
OO77 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO77/tree30_change_name/txt_files/OO77_bootstrap_pvalues_merged_separated_c5.txt")
OO99 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO99/tree30_change_name/txt_files/OO99_bootstrap_pvalues.txt")
OO100 = read.delim("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO100/tree30_change_name/txt_files/OO100_bootstrap_pvalues_remove_HLA.txt")


All_pvalue = rbind(OO77, OO99, OO100)
All_pvalue$padjust =  p.adjust(All_pvalue $p_value, method = "BH")

pvalue_results = All_pvalue  %>% filter(patient == "OO99")
pvalue_result_wider = pvalue_results %>% pivot_wider(id_cols = c("clone_mut"),names_from = "clone_cell", values_from = "padjust") %>% as.data.frame()

rownames(pvalue_result_wider) = pvalue_result_wider$clone_mut 
pvalue_result_wider = pvalue_result_wider[,2:ncol(pvalue_result_wider)]

pheatmap(pvalue_result_wider , colorRampPalette(brewer.pal(n = 11, name ="RdBu")[c(2,4,5,8)])(4), 
         breaks = c(0,0.005,0.01,0.05,0.1), legend_breaks = c(0,0.005,0.01,0.05,0.1),
         show_colnames = T, show_rownames = T,fontsize = 8,
         cluster_cols = F,cluster_rows = F,
         gaps_col = c(2,7), gaps_row = c(2,7),
         filename = file.path(plot_dir_boostrap,"pvalue_all.pdf"),
          width = 4, height = 4)

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
  mutations = Combined_all %>% filter( treeCLUSTER == clone) %>% pull(mutation_id) %>% unique()

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

clones = setdiff(unique(Combined_all$treeCLUSTER), c("7","10"))

Threshold_test_intra_clone <- function(threshold){
df = purrr::map_dfr(clones,~Intra_Clone_Dist_score(.x, threshold = threshold ))
df$Clone = factor(as.character(df$Clone), levels = clones)
df$threshold = threshold
df
}

thresholds = 6
All_thresholds_intra_clone = purrr::map_dfr(thresholds,Threshold_test_intra_clone)


for (clone in clones){
Clone  = All_thresholds_intra_clone %>% filter(Clone == clone) %>% pull(dist)
nmuts = All_thresholds_intra_clone %>% filter(Clone == clone) %>% pull(n_muts) %>% unique()
png(file.path(plot_Intra_Clone_score_dir,paste0("Clone_", clone,".png")))
hist(Clone , breaks =25, main = paste0("Clone ",clone,"-",nmuts), cex.main = 3, xlab = "pairwise intra clone distance" , cex.lab = 2)
dev.off()
}


## Criteria: Mutation need to co-occur with at least 1 other mutation in the same branch | not co-occur with mutation from other branch
tree = list( c(2,3,4),
             c(6,7,8,9,10),
             c(11),
             c(12))


## each mut 
Filter_criteria_each_mut <- function(i){
mut_i = rownames(mutation_matrix)[i]
mut_clone = anno_mut[mut_i,"Clone"]

#Find same branch and different branch
for (c in seq_along(tree)){
   if (sum(mut_clone == tree[[c]]) >0 ) {
    same_branches = tree[[c]]
    diff_branches = tree[-c] %>% unlist()
    }
    }
#Find mutations that co-occur with this mut
result = TRUE
#same branch
same_branch_mut = anno_mut %>% filter(Clone %in% same_branches ) %>% rownames() %>% setdiff(mut_i)
same_branch_cooccur_sum = 0
for (mut_j in same_branch_mut ){
j = which(rownames(mutation_matrix)== mut_j)
mat_mut_ij = mutation_matrix[c(i,j),]
co_occur_same = ifelse(sum(colSums(mat_mut_ij) >1) > 0, T, F)
same_branch_cooccur_sum = same_branch_cooccur_sum + co_occur_same
   }
#different branch
diff_branch_mut = anno_mut %>% filter(Clone %in% diff_branches ) %>% rownames() 
diff_branch_cooccur_sum = 0
for (mut_k in diff_branch_mut ){
k = which(rownames(mutation_matrix) == mut_k)
mat_mut_ik = mutation_matrix[c(i,k),]
co_occur_diff = ifelse(sum(colSums(mat_mut_ik) >1) > 0, T, F)
diff_branch_cooccur_sum = diff_branch_cooccur_sum + co_occur_diff
   }
if ((same_branch_cooccur_sum == 0 ) & (diff_branch_cooccur_sum >= 2) ) {result = FALSE}
result
}

Filter_results <- purrr::map(1:nrow(mutation_matrix), Filter_criteria_each_mut) %>% unlist()
table(Filter_results )

mut_remove = rownames(anno_mut)[which(Filter_results == F)]
mut_remove_df = annotation %>% filter(mutation_id %in% mut_remove)

write.table(mut_remove_df, file.path(txt_files, paste0(patient,"mut_co_occur_remove.txt")), sep = "\t",quote = F)
