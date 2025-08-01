
# Single-cell RNA-seq analysis - integration
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(patchwork)
library(RColorBrewer)


#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path<-file.path(main_path,'rds',"three_patients")
plot_3_path<-file.path(main_path,'plots_3_patients')
plot_path_results <- file.path(plot_3_path,'after_qc')
dir.create(plot_path_results)
#Figures data path
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)




functions_path <- '/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/downstream_analysis'

source(file.path(functions_path,'functions.R'))


#Load input file

seurat_integrated<-readRDS(file.path(rds_3_path,  paste0("three_patient_integrated_cca_ref.rds")))


#change sample names
meta <- seurat_integrated@meta.data
meta <- meta %>%
  mutate(Sample = case_when(
    grepl("-no", orig.ident) ~ sub("-no", " Untreated", orig.ident),
    grepl("-vitro", orig.ident) ~ sub("-vitro", " Ex vivo", orig.ident),
    grepl("-vivo", orig.ident) ~ sub("-vivo", " In vivo", orig.ident),
    TRUE ~ orig.ident
  ))
seurat_integrated@meta.data <- meta

#set color
color_palette <- c(brewer.pal(12, "Set3")[1:8],  brewer.pal(8, "Set2")[4])
#Plot number of cells
meta <- seurat_integrated@meta.data 
sample_order = c("OO77 Untreated"  , "OO77 Ex vivo" ,  "OO77 In vivo" ,
                    "OO99 Untreated"  ,    "OO99 Ex vivo" , "OO99 In vivo" ,
                    "OO100 Untreated" , "OO100 Ex vivo" , "OO100 In vivo" )

meta$Sample <- factor(meta$Sample, levels = sample_order)
meta_sum <- meta %>% group_by(Sample) %>% summarize(n_cells= n())
p1 = meta_sum %>% ggplot(aes(x = Sample , y= n_cells, fill = Sample )) +
             geom_bar(stat = "identity", width = 0.6)+
             geom_text(aes(label = n_cells), vjust = -0.5, size = 10) +
             ylim(0,9000)+
             scale_fill_manual(values = color_palette  )+
             ylab("Number of cell") +
             xlab("Sample")+ 
             theme_classic() +
             theme(legend.position= "None", 
               legend.title = element_text(size = 28, color = "black"),
                legend.text = element_text(size = 28, color = "black"),
                axis.text.x = element_blank(),
                axis.title.x.bottom = element_blank(),
                plot.title = element_text(size = 28, color = "black"),
                axis.text.y.left = element_text(size = 28, color = "black"),
                axis.title.y.left = element_text(size = 28, color = "black"))
             # + ggtitle("Number of cell")
p1             
ggsave(file.path(plot_path_results,"nCells.pdf"), width = 14, height = 5)          

data = p1$data
write.table(data, file.path(Figures_data, "Figure_Sup_1A_number_of_cell_per_sample.txt"),sep = "\t", quote = F) 

#Boxplot
options(scipen = 999)  
p2 = meta %>% ggplot(aes(x = Sample, y =  nCount_RNA,group = Sample,  fill =  Sample))+
         geom_boxplot()+
         scale_fill_manual(values =color_palette )+
         ylim(0,130000)+
          xlab("Sample")+ 
          ylab("Number of UMI / cell") +
          theme_classic() +
             theme(legend.position= "None", 
                axis.text.x = element_blank(),
                axis.title.x.bottom = element_blank(),
                plot.title = element_text(size = 28),
                   axis.text.y.left = element_text(size = 28, color = "black"),
                   axis.title.y.left = element_text(size = 28, color = "black")) 
           # + ggtitle("Number of UMI/cell")
ggsave(file.path(plot_path_results,"nUMI_boxplot.pdf"), width = 12, height = 5)

data = p2$data
write.table(data, file.path(Figures_data, "Figure_Sup_1A_nUMI_per_cell_per_sample.txt"),sep = "\t", quote = F) 


p3 =  meta %>% ggplot(aes(x = Sample, y =  nFeature_RNA,group = Sample,  fill =  Sample))+
         geom_boxplot()+
         scale_fill_manual(values = color_palette  )+
         ylim(0,15000)+
          xlab("Sample")+ 
          ylab("Number of gene/cell") +
          theme(legend.position= "None")+
          theme_classic() +
             theme(legend.position= "None", 
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 28, color = "black"),
                axis.title.x.bottom = element_blank(),
                plot.title = element_text(size = 38, color = "black"),
                axis.text.y.left = element_text(size = 28, color = "black"),
                axis.title.y.left = element_text(size = 28, color = "black")) 
         # + ggtitle("Number of gene/ cell")

p1 + p2 + p3 + plot_layout(nrow = 3)
ggsave(file.path(plot_path_results,"combined_cell_UMI_gene.pdf"), width = 14, height = 15)

data = p3$data
write.table(data, file.path(Figures_data, "Figure_Sup_1A_nGene_per_cell_per_sample.txt"),sep = "\t", quote = F) 


#Calculate median
#median number of cell per sample
median(meta_sum$n_cells)
#median number of UMNI per cell
median(meta$nCount_RNA)
#median number of gene per cell
median(meta$nFeature_RNA)