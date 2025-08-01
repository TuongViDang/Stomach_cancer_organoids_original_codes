library(tidyverse)
library(pheatmap)
library(RColorBrewer)



main_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA")
plot_dir = file.path(main_dir,"plots", "COPY_NUMBER")
dir.create(plot_dir, recursive = T)


Driver_cp_OO77 = read.delim("/group/poetsch_projects/poetsch_sc/CANOPY/OO77/Drivers_copy_number.txt")
Driver_cp_OO100 = read.delim("/group/poetsch_projects/poetsch_sc/CANOPY/OO100/Drivers_copy_number.txt")
Driver_cp_OO99 = read.delim("/group/poetsch_projects/poetsch_sc/CANOPY/OO99/Drivers_copy_number.txt")

mM_df = Driver_cp_OO77[,c(1,4:9)] %>% left_join(Driver_cp_OO99[,c(1,4:11)], by= "driver_gene") %>% 
                                             left_join(Driver_cp_OO100[,c(1,4:9)], by= "driver_gene") 

mM_matrix = as.matrix(mM_df[,c(2,5,3,6,4,7,8,12,9,13,10,14,11,15,16,19,17,20,18,21)])
rownames(mM_matrix) = mM_df[,1]
col_anno_patient = data.frame(patient = c(rep("OO77",6),rep("OO99",8),rep("OO100",6)), row.names = colnames(mM_matrix))
# Additional data frame for column annotation (sample)
col_anno_sample <- data.frame(
  sample = c( rep("U",2), rep("E",2),rep("I",2), rep("U",4), rep( "E",2),  rep("I",2), rep( "U",2), rep( "E",2), rep( "I",2) ),
  row.names = colnames(mM_matrix)
)

# Merge both annotations into a single data frame
col_anno <- cbind(col_anno_sample, col_anno_patient)


pheatmap::pheatmap(mM_matrix, cluster_col = F, cluster_row = T, gaps_col = c(2,4,6,8,10,12,14,16,18), show_colnames = F, annotation_col = col_anno,
color = c("white",colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(9)),breaks = seq(0,10,1),
         filename = file.path(plot_dir,"mM_heatmap.pdf"), width = 13, height=10)

