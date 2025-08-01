
#singularity exec /group/biotec_poetsch/alex/r431_mantic.sif R
library(cloneMap)
library(tidyverse)
library(patchwork)
  

#Conipher path
CONIPHER_DIR="/group/poetsch_projects/poetsch_sc/Conifer/tree30_change_name"
txt_dir=file.path(CONIPHER_DIR,"txt_files")
dir.create(txt_dir, recursive = T)
plot_dir=file.path(CONIPHER_DIR,"FishPlot_CloneMap")
dir.create(plot_dir, recursive = T)

#Path to save Figures data
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)


patient = "OO77"
tree_df = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/tree_30_final/Trees/allTrees.txt"),comment.char = "#",header=F)
colnames(tree_df) = c("parent","clusterID")
#Change name of clone
# Define mapping for replacements
clusterID_map <- c(`10` = 1, `1` = 3, `6` = 4, `2` = 5, `5` = 6, `11` = 7, `9` = 8, `3` = 2)


# Apply transformations
tree_df <- tree_df %>%
  mutate(
    clusterID = ifelse(clusterID %in% names(clusterID_map), clusterID_map[as.character(clusterID)], clusterID),
    parent = ifelse(parent %in% names(clusterID_map), clusterID_map[as.character(parent)], parent)
  )

Find_direct_children <- function(cluster,tree){
    direct_children = tree %>% filter(parent == cluster) %>% pull(clusterID)
    direct_children
}

Find_all_children <- function(cluster,tree) {
  all_children <- vector()
  direct_children <- Find_direct_children(cluster,tree=tree)
  all_children <- c(all_children, direct_children)
  for (child in direct_children) {
    grandchildren <- Find_all_children(child, tree = tree)
    all_children <- c(all_children, grandchildren)
  }
  return(all_children)
}

Clones_all = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/tree_30_final/Trees/cloneProportionsMinErrorTrees.txt"))
#Change name of clone
Clones_all <- Clones_all%>%
  mutate(
    clusterID = ifelse(clusterID %in% names(clusterID_map), clusterID_map[as.character(clusterID)], clusterID)
  )
Clones_all <- Clones_all  %>% left_join(tree_df,by = "clusterID")

New_proportion <- function(cluster,sample,clones_df,tree){
    children = Find_all_children(cluster,tree = tree)
    children = c(cluster,children)
    children_in_sample = clones_df[(clones_df[sample] > 0 ) & clones_df$clusterID %in% children,]
    Sum_Pro_in_Sample = sum(children_in_sample[sample])
    return(Sum_Pro_in_Sample)
    
}

tumor_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_tumor"),clones_df = Clones_all, tree=tree_df)) %>% unlist()
ex_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_ex"),clones_df = Clones_all, tree=tree_df)) %>% unlist()
in_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_in"),clones_df = Clones_all, tree=tree_df)) %>% unlist()

Clones_all_new = data.frame(clusterID = Clones_all$clusterID, parent= Clones_all$parent,tumor = tumor_new, ex_vivo = ex_new, in_vivo = in_new)



write.table(Clones_all_new,file.path(txt_dir,paste0("Clone_proportion_include_children_",patient,".txt")), sep = "\t",quote = F)

#and providing both timepoints to label and a plot title
cols <- c("darkkhaki", "dodgerblue4", "lightgreen", "mediumseagreen", "lightpink",  "lightgoldenrod1","cadetblue2","plum") 
names(cols) <- Clones_all_new$clusterID

#Untreated 
CCF_Untreated <- Clones_all_new[,c(1,3)]
colnames(CCF_Untreated) <- c("clones","CCF")
CCF_Untreated$CCF <- round(CCF_Untreated$CCF/100,2)
CCF_Untreated <- CCF_Untreated %>% filter(CCF > 0)




pdf(file.path(plot_dir,paste0("CloneMap_",patient,"_Untreated.pdf")),width= 8, height = 6)
cloneMap(tree_df, CCF_Untreated, clone.cols=cols)
dev.off()

#Ex vivo 
CCF_Ex <- Clones_all_new[,c(1,4)]
colnames(CCF_Ex) <- c("clones","CCF")
CCF_Ex$CCF <- round(CCF_Ex$CCF/100,2)
CCF_Ex <- CCF_Ex %>% filter(CCF > 0)

pdf(file.path(plot_dir, paste0("CloneMap_",patient,"_EX.pdf")),width= 8, height = 6)
cloneMap(tree_df, CCF_Ex, clone.cols=cols)
dev.off()

#In vivo 
CCF_In <- Clones_all_new[,c(1,5)]
colnames(CCF_In) <- c("clones","CCF")
CCF_In$CCF <- round(CCF_In$CCF/100,2)
CCF_In <- CCF_In %>% filter(CCF > 0)
CCF_In[1,2] <- 0.999
df = data.frame( clones = 2, CCF = 0.001)
CCF_In <-  rbind(CCF_In , df)
pdf(file.path(plot_dir,paste0("CloneMap_",patient,"_IN.pdf")),width= 8, height = 6)
cloneMap(tree_df, CCF_In, clone.cols=cols)
dev.off()


###### Recent expanding clone  ######


Observed_CCF_clone_proportion_plots <- function(Clones_all_new, Clones_all, value_column, sample_name, plot_dir, patient, colors ) {
 
  # Filter and arrange for the specific condition
  Clones_all_new_sample <- Clones_all_new %>% 
    filter(!!sym(value_column) > 0) %>% 
    arrange(desc(!!sym(value_column)))
  
  Clones_all_new_sample$clusterID <- as.character(Clones_all_new_sample$clusterID)
   names(cols) <- Clones_all_new$clusterID
  colors_filtered <- colors[Clones_all_new_sample$clusterID]
  Clones_all_new_sample$clusterID <- factor(Clones_all_new_sample$clusterID, levels = sort(Clones_all_new_sample$clusterID))
  
  # Plot for the specific condition
  p1 <- Clones_all_new_sample %>% 
    ggplot(aes_string(x = "clusterID", y = value_column, fill = "clusterID")) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_manual(values = colors_filtered) + 
    xlab("Clone")+
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_blank(),
          axis.title.y.left = element_blank(), 
          axis.text.y.left = element_text(size = 32, color = "black"),
          axis.text.x.bottom = element_text(size = 32 , color = "black"),
          axis.title.x.bottom = element_text(size = 32))
  
  #ggsave(file.path(plot_dir, paste0(patient, "_Observed_CCF_", value_column, ".pdf")), plot = p1)
  
  # Pivot longer and filter for clone proportion
  Clones_all_longer <- Clones_all %>% 
    pivot_longer(cols = 1:3, names_to = "Sample", values_to = "Proportion")
  
  Clones_all_longer_sample <- Clones_all_longer %>% 
    filter((Sample == sample_name) & (Proportion > 0)) %>% 
    arrange(desc(Proportion))
  
  Clones_all_longer_sample$clusterID <- as.character(Clones_all_longer_sample$clusterID)
  colors_prop <- colors[Clones_all_longer_sample$clusterID]
  Clones_all_longer_sample$clusterID <- factor(Clones_all_longer_sample$clusterID, levels = rev(Clones_all_longer_sample$clusterID))
  
  # Plot for clone proportion
  label_data <- Clones_all_longer_sample %>%
  group_by(Sample) %>%
  mutate(
    percent = Proportion / sum(Proportion),
    label = sprintf("%.0f", percent * 100),  
    pos = cumsum(Proportion) - 0.5 * Proportion
  )

  p2 <- ggplot(label_data, aes(x = Sample, y = Proportion, fill = clusterID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +
  geom_text(aes(y = pos, label = label), color = "black",hjust = 2, size = 11) +
  scale_fill_manual(values = colors_prop) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_blank()
  )

  #ggsave(file.path(plot_dir, paste0(patient, "_Clone_proportion_", value_column, ".pdf")), plot = p2)
  
  # Combine plots
  combined_plot <- p1 + p2 + plot_layout(widths = c(4, 3)) +
    plot_annotation(title = "Clone proportion (%)",theme = theme(plot.title = element_text(size = 32)))
  
  ggsave(file.path(plot_dir, paste0(patient, "_Observed_CCF_Clone_proportion_", value_column, ".pdf")), plot = combined_plot, width = 6.5, height = 4, scale = 1)

Clones_all_longer_sample 
}

#Untreated 
Untreated = Observed_CCF_clone_proportion_plots (
  Clones_all_new = Clones_all_new,
  Clones_all = Clones_all, 
  value_column = "tumor", 
  sample_name = "OO77_tumor", 
  plot_dir = plot_dir, 
  patient = patient, 
  colors = cols
)
#Ex vivo
Ex_vivo = Observed_CCF_clone_proportion_plots(
  Clones_all_new = Clones_all_new,
  Clones_all = Clones_all, 
  value_column = "ex_vivo", 
  sample_name = "OO77_ex", 
  plot_dir = plot_dir, 
  patient = patient, 
  colors = cols
)
#In vivo
In_vivo = Observed_CCF_clone_proportion_plots(
  Clones_all_new = Clones_all_new,
  Clones_all = Clones_all, 
  value_column = "in_vivo", 
  sample_name = "OO77_in", 
  plot_dir = plot_dir, 
  patient = patient, 
  colors = cols
)

#Save data
write.table(Clones_all_new,file.path(Figures_data, paste0("Figure_3_bulkWGS_Clone_evolution_",patient,".txt")), sep = "\t", quote = F)



#Reccent expanding clone (not use)
Recent_expanding_clone <- function(clones_df) {
   terminal_nodes = clones_df$clusterID[ !(clones_df$clusterID %in% tree_df$parent)]
   terminal_df <- clones_df %>% filter(clusterID %in% terminal_nodes)
   last_expand_clone_df <- terminal_df[which.max(terminal_df$Proportion),]     
   last_expand_clone_df    }

df_list = list(Untreated, Ex_vivo, In_vivo)
recent_expanding_clones <- purrr::map_dfr(df_list, Recent_expanding_clone)


recent_expanding_clones$Sample <- factor(recent_expanding_clones$Sample, levels = c("OO77_tumor", "OO77_in", "OO77_ex"))
color_clones <- cols[c("6","26","32")]
names(color_clones) <-  c("OO77_tumor", "OO77_in", "OO77_ex")

recent_expanding_clones %>% ggplot(aes(x = Sample, y = Proportion, fill = Sample)) + geom_bar(stat = "identity", width = 0.6) +
 scale_fill_manual(values = color_clones) +
  scale_y_continuous(limits = c(0,100))+
 theme_minimal()+ 
 theme(legend.position = "none", axis.title.x.bottom = element_blank(),axis.text.x.bottom = element_text(size = 15)) +
  scale_x_discrete(labels = c("OO77_tumor" = "Untreated","OO77_ex" = "Ex vivo", "OO77_in" = "In vivo" )) +
  ylab("Tumor proportion of the largest terminal clone (%)")


ggsave(file.path(plot_dir, paste0(patient, "_Recent_expanding_clone.pdf")), width = 4, height = 4, scale = 1)


                            