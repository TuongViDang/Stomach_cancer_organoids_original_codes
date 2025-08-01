## Load libraries
library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(pheatmap)
library(CellChat)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(magrittr)
library(ggpubr)

#Specify the path
main_path <-'/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_path,"Liana")
dir.create(plot_result_path)

#Figures data path
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)



source("/group/poetsch_users/poetsch_vi/liana/R/liana_plot.R")

seurat_integrated <- readRDS(file.path(rds_3_path, "clustered_0.35_integrated_seurat_ref.rds"))
cluster_sizes = as.data.frame(table(seurat_integrated@meta.data$orig.ident,seurat_integrated@meta.data$seurat_clusters ))
colnames(cluster_sizes) = c("orig.ident", "cluster","count")



context_df_dict <- readRDS(file.path(rds_3_path ,"Liana_context_df_dict.rds"))

untreateds = c(1,4,7)

dict_sub_df <- bind_rows(context_df_dict[untreateds])
# Group by unique combinations of source, target, ligand.complex, and receptor.complex
# Then sum aggregate_rank and mean_rank
summarized_data <- dict_sub_df %>%
  group_by(source, target, ligand.complex, receptor.complex) %>%
  summarize(
    total_aggregate_rank = sum(aggregate_rank, na.rm = TRUE),
    total_mean_rank = sum(mean_rank, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(total_aggregate_rank)



untreated_df <- do.call(rbind, context_df_dict[c(1, 4, 7)])

pdf(file.path(plot_result_path,"ligand_receptors_c9_untreated.pdf"), width = 10, height = 9)
untreated_df  %>%
  liana_dotplot(source_groups = c("9"),
                ntop = 10)
dev.off()



vitro_df <- do.call(rbind, context_df_dict[c(2, 5, 8)])
pdf(file.path(plot_result_path,"ligand_receptors_c9_vitro.pdf"), width = 10, height = 9)
vitro_df  %>%
  liana_dotplot(source_groups = c("9"),
                ntop = 10)
dev.off()



vivo_df <- do.call(rbind, context_df_dict[c(3, 6, 7)])
pdf(file.path(plot_result_path,"ligand_receptors_c9_vivo.pdf"), width = 10, height = 9)
vivo_df  %>%
  liana_dotplot(source_groups = c("9"),
                ntop = 10)
dev.off()



#Automatic subset
sets = list(Untreated = c(1,4,7), vitro = c(2,5,8),vivo = c(3,6,9))

for (index in 1:3){
sample = names(sets)[index]
dict_sub <- context_df_dict[sets[[index]]]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))
top_100$source = factor(top_100$source, levels = 1:9)
top_100$target = factor(top_100$target, levels = 1:9)
df = as.matrix(table(top_100$source,top_100$target))

pheatmap(df,show_rownames=T,show_colnames=T,cluster_col=F,display_numbers=F,cluster_row=F,fontsize=17,
color = colorRampPalette(brewer.pal(n = 7, name ="OrRd")[c(1,3,6,7)])(80),breaks = seq(0,40,0.5),
filename=file.path(plot_result_path,paste("top_100_",sample,".pdf")))

#save Figure 
write.table(df, file.path(Figures_data, paste0("Figure_6E_heatmap_sender_receiver_",sample, ".txt")), sep = "\t", quote = F)

}



#Top Ligand-Receptor
for (index in 1:3){
sample = names(sets)[index]
dict_sub <- context_df_dict[sets[[index]]]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))


#ligand
top_100_ligand <- top_100 %>% group_by(ligand.complex) %>% summarize(count = n(),cluster=source,receptor.complex=receptor.complex) %>% arrange(desc(count)) 
top_100_ligand$cluster <- factor(top_100_ligand$cluster,levels= 1:9)
top_lig = unique(top_100_ligand$ligand.complex ) %>% head(10)
top_100_ligand <- top_100_ligand %>% filter(ligand.complex %in% top_lig)
top_100_ligand$ligand.complex <- factor(top_100_ligand$ligand.complex,levels= rev(top_lig))

#
top_100_receptor <- top_100 %>% group_by(receptor.complex) %>% summarize(count = n(),cluster=source) %>% arrange(desc(count)) 
top_100_receptor$cluster <- factor(top_100_receptor$cluster,levels=0:13)
top_rec = unique(top_100_receptor$receptor.complex ) %>% head(10)
top_100_receptor  <- top_100_receptor %>% filter(receptor.complex %in% top_rec)
top_100_receptor$receptor.complex <- factor(top_100_receptor$receptor.complex,levels= rev(top_rec))

#Set color
custom_palette <- c(  "#80B1D3", "limegreen", "#FCCDE5" ,
                  "#FDB462","#D9D9D9", "#CCEBC5" ,
                 "#BC80BD" , "khaki1",  "brown1")   

top_100_ligand$count = rep(1, nrow(top_100_ligand))
top_100_receptor$count = rep(1, nrow(top_100_receptor))
# Create plot for top 200 ligands
p1 <- top_100_ligand  %>%
  ggplot(aes(x = ligand.complex, y = count, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", width=0.8) +
  scale_fill_manual(values = custom_palette) + theme_bw()+
  coord_flip() +
  labs(title = "Top Ligands") + 
  theme(legend.text = element_text(size = 28, color = "black"),
       axis.title.y.left = element_blank(),
       axis.text.y.left = element_text(size = 28, color = "black"),
       axis.text.x.bottom = element_text(size = 28,color = "black" ),
        axis.title.x.bottom = element_text(size = 28, color = "black"))
data <- p1$data 
write.table(data, file.path(Figures_data, paste0("Figure_6G_top_ligands_",sample, ".txt")), sep = "\t", quote = F)

# Create plot for top 100 receptors
p2 <- top_100_receptor   %>%
  ggplot(aes(x = receptor.complex, y = count, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack",width=0.8) +
  scale_fill_manual(values = custom_palette) +  theme_bw()+
  coord_flip() +
  labs(title = "Top Receptors") +
  theme( legend.text = element_text(size = 28, color = "black"),
       axis.text.y.left = element_text(size = 28, color = "black"),
       axis.text.x.bottom = element_text(size = 28, color = "black"),
       axis.title.y.left = element_blank(),
       axis.title.x.bottom = element_text(size = 28, color = "black"))

# Combine plots and share legend

ggpubr::ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="right")
ggsave(file.path(plot_result_path,paste0(sample,"_top_100_20_Ligand_Receptor_barplot.pdf")),width=15,height=8)



#Top Ligand and their interaction partners
set2_palette <- brewer.pal(8, "Set2")
set3_palette <- brewer.pal(12, "Set3")
combined_palette <- c(set2_palette, set3_palette, brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
p3 <- top_100_ligand  %>%
  ggplot(aes(x = ligand.complex, y = count, fill = receptor.complex)) +
  geom_bar(stat = "identity", position = "stack",width=0.8) + 
  scale_fill_manual(values = combined_palette) + theme_bw()+
  coord_flip() +
  labs(title = "Top Ligands partners") + 
   theme( legend.text =  element_text(size = 28, color = "black"),
        axis.text.y.left = element_text(size = 28, color = "black"),
       axis.text.x.bottom = element_text(size = 28, color = "black"),
       axis.title.y.left = element_blank(),
       axis.title.x.bottom = element_text(size = 28, color = "black"))


data <- p3$data 
write.table(data, file.path(Figures_data, paste0("Figure_Sup_9A_top_ligands_Rpartners_",sample, ".txt")), sep = "\t", quote = F)


ggpubr::ggarrange(p1, p3, ncol=2, legend="right")
ggsave(file.path(plot_result_path,paste0(sample,"top_100_20_Ligand_Cluster_Partner_barplot.pdf")),width=15,height=8)
}



####### Plot receptor with uniformed colors ######
# Initialize variables to store unique receptors
all_receptors <- c()

# First loop to gather all unique receptors across both samples
for (index in c(1, 3)) {
sample = names(sets)[index]
dict_sub <- context_df_dict[sets[[index]]]
top_100 <- purrr::map_dfr(dict_sub, ~ head(.x,100))
top_100_ligand <- top_100 %>% group_by(ligand.complex) %>% summarize(count = n(),cluster=source,receptor.complex=receptor.complex) %>% arrange(desc(count)) 
top_100_ligand$cluster <- factor(top_100_ligand$cluster,levels= 1:9)
top_lig = unique(top_100_ligand$ligand.complex ) %>% head(10)
top_100_ligand <- top_100_ligand %>% filter(ligand.complex %in% top_lig)
top_100_ligand$ligand.complex <- factor(top_100_ligand$ligand.complex,levels= rev(top_lig))
top_rec <- unique(top_100_ligand$receptor.complex) 
all_receptors <- union(all_receptors, top_rec)
}

# Create a consistent color palette for all receptors
set2_palette <- brewer.pal(8, "Set2")
set3_palette <- brewer.pal(12, "Set3")
set1_palette <- brewer.pal(8, "Set1")
dark2_palette <- brewer.pal(8, "Dark2")
combined_palette <- c(set2_palette, set3_palette, set1_palette,dark2_palette )

# Ensure that we have enough colors for all receptors
receptor_colors <- setNames(combined_palette[1:length(all_receptors)], all_receptors)


# Initialize a list to store plots
plots <- list()

# Second loop to generate plots with consistent colors
for (index in c(1, 3)) {
  sample = names(sets)[index]
  dict_sub <- context_df_dict[sets[[index]]]
  top_100 <- purrr::map_dfr(dict_sub, ~ head(.x, 100))

  # Process ligands
  top_100_ligand <- top_100 %>% 
    group_by(ligand.complex) %>% 
    summarize(count = n(), cluster = source, receptor.complex = receptor.complex) %>% 
    arrange(desc(count))
  
  top_100_ligand$cluster <- factor(top_100_ligand$cluster, levels = 0:13)
  top_lig <- unique(top_100_ligand$ligand.complex) %>% head(10)
  top_100_ligand <- top_100_ligand %>% 
    filter(ligand.complex %in% top_lig)
  top_100_ligand$ligand.complex <- factor(top_100_ligand$ligand.complex, levels = rev(top_lig))

  # Apply consistent colors based on the common receptor set
  top_100_ligand$receptor.complex <- factor(top_100_ligand$receptor.complex, levels = all_receptors)
  top_100_ligand$count = rep(1, nrow(top_100_ligand))

  # Plot top ligands and their interaction partners (p3)
  p3 <- top_100_ligand %>%
    ggplot(aes(x = ligand.complex, y = count, fill = receptor.complex)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = receptor_colors) + 
    theme_bw() +
    coord_flip() +
    labs(title = paste0("Top Ligands Partners (", sample, ")")) + 
    theme(
      legend.position = "right", 
      legend.title = element_blank(),
      legend.text =  element_text(size = 20, color = "black"),
      axis.text.y.left = element_text(size = 20, color = "black"),
      axis.title.y.left = element_blank(),
      axis.text.x.bottom = element_text(size = 20, color = "black"),
      axis.title.x.bottom = element_text(size = 20, color = "black")
    )
  
  # Store the plot in the list
  plots[[index]] <- p3
}

# Combine the two plots into one figure, arranged in one row
combined_plot <- ggpubr::ggarrange(plots[[1]], plots[[3]], ncol = 2, common.legend = FALSE, legend = "right")

# Save the combined plot
ggsave(file.path(plot_result_path, "combined_top_100_Ligand_Partner_barplot_legend.pdf"), plot = combined_plot, width = 22, height = 10)


plots[[1]]
ggsave(file.path(plot_result_path, "Untreated_combined_top_100_Ligand_Partner_barplot_legend.pdf"), width = 9, height = 8)

plots[[3]]
ggsave(file.path(plot_result_path, "Vivo_combined_top_100_Ligand_Partner_barplot_legend.pdf"), width = 10, height = 8)

### CELL CHAT
#functions
cellchat_vidualization<-function(sample_name,topn){

  context_df_dict_top <- lapply(context_df_dict,function(df){head(df,topn)}) 
  sample<-context_df_dict_top[[sample_name]]
count_data <- sample %>%
  group_by(source, target) %>%
  tally()

wider_data <- count_data %>%
  pivot_wider(names_from = target, values_from = n, values_fill = 0)


rowname<- paste0("c",wider_data$source)
wider_data<-wider_data[,2:ncol(wider_data)]

rownames(wider_data)<-rowname
colnames(wider_data)<-as.factor(paste0("c",colnames(wider_data)))

#add missing cluster
levels<-c("c1","c2","c3","c4","c5","c6","c7","c8","c9")
colors <- custom_palette

new_row=rep(0,ncol(wider_data))
for (i in levels){
  if (!(i %in% rownames(wider_data))){
    wider_data=rbind(wider_data,new_row)
    rownames(wider_data)[nrow(wider_data)]=i
  }
}
rowname=rownames(wider_data)

new_col=rep(0,nrow(wider_data))
for (i in levels){
  if (!(i %in% colnames(wider_data))){
    wider_data=cbind(wider_data,new_col)
  colnames(wider_data)[ncol(wider_data)]=i
  }
}

rownames(wider_data)<-rowname

groupSize <- cluster_sizes %>% filter(orig.ident == sample_name)%>%pull(count)
names(groupSize)<-levels
matrix<-as.matrix(wider_data)

matrix<-matrix[levels,levels]

pdf(file.path(plot_result_path, paste0(sample_name,"_circle_chart_top_",topn,".pdf")),width = 4,height = 4)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(matrix, 
                color.use = custom_palette,
                vertex.weight = groupSize, 
                vertex.label.cex = 1,
                vertex.label.color  = 'black' , 
                weight.scale = T, 
                label.edge= F, 
                edge.width.max = 7,
                title.name = sample_name)
dev.off()

write.table(matrix, file.path(Figures_data,paste0( "Figure_6F_circle_plot_top_50_interaction_", sample_name, ".txt")), sep = "\t", quote =F )


#Buble Plot

sample<-sample%>%mutate(interaction=paste0(ligand.complex,"-",receptor.complex))

#for senders
sender<-sample%>%select(source,aggregate_rank,interaction)%>%group_by(source,interaction)%>%
  summarise(count=n(),avg_aggre_rank=1/mean(aggregate_rank))

#for receiver
receiver<-sample%>%select(target,aggregate_rank,interaction)%>%group_by(target,interaction)%>%
  summarise(count=n(),avg_aggre_rank=1/mean(aggregate_rank))

lower_count<-min(range(sender$count)[1],range(receiver$count)[1])
higher_count<-max(range(sender$count)[2],range(receiver$count)[2])
lower<-min(range(sender$avg_aggre_rank)[1],range(receiver$avg_aggre_rank)[1])
higher<-max(range(sender$avg_aggre_rank)[2],range(receiver$avg_aggre_rank)[2])

p1=ggplot(sender, aes(x = source, y =interaction, size = count, color = avg_aggre_rank)) +
  geom_point() +
  scale_size_continuous(limits = c(lower_count,higher_count),breaks=seq(lower_count,higher_count,1))+
  scale_x_discrete()+
  labs(title = paste0(sample_name,"_Senders"),
       x = "Sender",
       y = "Interaction",
       size = "count",
       color = "significant") +
  scale_color_continuous(
    breaks = c(lower,higher),#ne the custom breaks for color legend
    labels = c("Low",  "High"),  # Optional labels for the breaks
    limits = c(lower,higher))+# the color scale limits
  theme_minimal()

  p2=ggplot(receiver, aes(x = target, y =interaction, size = count, color = avg_aggre_rank)) +
    geom_point() +
    scale_size_continuous(limits = c(lower_count,higher_count),breaks=seq(lower_count,higher_count,1))+
    scale_x_discrete()+
    labs(title = paste0(sample_name,"_Receivers"),
         x = "Receiver",
         y = "Interaction",
         size = "count",
         color = "significant") +
    scale_color_continuous(
      breaks = c(lower,higher),#ne the custom breaks for color legend
      labels = c("Low",  "High"),  # Optional labels for the breaks
      limits = c(lower,higher))+# the color scale limits
    theme_minimal()
ggarrange(p1,p2,ncols=2,common.legend = T,label.y = "interactions")
# ggsave(filename=paste0(sample_name,"_buble_chart_top_",topn,".pdf"),path=plot_result_path,width = 18000,height = 15000,dpi=1500,units="px")
}


purrr::map(names(context_df_dict),~cellchat_vidualization(.,topn= 50))

