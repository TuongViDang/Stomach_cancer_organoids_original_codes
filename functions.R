##READ DATA
read_scRNAseq <- function(data_path, patients, min.cells = 3, min.features = 200) {
  all_patients <- vector("list", length(patients))
  
  for (i in seq_along(patients)) {
    counts <- Read10X(data.dir = file.path(data_path, patients[i]))   
    all_patients[[i]] <- CreateSeuratObject(counts = counts, project = patients[i], min.cells = min.cells, min.features = min.features)
  }
  
  return(all_patients)
}

#QC
#### CALCULATE QC METRICS ####

qc_scRNAseq <- function(object){
  
  # number of genes detected per UMI: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
  object$log10GenesPerUMI <- log10(object$nFeature_RNA) / log10(object$nCount_RNA)
  
  # mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes
  object <- PercentageFeatureSet(object = object, pattern = "^MT-", col.name = 'mitoRatio')
  object <- PercentageFeatureSet(object = object, pattern = "^RP[SL]", col.name = 'riboRatio')
  
  object$mitoRatio <- object$mitoRatio / 100
  object$riboRatio <- object$riboRatio / 100
  
  return(object)
  
}


###VIOLIN PLOTS###

plot_violin_scRNAseq<-function(object,feature,y_axis_name,plot_path,plot_name,width,height,dpi,units,device
){
  library(ggplot2)
  library(cowplot)
  
  vlnplots <-lapply(seq_along(features),function(i){
    VlnPlot(object = object, group.by = "orig.ident", features = features[i], pt.size = 0.00001) +
      ylab(y_axis_name[i])+
      theme(legend.position = 'none',
            axis.title.x = element_blank())
  })
  
  vln_qc <- plot_grid(plotlist = vlnplots, labels = '', ncol = 5)
  vln_com<-VlnPlot(object,group.by = 'orig.ident',features='log10GenesPerUMI',pt.size=0.00001)+
    ylab('complexity')+
    theme(legend.position = 'none',
          axis.title.x = element_blank())
  vln_com_grid<-plot_grid(NULL,vln_com,NULL,nrow=1,rel_widths = c(1,2,1))
  vln_qc_final<-plot_grid(vln_qc,vln_com_grid,nrow=2,rel_heights = c(2,1),rel_widths = c(2,1))
  
  ggsave(filename =plot_name,path = plot_path,plot = vln_qc_final, width = width, height = height, units = units, dpi = dpi, device = device)
}

#NUMBER OF CELLS PLOT

plot_number_cells<-function(metadata,plot_name,plot_path,width = width, height = height, units = units, dpi = dpi){
  theme_set(theme_classic())
  
  ncells <- ggplot(metadata, aes(x = orig.ident, fill = orig.ident)) +
    geom_bar(show.legend = FALSE) +
    geom_text(stat='count',aes(label=stat(count)),vjust=-0.5)+
    scale_y_continuous(breaks = seq(0, 5000, 1000)) +
    ylab("Number of Cells")+
    labs(title='Number of cells in each organoid')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 16)) 
  
  
  ggsave(filename = plot_name,path=plot_path, plot = ncells,width = width, height = height, units = units, limitsize = FALSE, dpi = dpi)  
}


#DENSITY PLOT
Plot_density_qc<-function(metadata, umi_cut, gene_cut, mt_cut, log10_cut 
                          , plot_name , plot_path, width , height , dpi, units
){
  plot_density <- function(data, x, xlab, vline = NULL) {
    ggplot(data, aes(color = orig.ident,  x = !!rlang::enquo(x), fill = orig.ident)) +
      geom_density(alpha = 0.2) +
      xlab(xlab) +
      ylab("Cell density") +
      geom_vline(xintercept = vline,linetype='dashed') +
      theme(legend.title = element_blank())
    #+geom_text(aes(x = vline, y = 0, label = vline), vjust = 1.5, hjust = -0.2) 
    #+theme(legend.position = 'top',legend.title = element_blank())
  }
  
  
  umi_transcript <- plot_density(metadata, x= nCount_RNA, xlab = 'UMI count',vline =  umi_cut)
  #ggsave(filename = 'UMI_count.pdf',path=plot_path_before_qc, plot = umi_transcript,  width = 9600, height = 9600, units = 'px', limitsize = FALSE, dpi = 1500)
  
  genes_per_cells <- plot_density(metadata, x=nFeature_RNA, xlab='Gene count', vline= gene_cut)
  #ggsave(filename = 'Feature_count.pdf', plot = genes_per_cells, path = plot_path_before_qc, width = 9600, height = 9600, units = 'px', limitsize = FALSE, dpi = 1500)
  
  mt_genes <- plot_density(metadata, x=mitoRatio, xlab = 'Proportion of mitochondrial reads',vline=mt_cut)
  #ggsave(filename = 'Mitochondria_proportion', plot = mt_genes, path = plot_path_before_qc, width = 9600, height = 9600, units = 'px', limitsize = FALSE, dpi = 1500)
  
  complexity <- plot_density(metadata, log10GenesPerUMI, 'Complexity',  log10_cut)
  #ggsave(filename = 'complexity.pdf', plot = complexity, path = plot_path_before_qc, width = 9600, height = 9600, units = 'px', limitsize = FALSE, dpi = 1500)
  
  r_genes <- plot_density(metadata, riboRatio, 'Proportion of ribosomal reads')
  #ggsave(filename = paste0('r_genes_', plot_name), plot = r_genes, path = plot_path_before_qc, width = 9600, height = 9600, units = 'px', limitsize = FALSE, dpi = 1500)
  
  library(ggpubr)
  distr_qc <- ggarrange(umi_transcript, genes_per_cells,mt_genes,r_genes,common.legend = T)
  ggsave(filename=plot_name,path=plot_path,plot=distr_qc,width = width, height = height, units = units, dpi = dpi)
}


#Plot contribution cell cycle phase for each patient
plot_patients_contribution_phase <- function (srat) { 
## take an integrated Seurat object, plot distributions over orig.ident
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

count_table <- table(srat@meta.data$Phase, srat@meta.data$orig.ident)
count_mtx   <- as.data.frame.matrix(count_table)
count_mtx$Phase <- rownames(count_mtx)
melt_mtx    <- melt(count_mtx)
colnames(melt_mtx)[2] <- "patients"
melt_mtx$patients <- as.factor(melt_mtx$patients)


Group_size   <- aggregate(value ~ patients, data = melt_mtx, FUN = sum)
sorted_labels<-sort(unique(melt_mtx$patients),decreasing = T)
Group_size$patients<-factor(Group_size$patients,levels = sorted_labels)
melt_mtx$patients <- factor(melt_mtx$patients,levels = sorted_labels)
melt_mtx$Phase<-factor(melt_mtx$Phase,levels = c('G1','S','G2M'))

p1 <- ggplot(Group_size, aes(y= patients,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
  theme_bw() + scale_x_log10() + xlab("Cells per samples, log10 scale") + ylab("")
p2 <- ggplot(melt_mtx,aes(x=patients,y=value,fill=Phase)) + 
  geom_bar(position="fill", stat="identity") + scale_x_discrete(labels= c("OO100-no"="OO100 untreated","OO100-vivo"="OO100 in vivo","OO99-no"="OO99 untreated","OO99-vivo"="OO99 in vivo","OO77-no"="OO77 untreated","OO77-vivo"="OO77 in vivo")) +
  coord_flip() + 
  scale_fill_brewer(palette = "Set2") +
  guides(fill = guide_legend(reverse = TRUE))+
  ylab("Fraction of cells in each sample") + xlab("") + theme_minimal() +theme(legend.position="top", legend.title=element_blank())

p2 + p1 + plot_layout(widths = c(5,1))
}



#Plot cluster contribution  for each patient
plot_cluster_contribution_each_sample <- function (srat) { 
## take an integrated Seurat object, plot distributions over orig.ident
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$Sample)
count_mtx   <- as.data.frame.matrix(count_table)
count_mtx$seurat_clusters <- rownames(count_mtx)
melt_mtx    <- melt(count_mtx)
colnames(melt_mtx)[2] <- "Sample"

ordered_samples = rev(c("OO77 Untreated",  "OO77 Ex vivo", "OO77 In vivo",
                   "OO99 Untreated", "OO99 Ex vivo", "OO99 In vivo",
                  "OO100 Untreated", "OO100 Ex vivo", "OO100 In vivo" ))

melt_mtx$Sample <- factor(melt_mtx$Sample, levels =ordered_samples )
melt_mtx$seurat_clusters <-  factor(melt_mtx$seurat_clusters, levels = rev(1:9))

Group_size   <- aggregate(value ~ Sample, data = melt_mtx, FUN = sum)


p1 <- ggplot(Group_size, aes(y= Sample ,x = value)) + 
  geom_bar(position="dodge", stat="identity",fill = "grey60", width = 0.6)  + 
  theme_bw() + 
  scale_x_log10() + 
  xlab("Cells per sample, log10") + 
  ylab("") +
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15, color = "black"),
         axis.title.y.left = element_blank(),
         axis.text.y.left = element_blank(),
         axis.title.x.bottom = element_text(size = 15, color = "black"),
         axis.text.x.bottom = element_text(size = 15, color = "black"))
  

# Define the color palettes
custom_palette <- c(  "#80B1D3", "limegreen", "#FCCDE5" ,
                  "#FDB462","#D9D9D9", "#CCEBC5" ,
                 "#BC80BD" , "khaki1",  "brown1")   
names(custom_palette ) <- 1:9              

p2 <- ggplot(melt_mtx,aes(x= Sample,y=value,fill= seurat_clusters)) + 
  geom_bar(position="fill", stat="identity", width = 0.6) + theme_bw() + coord_flip() + 
  scale_fill_manual(values = custom_palette) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1))+
  ylab("Fraction of cells in each sample") + xlab("Sample") + 
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15, color = "black"),
        axis.title.y.left = element_text(size = 15, color = "black"),
         axis.text.y.left = element_text(size = 15, color = "black"),
         axis.title.x.bottom = element_text(size = 15, color = "black"),
         axis.text.x.bottom = element_text(size = 15, color = "black"))

p2 + p1 + plot_layout(widths = c(3.5,1.5))
}










#Plot contribution each patient to each cluster



plot_integrated_clusters = function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$Sample)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "Sample"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) +
       geom_bar(position="dodge", stat="identity",fill = "grey60") + 
       theme_bw() + scale_x_log10() + 
       xlab("Cells per cluster, log10 scale") + 
       ylab("")+
       theme(axis.title.y.left = element_text(size = 15, color = "black"),
             axis.text.y.left = element_text(size = 15, color = "black"),
             axis.title.x.bottom = element_text(size = 15, color = "black"),
             axis.text.x.bottom = element_text(size = 15, color = "black"))

  
  #change the name of sample
  melt_mtx <- melt_mtx %>%
  mutate(Sample = case_when(
    grepl("-no", Sample ) ~ sub("-no", " Untreated", Sample),
    grepl("-vitro", Sample) ~ sub("-vitro", " Ex vivo", Sample),
    grepl("-vivo", Sample) ~ sub("-vivo", " In vivo", Sample),
    TRUE ~ Sample
  ))
  melt_mtx$Sample <- as.factor(  melt_mtx$Sample )
  #color
  color_palette <- c(brewer.pal(12, "Set3")[1:8],  brewer.pal(8, "Set2")[4])
  names(  color_palette  ) <- levels( melt_mtx$Sample )
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill= factor(Sample, levels = rev(levels(Sample)) ))) + 
    scale_fill_manual(values =  color_palette )+
    guides(fill = guide_legend(nrow = 3, reverse = TRUE)) +
    geom_bar(position="fill", stat="identity") +
    theme_bw() + 
    coord_flip() + 
    ylab("Fraction of cells in each sample") +
    xlab("Cluster number") + 
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.text = element_text(size = 15, color = "black"),
          axis.title.y.left = element_text(size = 15, color = "black"),
          axis.text.y.left = element_text(size = 15, color = "black"),
          axis.title.x.bottom = element_text(size = 15, color = "black"),
          axis.text.x.bottom = element_text(size = 15, color = "black"))

   
  p2 + p1 + plot_layout(widths = c(3,1))

  
}




#Plot contribution each Phase to each cluster

plot_integrated_clusters_Phase = function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$Phase)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    #scale_fill_brewer(palette = "Set2") +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  p2 + p1 + plot_layout(widths = c(3,1))
  
  
}




