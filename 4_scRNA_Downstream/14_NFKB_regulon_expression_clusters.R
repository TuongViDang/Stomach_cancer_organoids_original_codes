
## Load libraries
library(Seurat)
library(tidyverse)
library(HGNChelper)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(RColorBrewer)


#Specify the path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_path,"NFKB1_regulon_expression_clusters")
dir.create(plot_result_path)

#Load file
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))
DefaultAssay(seurat_integrated) <- "RNA"
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
seurat_integrated <- SCTransform(seurat_integrated, 
                                 vst.flavor = "v1",
                                  variable.features.n = nrow(seurat_integrated ))

SCENIC_dir = "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCENIC/"
setwd(SCENIC_dir )
regulons <- read.csv("results/regulons.csv", stringsAsFactors = FALSE)
regulons_clean <- regulons[-1, ]
# Rename the columns for clarity
colnames(regulons_clean) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue",  "OrthologousIdentity"  ,
                              "Annotation"   ,         "Context"       ,        "TargetGenes"    , "RankAtMax" )

# Convert numeric columns (AUC, NES, MotifSimilarityQvalue) from character to numeric
regulons_clean$AUC <- as.numeric(regulons_clean$AUC)
regulons_clean$NES <- as.numeric(regulons_clean$NES)
regulons_clean$MotifSimilarityQvalue <- as.numeric(regulons_clean$MotifSimilarityQvalue)

# View cleaned dataframe
head(regulons_clean)
regulons_clean = regulons_clean  %>% mutate( ID = paste0(TF, "_", MotifID))

#NFKB1 regulon
NFKB1_regulons <- regulons_clean[grepl("NFKB1", regulons_clean$TF, ignore.case=TRUE), ]

cytokines <- c("CXCL1", "CXCL2", "CXCL3", "CXCL8")
NFKB1_cytokines_regulons <- NFKB1_regulons[grepl(paste(cytokines, collapse="|"), NFKB1_regulons$TargetGenes), ]
NFKB1_cytokines_regulons_direct <- NFKB1_cytokines_regulons %>%
  filter(Annotation == "gene is directly annotated")

#NFKB1
json_string <- NFKB1_cytokines_regulons_direct$TargetGenes[2]

# Remove leading and trailing square brackets
json_clean <- gsub("^\\[|\\]$", "", json_string)

# Replace (' and ') with c(" and ")
json_clean <- gsub("\\('", 'c("', json_clean)
json_clean <- gsub("'\\)", '")', json_clean)

# Replace ', ' separator with ', '
json_clean <- gsub("', ", '", ', json_clean)
# Evaluate each element
regulon_list <- eval(parse(text = paste0("list(", json_clean, ")")))

# Convert to data frame
regulon_df <- do.call(rbind, lapply(regulon_list, function(x) data.frame(Gene = x[1], Weight = as.numeric(x[2]))))
rownames(regulon_df) <- NULL

regulon_df <- regulon_df %>% arrange(Weight)

plot <- DotPlot(seurat_integrated,features= regulon_df$Gene) + coord_flip()
plot
ggsave(filename="Dotplot_NFKB1_regulon.pdf",path=plot_result_path,width= 8,height= 8)

data = plot$data

# Bar plot of Log2FC
data_wide <- data %>%
  pivot_wider(
    id_cols = features.plot,
    names_from = id,
    values_from = avg.exp,
    values_fill = 0
  ) %>%
  rename(gene = features.plot)

data_wide = data_wide %>%  mutate(others = rowMeans(select(., `1`:`8`)),
                                  Log2FC = log2(`9`/others)) %>%
                           arrange(desc(Log2FC))


data_wide %>% ggplot(aes(x = reorder(gene, Log2FC), y = Log2FC, fill = Log2FC)) + 
              geom_bar(stat = "identity") +
              coord_flip() +
  labs(
       x = "Gene",
       y = "Log2FC(Cluster 9/Others)") +
  scale_fill_gradient2(
    low = "white",     # low values
    high = "red",     # high values
    midpoint = 0      # centered at zero for log fold-change
  ) +
  theme_classic()+
  theme(legend.position = "None",
    plot.title = element_text(size = 33, color = "black" ),
    axis.title.y.left = element_text(size = 33, color = "black" ),
    axis.text.y.left = element_text(size = 33,  color = "black"),
    axis.text.x.bottom = element_text(size = 33,  color = "black"),
    axis.title.x.bottom = element_text(size = 33,  color = "black"))
ggsave(filename="Barplot_NFKB1_regulon_expression_clusters.pdf",path=plot_result_path,width= 13.8,height= 15)   


#reorder the Dotplot with the same gene order as in the barplot
gene_order = rev(data_wide$gene)
DotPlot(seurat_integrated,features= gene_order ) + coord_flip()
ggsave(filename="Dotplot_NFKB1_regulon_reordered.pdf",path=plot_result_path,width= 8,height= 8)
