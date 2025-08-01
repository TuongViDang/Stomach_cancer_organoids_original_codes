library(tidyverse)

#Figures data path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)



SCENIC_dir = "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCENIC/"
setwd(SCENIC_dir )

# Read while skipping the first row
regulons <- read.csv("results/regulons.csv", skip=1, stringsAsFactors=FALSE, header=TRUE)

regulons[1,]
# Remove the first row (header within data)
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


# Subset for TFs of interest
keywords <- c("NFKB", "REL")


# Use grepl to find matches in the 'TF' column
NFKB_regulons <- regulons_clean[grepl(paste(keywords, collapse="|"), regulons_clean$TF, ignore.case=TRUE), ]
NFKB1_regulons <- regulons_clean[grepl("NFKB1", regulons_clean$TF, ignore.case=TRUE), ]


# View
head(NFKB_regulons)


# Find which regulons contain CXCL1 or CXCL3
cytokines <- c("CXCL1", "CXCL2", "CXCL3", "CXCL8")
NFKB_cytokines_regulons <- NFKB_regulons[grepl(paste(cytokines, collapse="|"), NFKB_regulons$TargetGenes), ]

NFKB1_cytokines_regulons <- NFKB1_regulons[grepl(paste(cytokines, collapse="|"), NFKB1_regulons$TargetGenes), ]
NFKB1_cytokines_regulons_direct <- NFKB1_cytokines_regulons %>%
  filter(Annotation == "gene is directly annotated")


#NFKB2_cytokines_regulons <- NFKB2_regulons[grepl(paste(cytokines, collapse="|"), NFKB2_regulons$TargetGenes), ]
#NFKB2_cytokines_regulons_direct <- NFKB2_cytokines_regulons %>%
#  filter(Annotation == "gene is directly annotated")


TF_kept = NFKB_cytokines_regulons$TF %>% unique()


# Read AUC matrix
auc_matrix <- read.csv("results/auc_mtx.csv")
colnames(auc_matrix) = sub("\\..*","", colnames(auc_matrix))


#only keep NFKB 
auc_matrix <- auc_matrix[,TF_kept]

write.table(auc_matrix,"results/auc_mtx_NFKB.tsv", sep = "\t", quote = F)

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


regulon_df$highlight <- ifelse(rank(-regulon_df$Weight) <= 4, "Top 4", "Others")

plot <- ggplot(regulon_df, aes(x = reorder(Gene, Weight), y = Weight, fill = highlight)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "NFKB1 Regulon",
       x = "Gene",
       y = "Weight") +
  scale_fill_manual(values = c("Top 4" = "red", "Others" = "grey")) +
  theme_classic()+
  theme(legend.position = "None",
    plot.title = element_text(size = 33, color = "black" ),
    axis.title.y.left = element_text(size = 33, color = "black" ),
    axis.text.y.left = element_text(size = 33,  color = "black"),
    axis.text.x.bottom = element_text(size = 33,  color = "black"),
    axis.title.x.bottom = element_text(size = 33,  color = "black"))
plot 
ggsave("NFKB1_regulons.pdf", width= 13.8, height = 15, scale =1)

data <- plot$data
head(data)
write.table(data, file.path(Figures_data,"Figure_6A_NFKB1_regulon.txt" ), sep = "\t", quote = F)



#NFKB2
json_string <- NFKB2_cytokines_regulons_direct$TargetGenes[1]

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


regulon_df$highlight <- ifelse(rank(-regulon_df$Weight) <= 4, "Top 4", "Others")

 ggplot(regulon_df, aes(x = reorder(Gene, Weight), y = Weight, fill = highlight)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "NFKB2 Regulon",
       x = "Gene",
       y = "Weight") +
  scale_fill_manual(values = c("Top 4" = "red", "Others" = "grey")) +
  theme_classic()+
  theme(legend.position = "None",
    plot.title = element_text(size = 20, color = "black" ),
    axis.title.y.left = element_text(size = 20, color = "black" ),
    axis.text.y.left = element_text(size = 20,  color = "black"),
    axis.text.x.bottom = element_text(size = 20,  color = "black"),
    axis.title.x.bottom = element_text(size = 20,  color = "black"))
ggsave("NFKB2_regulons.pdf", width= 8, height = 9.5)
