library(tidyverse)
library(ggpubr)
library(vegan)


CONIPHER_DIR="/group/poetsch_projects/poetsch_sc/Conifer"
Integrate_dir="/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/ITH"
plot_dir=file.path(Integrate_dir,"plot")
dir.create(plot_dir, recursive = T)


#Figures data path
main_path <-'/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)

#####  VI's METHOD ###### 
#JACCARD INDEX
#Number of clones at the snapshot
#OO77
# Clone proportion in scRNA-seq
OO77_sc_proportion_dir=file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO77","tree30_change_name","txt_files")
OO77_clone_pro = read.delim(file.path(OO77_sc_proportion_dir, "Clone_proportion.txt"))
OO77_clone_pro  = OO77_clone_pro %>% mutate(sample = case_when(sample == "OO77-no" ~ "OO77_tumor",
                                                              sample == "OO77-vitro" ~ "OO77_ex",
                                                              sample == "OO77-vivo" ~ "OO77_in"))

OO77_clone_pro = OO77_clone_pro %>% pivot_wider(names_from = "sample", values_from = "proportion",  values_fill = list(proportion = 0))
colnames(OO77_clone_pro)[1] = "clusterID"

#conifer tree structure
OO77_tree_df = read.delim(file.path(CONIPHER_DIR,"OO77/output/tree_30_final/Trees/allTrees.txt"),comment.char = "#",header=F)
colnames(OO77_tree_df) = c("parent","clusterID")

#we need to change cluster name
# Define mapping for replacements. Here make a trick: change 6 to 5a, 2 to 5b
clusterID_map = c(`10` = "1", `1` = "3", `6` = "5a", `2` = "5b", `5` = "6", `11` = "7", `9` = 8, `3` = 2)

#change clone name 
OO77_tree_df = OO77_tree_df %>%
  mutate(
    clusterID = ifelse(clusterID %in% names(clusterID_map), clusterID_map[as.character(clusterID)], clusterID),
    parent = ifelse(parent %in% names(clusterID_map), clusterID_map[as.character(parent)], parent)
  )
table(OO77_tree_df$clusterID)




#OO99
# Clone proportion in scRNA-seq
OO99_sc_proportion_dir=file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO99","tree30_change_name","txt_files")
OO99_clone_pro = read.delim(file.path(OO99_sc_proportion_dir, "Clone_proportion.txt"))
OO99_clone_pro  = OO99_clone_pro %>% mutate(sample = case_when(sample == "OO99-no" ~ "OO99_tumor",
                                                              sample == "OO99-vitro" ~ "OO99_ex",
                                                              sample == "OO99-vivo" ~ "OO99_in"))
OO99_clone_pro = OO99_clone_pro %>% pivot_wider(names_from = "sample", values_from = "proportion",  values_fill = list(proportion = 0))
colnames(OO99_clone_pro)[1] = "clusterID"


#conifer tree structure
OO99_tree_df = read.delim(file.path(CONIPHER_DIR,"OO99/output/tree_30/Trees/allTrees.txt"),comment.char = "#",header=F)
colnames(OO99_tree_df) = c("parent","clusterID")

#we need to change cluster name
clusterID_map <- c(`30` = 2, `8` = 3, `2` = 4, `4` = 5, `12` = 6, `15` = 7, 
                   `13` = 8, `23` = 9, `29` = 10, `5` = 11, `14` = 12)

#change clone name 
OO99_tree_df = OO99_tree_df %>%
  mutate(
    clusterID = ifelse(clusterID %in% names(clusterID_map), clusterID_map[as.character(clusterID)], clusterID),
    parent = ifelse(parent %in% names(clusterID_map), clusterID_map[as.character(parent)], parent)
  )
table(OO99_tree_df$clusterID)

#OO100
# Clone proportion in scRNA-seq
OO100_sc_proportion_dir=file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA/OO100","tree30_change_name","txt_files")
OO100_clone_pro = read.delim(file.path(OO100_sc_proportion_dir, "Clone_proportion.txt"))
OO100_clone_pro  = OO100_clone_pro %>% mutate(sample = case_when(sample == "OO100-no" ~ "OO100_tumor",
                                                              sample == "OO100-vitro" ~ "OO100_ex",
                                                              sample == "OO100-vivo" ~ "OO100_in"))
OO100_clone_pro = OO100_clone_pro %>% pivot_wider(names_from = "sample", values_from = "proportion",  values_fill = list(proportion = 0))
colnames(OO100_clone_pro)[1] = "clusterID"


#conifer tree structure
OO100_tree_df = read.delim(file.path(CONIPHER_DIR,"OO100/output/tree_30_final/Trees/allTrees.txt"),comment.char = "#",header=F)
colnames(OO100_tree_df) = c("parent","clusterID")


#we need to change cluster name
clusterID_map <- c(`14` = 2, `5` = 3, `7` = 4, `23` = 5, `4` = 6, `18` = 7, 
                   `2` = 8, `19` = 9, `22` = 10, `6` = 11, `11` = 12)

#change clone name 
OO100_tree_df = OO100_tree_df %>%
  mutate(
    clusterID = ifelse(clusterID %in% names(clusterID_map), clusterID_map[as.character(clusterID)], clusterID),
    parent = ifelse(parent %in% names(clusterID_map), clusterID_map[as.character(parent)], parent)
  )
table(OO100_tree_df$clusterID)


# Function to find the lineage of a clone
find_lineage <- function(clusterID, tree_df) {
  lineage <- c(clusterID)
  while(TRUE) {
    parent <- tree_df$parent[tree_df$clusterID == clusterID]
    if(length(parent) == 0) break
    lineage <- c(lineage, parent)
    clusterID <- parent
  }
  return(lineage)
}

# Function to find the most recent common ancestor (MRCA)
find_mrca <- function(clone1, clone2, tree_df) {
  lineage1 <- find_lineage(clone1, tree_df)
  lineage2 <- find_lineage(clone2, tree_df)
  # Find the common ancestor by intersecting the two lineages
  mrca <- intersect(lineage1, lineage2)[1] # The first common node
  return(mrca)
}

# Function to calculate the distance to the MRCA
#based on tree step
calculate_distance_tree <- function(clone1, clone2, tree_df) {
  mrca <- find_mrca(clone1, clone2, tree_df)
  lineage1 <- find_lineage(clone1, tree_df)
  lineage2 <- find_lineage(clone2, tree_df)
  # Evolutionary distance is the number of steps to MRCA normalized by 100*total edges  
  all_nodes <- unique(c(tree_df$parent, tree_df$clusterID))
  num_edges <- length(all_nodes)-1
  distance1 <- (which(lineage1 == mrca) - 1)/num_edges
  distance2 <- (which(lineage2 == mrca) - 1)/num_edges
  total_distance <- (distance1 + distance2)
  return(total_distance)
}
# Function to calculate the weighted distance between two clones
calculate_weighted_distance <- function(clone1, clone2, clone_proportion, tree_df, calculate_distance_func) {
  distance <- calculate_distance_func(clone1, clone2, tree_df)
  proportion1 <- clone_proportion$proportion[clone_proportion$clusterID == clone1]
  proportion2 <- clone_proportion$proportion[clone_proportion$clusterID == clone2]
  #Find lineage
  lineage1 <- find_lineage(clone1, tree_df)
  lineage2 <- find_lineage(clone2, tree_df)
  # Weighted distance as distance times the sum of the proportions
  weighted_distance <- distance * (proportion1 + proportion2 ) 
  return(weighted_distance)
}



# Define the function to calculate ITH for a given sample
calculate_ITH <- function(sample,patient, clone_proportion_df, tree_df, calculate_distance_func) {
  # Filter the DataFrame based on the sample
  df  <- clone_proportion_df %>%
    select(!!rlang::sym(sample), clusterID) %>%
    filter(!!rlang::sym(sample) > 0)
  
  colnames(df)[1] <- "proportion"
  clone_ids <- df$clusterID
  n <- length(clone_ids)
  
  # Initialize the distance matrix
  distance_matrix <- matrix(0, n, n)
  rownames(distance_matrix) <- clone_ids
  colnames(distance_matrix) <- clone_ids
  
  # Compute the distance matrix
  for (i in 1:n) {
    for (j in 1:n) {
      clone1 <- clone_ids[i]
      clone2 <- clone_ids[j]
      distance_matrix[i, j] <- calculate_weighted_distance(clone1, clone2, clone_proportion = df, tree_df = tree_df, calculate_distance_func = calculate_distance_func)
    }
  }
  
  # Calculate ITH: Sum of the upper triangle without diagonal
  ITH <- sum(distance_matrix[upper.tri(distance_matrix, diag = FALSE)])/ n
  
  ITH_df = data.frame(patient = patient, Treatment = str_split(sample,"_")[[1]][2], ITH = ITH)
  ITH_df
}

samples = c("OO77_tumor","OO77_ex", "OO77_in")
OO77_ITH = purrr::map_dfr(samples,~calculate_ITH(.x, patient= "OO77",  clone_proportion_df = OO77_clone_pro, tree_df = OO77_tree_df, calculate_distance_func = calculate_distance_tree) )


samples = c("OO99_tumor","OO99_ex", "OO99_in")
OO99_ITH = purrr::map_dfr(samples,~calculate_ITH(.x, patient= "OO99",  clone_proportion_df = OO99_clone_pro, tree_df = OO99_tree_df, calculate_distance_func = calculate_distance_tree) )


samples = c("OO100_tumor","OO100_ex", "OO100_in")
OO100_ITH = purrr::map_dfr(samples,~calculate_ITH(.x, patient= "OO100",  clone_proportion_df = OO100_clone_pro, tree_df = OO100_tree_df ,  calculate_distance_func = calculate_distance_tree) )


#Weghted Shannon entropy

calculate_Entropy <- function(sample,patient, clone_proportion_df, tree_df){
df  <- clone_proportion_df %>%
    select(!!sym(sample), clusterID) %>%
    filter(!!sym(sample) > 0)
    colnames(df)[1] <- "proportion"

df <- df %>% mutate(Entropy = proportion/100 * log(proportion/100)  )
tumor_Entropy <- - sum(df$Entropy)
result = data.frame(patient = patient, Treatment = str_split(sample,"_")[[1]][2], ITH = tumor_Entropy)
result
}   

samples = c("OO77_tumor","OO77_ex", "OO77_in")
OO77_Entropy = purrr::map_dfr(samples, ~calculate_Entropy(.x, patient = "OO77", clone_proportion_df = OO77_clone_pro, tree_df = OO77_tree_df))

samples = c("OO99_tumor","OO99_ex", "OO99_in")
OO99_Entropy = purrr::map_dfr(samples, ~calculate_Entropy(.x, patient = "OO99", clone_proportion_df = OO99_clone_pro, tree_df = OO99_tree_df))

samples = c("OO100_tumor","OO100_ex", "OO100_in")
OO100_Entropy = purrr::map_dfr(samples, ~calculate_Entropy(.x, patient = "OO100", clone_proportion_df = OO100_clone_pro, tree_df = OO100_tree_df))

All = rbind(OO77_Entropy, OO99_Entropy, OO100_Entropy)

All$patient = factor(All$patient, levels = c("OO77", "OO99", "OO100"))
# Create a new column for labeling only tumor samples
#All$label <- ifelse(All$Treatment == "tumor", All$patient, NA)
#Ex vivo
All_ex <- All %>% filter(Treatment != "in")
All_ex$Treatment = factor(All_ex$Treatment, levels = c("tumor","ex"))

#t test
t_test <- t.test(ITH ~ Treatment, data = All_ex, paired = TRUE, alternative = "greater")
p_value <-  t_test$p.value



colors = c( "deepskyblue2","tomato","goldenrod3", "darkcyan","indianred1")

#plot
ggplot(All_ex, aes(x = Treatment, y =  ITH , group = patient)) +
  geom_point(aes(colour = Treatment),  shape = 16,size = 2, position = position_dodge(width = 0)) +
  geom_line(aes(colour = patient), linewidth = 2, alpha = 0.5, position = position_dodge(width = 0)) +
  scale_y_continuous(limits = c(0,max(All_ex$ITH + 1)),breaks = c(0,1,2,3))+
  scale_x_discrete(labels = c("ex"="Ex vivo","tumor" = "Untreated"))+
  scale_colour_manual(values = colors)+
  ylab("genetic ITH") +
  #geom_text(aes(label =  label), vjust = -0.5, position = position_dodge(width = 0.5), size = 6) +
  theme_classic() +
  theme( legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
         panel.spacing = unit(1, "lines"))  +
   # Add the p-value as annotation
  annotate("text", x = 1.5, y = 2.75, label = paste("p =", format.pval(p_value, digits = 3)), size = 8)

  ggsave("genetic_ITH_ex_vivo_Entropy.pdf",path=plot_dir,width= 650 ,height=650,dpi=150,units="px")

#In vivo
All_in <- All %>% filter(Treatment != "ex")
All_in$Treatment = factor(All_in$Treatment, levels = c("tumor","in"))

#t test
t_test <- t.test(ITH ~ Treatment, data = All_in, paired = TRUE, alternative = "greater")
p_value <-  t_test$p.value


colors = c("deepskyblue2","plum", "goldenrod3", "darkcyan","indianred1")

#Plot
ggplot(All_in, aes(x = Treatment, y = ITH , group = patient)) +
  geom_point(aes(colour = Treatment), position = position_dodge(width = 0)) +
  geom_line(aes(group = patient, colour = patient), linewidth = 2, alpha = 0.5, position = position_dodge(width = 0)) +
  scale_y_continuous(limits = c(-0.5,2.7), breaks = c(0,1,2,3))+
  scale_x_discrete(labels = c("in"="In vivo","tumor" = "Untreated"))+
  scale_colour_manual(values = colors)+
  ylab("genetic ITH") +
  #geom_text(aes(label = label), vjust = -0.5, position = position_dodge(width = 0.5), size = 6) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
         panel.spacing = unit(1, "lines"))  +
  # Add the p-value as annotation
  annotate("text", x = 1.5, y = 2.5, label = paste("p =", format.pval(p_value, digits = 3)), size = 8)

  ggsave("genetic_ITH_in_vivo_Entropy.pdf",path=plot_dir,width= 650,height=650,dpi=150,units="px")


#Save Figure data
write.table( All, file.path(Figures_data, "Figure_Sup_8A_Genetic_ITH.txt"), sep = "\t", quote = F)