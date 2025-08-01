
library(tidyverse)
patient = "OO100"
Main_dir = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient)
mut_file = read.delim(file.path(Main_dir,"output/tree_30_final/Trees/treeTable.tsv"))
table(mut_file$treeCLUSTER)

#filter cluster in tree 30
cluster_file = read.delim(file.path(Main_dir,"output/tree_30_final/Trees/clusterInfo.txt"))
tree_clusters = cluster_file  %>% filter(treeClust == TRUE) %>% pull(clusterID) %>% unique()
mut_file_tree = mut_file %>% filter(treeCLUSTER %in% tree_clusters) 
table(mut_file_tree$treeCLUSTER)
length(unique(mut_file_tree$mutation_id))

#filter mutation with artifact identified by Sigprofiler
mut_to_rm = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/",patient,"/tree_30/Filtered/vcf_filtered/mutation_to_rm_filtered.vcf"), header = F)
mut_to_rm = mut_to_rm[,c(1,2,4,5)]
colnames(mut_to_rm) = c("CHR","POS","REF","ALT")
mut_to_rm = mut_to_rm %>% mutate(mutation_id = paste0(":",sub("chr","",CHR),":",POS,":",REF,":",ALT))
mut_to_rm_list = mut_to_rm %>% pull(mutation_id)

#Check if they are detrimental
detrimental = read.delim(file.path(Main_dir,paste0(patient,"_detrimental_Sanika.txt")))
detrimental = detrimental %>% distinct(mutation_id, .keep_all=T) %>% select( mutation_id, Func_Annovar, Gene_Annovar ,Func_VEP, Gene_VEP)
detrimental_mut_id = detrimental %>% pull(mutation_id) %>% unique()
sum(mut_to_rm$mutation_id %in% detrimental_mut_id)
mut_retain = mut_to_rm_list[which(mut_to_rm$mutation_id %in% detrimental_mut_id)]
detrimental %>% filter(mutation_id %in% mut_retain)
mut_to_rm_list_final = setdiff(mut_to_rm_list, mut_retain)

mut_file_tree_filtered = mut_file_tree %>% filter(!(mutation_id %in% mut_to_rm_list_final ))
length(unique(mut_file_tree_filtered$mutation_id))

#save this file
write.table(mut_file_tree_filtered,file.path(Main_dir,"output/tree_30_final/Trees/treeTable_filtered_artifact_Sigprofiler.tsv"), sep = "\t", quote = FALSE )