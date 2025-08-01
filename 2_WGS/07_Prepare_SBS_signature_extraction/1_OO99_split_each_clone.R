library(tidyverse)
library(GenomicRanges)

#specify paths to files from Conipher tree
patient="OO99"
Main_dir = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient)
clone_specific_muts_dir=file.path(Main_dir,"tree30","each_clone","mutation_files")
dir.create(clone_specific_muts_dir,recursive=T)


#Load the tables
mut_file = file.path(Main_dir,"output/tree_30/Trees/treeTable.tsv")
cluster_file = file.path(Main_dir,"output/tree_30/Trees/clusterInfo.txt")
tree_file= file.path(Main_dir,"output/tree_30/Trees/allTrees.txt")


#only keep clusters in tree

cluster_df = read.delim(cluster_file)
cluster_df = cluster_df %>% filter(treeClust == "TRUE")

mut_df = read.delim(mut_file)
nrow(mut_df)
mut_df = mut_df %>% filter(treeCLUSTER %in% cluster_df$clusterID )
nrow(mut_df)
clone_list = names(table(mut_df$treeCLUSTER))

mut_df_distict = mut_df %>% distinct(mutation_id, .keep_all = T)
table(mut_df_distict$treeCLUSTER)
#Extract mutations belonging to each clone

clone_specific_muts = function(clone, out_dir){
muts = mut_df %>% filter(treeCLUSTER == clone) 
df <- unique(muts[, c("CHR", "POS", "REF", "ALT")])
df$CHR <- paste0("chr",df$CHR)
gr <- GRanges(seqnames = df$CHR,
              ranges = IRanges(start = df$POS, end = df$POS),
              strand = "*",
              REF = df$REF,
              ALT = df$ALT)
colnames(df) <- NULL
write.table(df,file.path(out_dir,paste0("Mutations_Clone_", clone ,".txt")),row.names=FALSE,quote=F,sep="\t")
gr
}
Granges_list = purrr::map(clone_list,~clone_specific_muts(.x, out_dir = clone_specific_muts_dir))
names(Granges_list) <- clone_list

#####RUN:Extract subsets of vcf using Subset_clone_specific_vcf.sh
#####RUN: SigProfiler with the vcf files




