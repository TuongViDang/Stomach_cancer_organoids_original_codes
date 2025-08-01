
library(CONIPHER)
library(tidyverse)
library(coin)



patient="OO100"

main_dir=file.path("/group/poetsch_projects/poetsch_sc/Conifer/",patient)
tree_input_tsv_loc=file.path(main_dir,"output/Clustering/",paste0(patient,".SCoutput.CLEAN_merged.tsv"))

snv = read.delim(tree_input_tsv_loc)
total_n_muts = length(unique(snv$mutation_id))

out_dir=file.path(main_dir,"output")
#min = 30
out_tree=file.path(out_dir,"tree_30_final")
dir.create(out_tree)


conipher_treebuilding(prefix = patient,
                      out_dir = out_tree,
                      input_tsv_loc = tree_input_tsv_loc,
                      min_ccf = 0.15,
                      min_cluster_size = 30 )
