library(CONIPHER)
library(tidyverse)
library(coin)


patient="OO99"
input_tsv_loc <- file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"data/OO99.All.mutect.ASCAT.header.conipher_input.tsv")
input_seg_tsv_loc <- file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"data/OO99.All.copynumber.caveman.seg.header.tsv")



#CLUSTERING

main_dir="/group/poetsch_projects/poetsch_sc/Conifer/OO99"

#default
out_dir=file.path(main_dir,"output")


conipher_clustering(case_id = "OO99",
             out_dir = out_dir,
             nProcs = 120,
             input_tsv_loc = input_tsv_loc,
             input_seg_tsv_loc = input_seg_tsv_loc)


#TREE BUILDING
tree_input_tsv_loc=file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/Clustering/",paste0(patient,".SCoutput.CLEAN.tsv"))

#default parameter
out_tree=file.path(out_dir,"tree")
dir.create(out_tree)


conipher_treebuilding(prefix = patient,
                      out_dir = out_tree,
                      input_tsv_loc = tree_input_tsv_loc,
                       min_cluster_size = 30)