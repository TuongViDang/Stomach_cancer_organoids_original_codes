library(tidyverse)
patient="OO100"
main_dir=file.path("/group/poetsch_projects/poetsch_sc/Conifer/",patient)
Clustering_output=file.path(main_dir,"output/Clustering/",paste0(patient,".SCoutput.CLEAN.tsv"))
file = read.delim(Clustering_output)
file = file %>% mutate( CLUSTER =  case_when(CLUSTER == 20 ~ 23,
                                              CLUSTER == 12 ~ 19,
                                              TRUE ~ CLUSTER))
file = file %>% filter(CLUSTER %in% c(1,14,5,7,23,4,18,2,19,22,6,11))
write.table(file,file.path(main_dir,"output/Clustering/",paste0(patient,".SCoutput.CLEAN_merged.tsv")), sep = "\t", quote = F )