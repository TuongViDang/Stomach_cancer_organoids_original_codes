library(tidyverse)
patient="OO77"
main_dir=file.path("/group/poetsch_projects/poetsch_sc/Conifer/",patient)
Clustering_output=file.path(main_dir,"output/Clustering/",paste0(patient,".SCoutput.CLEAN.tsv"))
file = read.delim(Clustering_output)
file = file %>% filter(CLUSTER %in% c(10,1,5,11,9,2,6,3))
write.table(file,file.path(main_dir,"output/Clustering/",paste0(patient,".SCoutput.CLEAN_final.tsv")), sep = "\t", quote = F )