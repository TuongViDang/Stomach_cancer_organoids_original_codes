library(tidyverse)
library(pheatmap)

patient="OO99"


main_dir = file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA", patient, "tree30_change_name")
plot_dir = file.path(main_dir,"plots","map_no")
dir.create(plot_dir, recursive = T)
rds_dir = file.path(main_dir,"rds")
dir.create(rds_dir, recursive = T)
txt_files = file.path(main_dir,"txt_files")
dir.create(txt_files , recursive = T)
#Add Annotation (Annovar + VEP) : see Mutagen script
#Annovar 
Annovar <- read.csv(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2",patient,"/Somatic/MuTect2/Annovar/",paste0(patient,".filtered.PASS.variants.annovar.hg38_multianno.csv")))
Annovar <- Annovar %>% mutate(mutation_id= paste0(":",sub("chr","",Annovar$Chr),":", Annovar$Start,":", Annovar$Ref, ":", Annovar$Alt),
                             Func_Annovar = paste0(Annovar$Func.refGene,"_",Annovar$ExonicFunc.refGene), 
                             Gene_Annovar = Gene.refGene) %>%                 
                        select(mutation_id,Func_Annovar,Gene_Annovar) 
#VEP 
VEP <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/",patient,"Somatic/MuTect2",paste0(patient,".filtered.PASS.vep_clean.vcf")), header = F )

colnames(VEP) <- c("Chr","Start","End","Ref","Alt","Func_VEP","Gene_VEP")
VEP <- VEP %>% mutate(mutation_id = paste0(":", sub("chr","",VEP$Chr),":", VEP$Start,":", VEP$Ref,":", VEP$Alt) ) %>% select(mutation_id,Func_VEP,Gene_VEP)

annotation <- left_join(Annovar,VEP, by = "mutation_id")



#check gene name lack and add VEP if necessary
sum(annotation$Gene_Annovar == "."  )
sum(annotation$Gene_VEP == "" )

annotation <- annotation %>%mutate(Gene_name=ifelse(Gene_Annovar == "." & Gene_VEP != "", Gene_VEP, Gene_Annovar ))

#check
sum(annotation$Gene_name == "." )    #fewer than Gene_Annovar 

#Load the file
Combined_all = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Integrate_WGS_scRNA", patient, "tree30","txt_files",paste0(patient,"_Combined_all.txt")))
table(Combined_all$treeCLUSTER)
#Change the name of clone

 Combined_all <- Combined_all %>% mutate(treeCLUSTER_new = case_when(treeCLUSTER == 30 ~ 2,
                                                                     treeCLUSTER == 8 ~ 3,
                                                                     treeCLUSTER == 2 ~ 4,
                                                                     treeCLUSTER == 4 ~ 5,
                                                                     treeCLUSTER == 12 ~ 6,
                                                                     treeCLUSTER == 15 ~ 7,
                                                                     treeCLUSTER == 13 ~ 8,
                                                                     treeCLUSTER == 23 ~ 9,
                                                                     treeCLUSTER == 29 ~ 10,
                                                                     treeCLUSTER == 5 ~ 11,
                                                                     treeCLUSTER == 14 ~ 12,
                                                                     TRUE ~ treeCLUSTER))
Combined_all$treeCLUSTER <- Combined_all$treeCLUSTER_new
Combined_all <- Combined_all %>% select(-treeCLUSTER_new)
table(Combined_all$treeCLUSTER)

#save this file
write.table(Combined_all, file.path(txt_files,paste0(patient,"_Combined_all.txt")),sep = "\t", quote = F)


#Load the file
Combined_all = read.delim(file.path(txt_files,paste0(patient,"_Combined_all.txt")))
table(Combined_all$treeCLUSTER)

###scRNA-seq expression ###
# Load libraries
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(RCurl)
library(cowplot)
library(patchwork)

#Specify the path
rds_all_samples_path = "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/rds"
seurat_integrated <- readRDS(file.path(rds_all_samples_path ,"clustered_0.8integrated_doublet_rm_filtered_all_patients_qc_cca.rds"))




#change the name of sample
meta <- seurat_integrated@meta.data
meta <- meta %>% mutate(patient = lapply(orig.ident,function(x){strsplit(x,"-")[[1]][1]}))
meta <- meta  %>% select(orig.ident,patient, Phase, seurat_clusters)
seurat_integrated@meta.data <- meta





#Plot any mutation
#OO99-no
treatment = "no"
organoid = paste0(patient,"-",treatment)
plot_dir_sub = file.path(plot_dir,organoid)
dir.create(plot_dir_sub)

seurat_sub = subset(seurat_integrated, subset = orig.ident == organoid )

#Before integration
DefaultAssay(seurat_sub)<-'RNA'
options(future.globals.maxSize = 4294967296)  # 4 GB

seurat_sub <- SCTransform(seurat_sub, vst.flavor = "v1",
                                 #method = "glmGamPoi",
                                 verbose = FALSE)

seurat_sub <- RunPCA(seurat_sub, npcs = 30, verbose = F)
seurat_sub <- RunUMAP(seurat_sub, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(seurat_sub, reduction = "umap", group.by = "orig.ident") 
ggsave(paste0("DimPlot_",organoid,".pdf"),path = plot_dir_sub,width = 10000,height = 8000,dpi = 1500,units = 'px')


#Phase
DimPlot(seurat_sub, reduction = "umap", group.by = "Phase") 
ggsave(paste0("DimPlot_",organoid,"_Phase.pdf"),path = plot_dir_sub,width = 10000,height = 8000,dpi = 1500,units = 'px')
#Cells belong to which clones
#skip clone 1 and 5
Combined_organoid = Combined_all %>% filter(sample == organoid)
Combined_organoid_CB = Combined_organoid %>% filter(mutated == TRUE) %>% group_by(CB_seurat) %>% summarize(clone = treeCLUSTER)

mutated_CB = as.data.frame(table(Combined_organoid_CB$CB_seurat,Combined_organoid_CB$clone)) 
colnames(mutated_CB) = c("CB_seurat","Clone","Freq")
mutated_CB = mutated_CB %>% pivot_wider(names_from = "Clone", values_from = "Freq")
head(mutated_CB)


clone_columns = c( "2", "3","4","6", "7" ,"8","9","10","11","12")
mutated_CB_clone = mutated_CB[,intersect(clone_columns,colnames(mutated_CB)) ]

Find_max_col = function(row) {
  max_value = max(row) # store the maximum value
  if (max_value > 0) {
    return(colnames(mutated_CB_clone)[which(row == max_value)][1])
  } else {
    return("UnI")
  }
}


mutated_CB_clone$Clone = apply(mutated_CB_clone, 1, Find_max_col)
mutated_CB$Clone =  mutated_CB_clone$Clone

table(mutated_CB$Clone)
#For clones of the same branch, priotize children over parents
mutated_CB <- mutated_CB %>% 
  mutate(Clone_pri = case_when(
    Clone %in% c("2","3") & `4` > 0 ~ "4",
    Clone == "2" & `3` > 0 ~ "3",
    Clone %in% c("6","7") & `8` > 0 ~ "8",
    Clone %in% c("6","7") & `9` > 0 ~ "9",
    Clone %in% c("6","7") & `10` > 0 ~ "10",
    Clone == "6" & `7` > 0 ~ "7",
    TRUE ~ Clone
  ))


mutated_CB$Clone = mutated_CB$Clone_pri
table(mutated_CB$Clone)
#save to seurat


meta = seurat_sub@meta.data
meta$CB_seurat = rownames(meta) 
meta = meta %>% left_join(mutated_CB[,c("CB_seurat","Clone")], by = "CB_seurat")
seurat_sub@meta.data = meta
rownames(seurat_sub@meta.data) = seurat_sub@meta.data$CB_seurat

seurat_sub$Clone = factor(seurat_sub$Clone , levels = clone_columns  )


seurat_sub_clone <- subset(seurat_sub, Clone %in% clone_columns )

colors = c( `2`= "aquamarine2", `3` = "moccasin" , `4` = "dodgerblue4", `6` ="yellow2",
            `7` = "lightseagreen", `8` = "lightsalmon4" , `9`= "plum" , 
             `10` = "tomato", `11` = "mediumseagreen", `12` = "purple")


DimPlot(seurat_sub_clone, reduction = "umap", group.by = "Clone" , cols = colors)  + labs(title = "Subclone")
ggsave(paste0("DimPlot_",organoid,"_Clone.pdf"),path = plot_dir_sub,width = 10000,height = 8000,dpi = 1500,units = 'px')                          

#save this seurat
saveRDS(seurat_sub, file.path(rds_dir,paste0(organoid,"_Clones.RDS")))

#load the file
seurat_sub = readRDS(file.path(rds_dir,paste0(organoid,"_Clones.RDS")))