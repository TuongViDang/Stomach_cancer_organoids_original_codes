# Load libraries
library(Seurat)
library(tidyverse)
library(liana)
library(scater)
library(magrittr)


#Specify the path
main_path <-'/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
rds_3_path  <-file.path(main_path,'rds', "three_patients")
txt_files_3 <-file.path(main_path,'txt_files', "three_patients")
plot_path <-file.path(main_path,'plots_3_patients')
plot_result_path <- file.path(plot_path,"Liana")
dir.create(plot_result_path)

#Load input file
#Specify the input argument
seurat_integrated <- readRDS(file.path(rds_3_path,'clustered_0.35_integrated_seurat_ref.rds'))

context_df_dict<-list()
for (sample_name in sort(unique(seurat_integrated$orig.ident))){
    sdata = subset(x = seurat_integrated, subset = orig.ident == sample_name)
    DefaultAssay(sdata)="RNA"
    options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10GB
    sdata<- NormalizeData(sdata,normalization.method = "LogNormalize")
    sdata.sce <- Seurat::as.SingleCellExperiment(sdata)
    liana_res <- liana_wrap(sce = sdata.sce,
                        idents_col='seurat_clusters',
                        assay='RNA', # specify the Seurat assay
                        assay.type = 'logcounts',
                        expr_prop=0.1,
                        verbose=T,
                        parallelize = T, workers = 50)

    liana_aggregate.magnitude <- liana_aggregate(liana_res = liana_res )
                                  #, aggregate_how='magnitude', verbose = F)

    # retain only the aggregate magnitude rank score
    # and format for input to liana_tensor_c2c function
    liana_aggregate.magnitude<-liana_aggregate.magnitude
    #colnames(liana_aggregate.magnitude)<-c('source', 'target', 'ligand.complex', 'receptor.complex', 'magnitude_rank')
    liana_aggregate.magnitude[['orig.ident']]<-sample_name

    context_df_dict[[sample_name]]<-liana_aggregate.magnitude
}

saveRDS(context_df_dict,file.path(rds_3_path,"Liana_context_df_dict.rds"))


