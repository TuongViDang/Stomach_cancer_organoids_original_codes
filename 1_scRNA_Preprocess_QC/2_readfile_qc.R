library(dplyr)
library(Seurat)
library(patchwork)
library(purrr)

# Specify the path
data_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/links_gastric_cancer_data'
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
plot_path=file.path(main_path,'plots')
plot_path_before_qc<-file.path(plot_path,'qc','before_qc')
plot_path_after_qc<-file.path(plot_path,'qc','after_qc')
rds_path<-file.path(main_path,'rds')
functions_path <- '/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/downstream_analysis'

source(file.path(functions_path,'functions.R'))


#READ DATA
patients<-list.files(data_path,full.names = F,pattern = '*')
patients

all_patients<-read_scRNAseq(data_path,patients = patients)  #function file
all_patients

#Merge all patiens

all_patients_qc <- merge(all_patients[[1]], c(all_patients[[2]],all_patients[[3]],all_patients[[4]],all_patients[[5]],
                                              all_patients[[6]],all_patients[[7]],all_patients[[8]],all_patients[[9]],
                                             all_patients[[10]],all_patients[[11]],all_patients[[12]],all_patients[[13]],
                                             all_patients[[14]],all_patients[[15]],all_patients[[16]],all_patients[[17]]),
                         add.cell.ids = patients)
                      
head(all_patients_qc@meta.data)
tail(all_patients_qc@meta.data)

#Add QC metrics

all_patients_qc<-qc_scRNAseq(all_patients_qc)   #function file

#Violin plots

features<-colnames(all_patients_qc@meta.data)[c(2,3,5,6)]
y_axis_name=c('UMI count','number of features','mito genes ratio','ribo genes ratio')


plot_violin_scRNAseq(object = all_patients_qc,feature = features,y_axis_name = y_axis_name
                    ,plot_path = plot_path_before_qc,plot_name='violin_qc_before_filter.pdf',width = 15000, height = 15000, dpi = 1500, units = 'px',device='pdf'
)      #function file


#Plot number of cell in each sample
metadata<-all_patients_qc@meta.data
colnames(metadata)

plot_number_cells(metadata,plot_name='number_of_cells_bf_qc.pdf',plot_path=plot_path_before_qc,width = 15000, height = 10000, units = 'px', dpi = 1500)  #function file


#Density Plot (perform on metadata of merged object

umi_cut <- 3000
mt_cut <- 0.15
gene_cut <- 500
log10_cut <- 0.8

Plot_density_qc(metadata,umi_cut=umi_cut,gene_cut=gene_cut, mt_cut=mt_cut, log10_cut=log10_cut,
                plot_name='Plot_density_before_qc.pdf' , plot_path=plot_path_before_qc, width=15000 , height=10000 , dpi=1500, units='px')  #function file

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
cor_umi_feature_mito<-metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point(size=0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500,linetype='dashed') +
  geom_hline(yintercept = 250,linetype='dashed') +
  facet_wrap(vars(orig.ident))

ggsave(filename='Correlation_umi_feature_mito_bf_qc.pdf',path = plot_path_before_qc,
       plot=cor_umi_feature_mito,width =3000,height = 3000,units = 'px')


#CELL FILTERING
filtered_all_patients_qc <- subset(x = all_patients_qc,
                                   subset = (nCount_RNA >= umi_cut) &
                                     (nFeature_RNA >= gene_cut) &
                                     (log10GenesPerUMI >= log10_cut) &
                                     (mitoRatio <= mt_cut)
                                   # & (riboRatio < 0.4)
)

#Produce those above plots for data after filtering
plot_violin_scRNAseq(object = filtered_all_patients_qc,feature = features,y_axis_name = y_axis_name
                     ,plot_path = plot_path_after_qc,plot_name='violin_qc_after_qc.pdf',width = 15000, height = 15000, dpi = 1500, units = 'px',device = 'pdf'
)

metadata_clean<- filtered_all_patients_qc@meta.data
dim(metadata_clean)

plot_number_cells(metadata_clean,plot_name='number_of_cells_after_qc.pdf',plot_path=plot_path_after_qc,width = 15000, height = 10000, units = 'px', dpi = 1500)
Plot_density_qc(metadata_clean,umi_cut,gene_cut,log10_cut,mt_cut,
                plot_name='Plot_density_after_qc.pdf' , plot_path=plot_path_after_qc, width=15000 , height=10000 , dpi=1500, units='px')


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
cor_umi_feature_mito<-metadata_clean %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point(size=0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500,linetype='dashed') +
  geom_hline(yintercept = 250,linetype='dashed') +
  facet_wrap(vars(orig.ident))

ggsave(filename='Correlation_umi_feature_mito_aft_qc.pdf',path = plot_path_after_qc,
       plot=cor_umi_feature_mito,width =3000,height = 3000,units = 'px')

# Create .RData object to load at any time
saveRDS(filtered_all_patients_qc, file=file.path(rds_path,'filtered_all_patients_qc.rds'))






