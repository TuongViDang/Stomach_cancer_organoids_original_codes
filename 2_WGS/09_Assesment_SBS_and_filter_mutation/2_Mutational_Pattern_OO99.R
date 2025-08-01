library(tidyverse)
library(MutationalPatterns)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)


#specify paths to files from Conipher tree
patient="OO99"
Main_dir = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient)
clone_specific_muts_dir=file.path(Main_dir,"tree30","each_clone")
dir.create(clone_specific_muts_dir,recursive=T)
plot_path= file.path(clone_specific_muts_dir, "plot")
dir.create(plot_path)

#Load the tables
mut_file = file.path(Main_dir,"output/tree_30/Trees/treeTable.tsv")
cluster_file = file.path(Main_dir,"output/tree_30/Trees/clusterInfo.txt")
tree_file= file.path(Main_dir,"output/tree_30/Trees/allTrees.txt")

#Run split_each_clone.R
#Now use MutationalPatterns

vcf_files <- list.files(file.path(clone_specific_muts_dir,"vcf"),full.names=TRUE)

vcf_names <- sub("_specific.mutation.vcf","",basename(vcf_files))
names(vcf_files) <- vcf_names
vcf_names_ordered <- c( "Clone_1" ,"Clone_30", "Clone_8" ,"Clone_2", "Clone_4","Clone_12","Clone_15","Clone_13","Clone_23", "Clone_29", "Clone_5","Clone_14")

# Match and reorder the vcf_files according to vcf_names_ordered
vcf_files_ordered <- vcf_files[match(vcf_names_ordered, vcf_names)]

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
grl <- read_vcfs_as_granges(vcf_files_ordered, 
                            vcf_names_ordered, 
                            ref_genome, 
                            type="snv", 
                            predefined_dbs_mbs=TRUE)



#12 types
muts <- mutations_from_vcf(grl[[1]])
head(muts, 12)

#change to 6 types
types <- mut_type(grl[[1]])
head(types, 12)

#Context
context <- mut_context(grl[[1]], ref_genome)
head(context, 12)
type_context <- type_context(grl[[1]], ref_genome)
lapply(type_context, head, 12)


#count mutation type occurrences for all VCF objects
type_occurrences <- mut_type_occurrences(grl, ref_genome)
type_occurrences

png(file.path(plot_path,"mutation_spectrum_all.png"))
plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE)
dev.off()

#each clone
#type_occurrences$total = rowSums(type_occurrences)

pdf(file.path(plot_path,"mutation_spectrum_each_clone.pdf"))
plot_spectrum(type_occurrences, CT = TRUE,  by = vcf_names_ordered, error_bars= 'none', indv_points = FALSE)
dev.off()

#Mut matrix
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)

pdf(file.path(plot_path,"96_profile_each_clone.pdf"),width = 15, height = 12)
plot_96_profile(mut_mat)
dev.off()

###After Filtered
#removed mutation

vcf_file <- "/group/poetsch_projects/poetsch_sc/Conifer/OO99/tree30/each_clone/Filtered/vcf_filtered/mutation_to_rm_filtered.vcf"
ref_genome="BSgenome.Hsapiens.UCSC.hg38"
grl <- read_vcfs_as_granges(vcf_file,
                            "removed", 
                            ref_genome, 
                            type="snv", 
                            predefined_dbs_mbs=TRUE)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

png(file.path(plot_path,"mutation_spectrum_removed.png"))
plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE)
dev.off()

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)

pdf(file.path(plot_path,"96_profile_removed.pdf"),width = 10, height = 5)
plot_96_profile(mut_mat)
dev.off()


#Filtered 
filtered_dir="/group/poetsch_projects/poetsch_sc/Conifer/OO99/tree30/each_clone/Filtered/vcf_filtered"
vcf_files <- list.files(filtered_dir,full.names=TRUE)
vcf_files <- vcf_files[1:12]
vcf_names <- sub("_filtered.vcf","",basename(vcf_files))
names(vcf_files) <- vcf_names
vcf_names_ordered <- c( "Clone_1" ,"Clone_30", "Clone_8" ,"Clone_2", "Clone_4","Clone_12","Clone_15","Clone_13","Clone_23", "Clone_29", "Clone_5","Clone_14")
vcf_files_ordered <- vcf_files[match(vcf_names_ordered, vcf_names)]

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
grl <- read_vcfs_as_granges(vcf_files_ordered, 
                            vcf_names_ordered, 
                            ref_genome, 
                            type="snv", 
                            predefined_dbs_mbs=TRUE)
type_occurrences <- mut_type_occurrences(grl, ref_genome)

png(file.path(plot_path,"mutation_spectrum_all_filtered.png"))
plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE)
dev.off()

#each clone
#type_occurrences$total = rowSums(type_occurrences)

pdf(file.path(plot_path,"mutation_spectrum_each_clone_filtered.pdf"))
plot_spectrum(type_occurrences, CT = TRUE,  by = vcf_names_ordered, error_bars= 'none', indv_points = FALSE)
dev.off()

#Mut matrix
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
head(mut_mat)

pdf(file.path(plot_path,"96_profile_each_clone_filtered.pdf"),width = 15, height = 12)
plot_96_profile(mut_mat)
dev.off()