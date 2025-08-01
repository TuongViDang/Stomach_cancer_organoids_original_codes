library(tidyverse)
library(ggpubr)

#Then run SigProfilerAssigment
#Import SigProfiler result
patient="OO99"
SigProfiler_dir=file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor",patient,"tree_30")
plot_dir=file.path(SigProfiler_dir,"plot")
dir.create(plot_dir)


SBS <- read.delim(file.path(SigProfiler_dir,"result_SBS/Assignment_Solution/Activities","Assignment_Solution_Activities.txt"))

keeps = which(colSums(SBS[,2:ncol(SBS)]) >0 )
SBS <- SBS[,c("Samples",names(keeps))]
SBS_long <- SBS %>% pivot_longer(cols=2:ncol(SBS),names_to = "Signature", values_to="Count")
SBS_long$Samples <- gsub("_specific","",SBS_long$Samples)
samples_order =  c( "Clone_1" ,"Clone_30", "Clone_8" ,"Clone_2", "Clone_4","Clone_12","Clone_15","Clone_13","Clone_23", "Clone_29", "Clone_5","Clone_14")
SBS_long$Samples <- factor(SBS_long$Samples, levels= samples_order )
sigs= rev(c( "SBS1",       "SBS5",     "SBS2",  "SBS13", "SBS3", "SBS4", "SBS6",  "SBS10a", "SBS16"    , "SBS17a",   "SBS17b",    "SBS18",   "SBS24" , "SBS29",   "SBS30",  "SBS31",  "SBS32",  "SBS35",  "SBS36",     "SBS39", "SBS40a",  "SBS41", "SBS51","SBS53",  "SBS54",  "SBS56",   "SBS57", "SBS58","SBS90", "SBS93"))
SBS_long$Signature <- factor(SBS_long$Signature, levels= sigs)

cols_SBS = rev(c("cadetblue3","cadetblue2","slateblue1","slateblue3","cyan2","yellow2","tomato","darkkhaki","mistyrose1","palegreen1","springgreen3","hotpink1","slategrey","darkolivegreen1", "bisque1","bisque2","bisque4","beige","azure","azure2","azure3","azure4","aliceblue",   rep("gray",5), "dimgray","lavenderblush1"))

SBS_long %>% ggplot(aes(x=Samples,y=Count,fill=Signature))+ geom_bar(stat="identity",position="stack",width=0.8) + scale_fill_manual(values=cols_SBS) + theme_classic() +
 ylab("Number of mutations in each signature") + guides(fill = guide_legend(reverse = TRUE))+
   theme(legend.title = element_blank(),
        axis.text.y.left = element_text(size=10),
        axis.text.x.bottom = element_text(size=10),
        axis.title.x.bottom = element_blank())
ggsave(file.path(plot_dir,"SBS_signature.pdf"), width = 14,height = 5 )

#Proportion
SBS_long %>% ggplot(aes(x= Samples,y=Count,fill=Signature))+ geom_bar(stat="identity",position="fill",width=0.6) + scale_fill_manual(values=cols_SBS) + theme_void() + 
 ylab("Number of mutations in each signature") + guides(fill = guide_legend(reverse = TRUE))+
   theme(legend.title = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x.bottom = element_text(size=8),
        axis.title.x.bottom = element_blank())   
ggsave(file.path(plot_dir,"SBS_signature_relative.pdf"), width = 14,height = 5)


# Calculate the proportion of each signature count per clone and total counts
SBS_long <- SBS_long %>%
  group_by(Samples) %>%
  mutate(
    Prop = Count / sum(Count),    # Proportion of each signature
    TotalCount = sum(Count)       # Total mutation count per clone
  ) %>%
  ungroup()


# Create the pie chart for each clone with size proportional to TotalCount
ggplot(SBS_long, aes(x = "", y = Prop, fill = Signature)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y") +  # Convert bar to pie chart
  facet_wrap(~ Samples, scales = "free", ncol = 3) +  # Facet to create one pie per clone
  scale_fill_manual(values = cols_SBS) +  # Custom SBS colors
  theme_void() +  # Remove axis and gridlines for a clean look
  labs(title = "Pie Charts of SBS Mutation Signatures per Clone") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 10)
  ) +
  #geom_text(aes(label = round(Prop * 100, 1)), position = position_stack(vjust = 0.5), size = 3) +  # Add labels for proportions
  scale_size_continuous(range = c(5, 30), trans = "sqrt") +  # Size based on total count
  guides(fill = guide_legend(reverse = TRUE))

ggsave(file.path(plot_dir,"SBS_signature_Pie_chart.pdf"), width = 10,height = 8 )


#Extract mutations belonging to artifact signature
artifact_SBS = c("SBS27", "SBS43","SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51","SBS52", "SBS53", "SBS54" , "SBS55" , "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", "SBS95")

#Manual test
Clone = "Clone_1"
SBS_long %>% filter( Samples ==  Clone , Count>0)
mut_contr = read.delim(file.path(SigProfiler_dir,"result_SBS/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities", paste0("Decomposed_Mutation_Probabilities_",Clone, "_specific.txt")))
range(mut_contr$SBS58 )
sum(mut_contr$SBS58 > 0.5)

#remove mutations has probability of SBS artifact >= 0.25
SBS_to_check = colnames(mut_contr)[5:90][colnames(mut_contr)[5:90] %in% artifact_SBS]

mutation_to_rm = c()
for (Clone in samples_order) {
mut_contr = read.delim(file.path(SigProfiler_dir,"result_SBS/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities", paste0("Decomposed_Mutation_Probabilities_",Clone, "_specific.txt")))
mut_contr <- mut_contr %>% mutate(mutation_id = paste0(":",paste0("chr",Chr),":",Pos,  gsub(">",":",   gsub(".\\[(.>.)\\].", ":\\1", MutationType)))) 
for (SBS in SBS_to_check){
    mut_rm_clone =  mut_contr$mutation_id[which(mut_contr[,SBS] >= 0.5)]
    n_mut = length(mut_rm_clone)
    mutation_to_rm  = c( mutation_to_rm ,mut_rm_clone)
    if (n_mut > 0){
    print(paste0(Clone,"_remove_",n_mut))}
}
}

#Filter vcf files
# Function to remove mutations from VCF and save the new VCF
OutputDir =  file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor",patient,"tree_30/Filtered/vcf_filtered")
dir.create(OutputDir, recursive = TRUE)
remove_mutations_from_vcf <- function(clone_name, mutation_to_rm) {
    # Read the original VCF file
    VcfDir = file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor",patient,"tree_30/vcf_SigProfiler")
    vcf_file_path <- file.path(VcfDir, paste0(clone_name, "_specific.mutation.vcf"))  # Path to the original VCF
    vcf_data <- read.table(vcf_file_path, header=FALSE, stringsAsFactors=FALSE)
    
    # Create a mutation ID for each row in the VCF (similar to mutation_id logic)
    vcf_data$mutation_id <- paste0(":", vcf_data$V1, ":", vcf_data$V2, ":", vcf_data$V4, ":", vcf_data$V5)
    
    # Filter out rows that contain mutations in the mutation_to_rm list
    mut_removed = intersect(vcf_data$mutation_id, mutation_to_rm)
    filtered_vcf <- vcf_data %>% filter(!(mutation_id %in% mut_removed))
    
    # Remove the mutation_id column
    filtered_vcf$mutation_id <- NULL
    n_row_removed = nrow(vcf_data) - nrow(filtered_vcf) 
    if(n_row_removed > 0 ){print(paste0("remove ",clone_name, ":",n_row_removed))}

    # Write the filtered data to a new VCF file
    output_vcf_path <- file.path(OutputDir, paste0(clone_name, "_filtered.vcf"))
    write.table(filtered_vcf, file=output_vcf_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    #collect mut remove 
    df = data.frame( mutation_id = mut_removed, clone = rep(clone_name, length(mut_removed)))
    df
}

# Loop through clones and process the VCF files
remove_muts = purrr::map_dfr(samples_order, ~remove_mutations_from_vcf(.x, mutation_to_rm = mutation_to_rm))

# Extract components from mutation_to_rm and create a data frame
mutation_data <- data.frame(
  CHR = sub("^:(.*?):.*", "\\1",  remove_muts$mutation_id),                # Extract Chromosome
  POS = sub("^.*:(\\d+):.*", "\\1", remove_muts$mutation_id),   # Extract Position
  CLONE = "removed" , # Assuming a single clone, change accordingly if needed
  REF = sub("^.*:(\\d+):([A-Z]):.*", "\\2", remove_muts$mutation_id),  # Extract Reference Allele
  ALT = sub("^.*:(\\d+):[A-Z]:(.*)", "\\2", remove_muts$mutation_id)  # Extract Alternate Allele
  
)

# Re-arrange columns
mutation_data <- mutation_data[, c("CHR", "POS", "CLONE", "REF", "ALT")]

# Write the data frame to a new file
write.table(mutation_data, file= paste0("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/",patient,"/tree_30/Filtered/vcf_filtered/mutation_to_rm_filtered.vcf"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
