library(tidyverse)
library(ggpubr)


#Import SigProfiler result
patient="OO100"
SigProfiler_dir=file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor",patient,"tree_30","Filtered")
plot_dir=file.path(SigProfiler_dir,"plot")
dir.create(plot_dir)


#Figures data path
main_path='/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis'
Figures_data <- file.path(main_path,'txt_files', "three_patients", "Figures_data")
dir.create(Figures_data)




SBS <- read.delim(file.path(SigProfiler_dir,"result_Extractor/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities","COSMIC_SBS96_Activities.txt"))
SBS <- SBS %>% filter(Samples != "mutation_to_rm_filtered")

keeps = which(colSums(SBS[,2:ncol(SBS)]) >0 )
SBS <- SBS[,c("Samples",names(keeps))]
SBS_long <- SBS %>% pivot_longer(cols=2:ncol(SBS),names_to = "Signature", values_to="Count")
SBS_long$Samples <- gsub("_filtered","",SBS_long$Samples)
samples_order = samples_order = c( "Clone_1" ,"Clone_14", "Clone_5" ,"Clone_7", "Clone_23","Clone_4","Clone_18","Clone_2", "Clone_19","Clone_22","Clone_6","Clone_11")
SBS_long$Samples <- factor(SBS_long$Samples, levels= samples_order )
table(SBS_long$Signature)

sigs= rev(c( "SBS1",       "SBS5",   "SBS18"))
SBS_long$Signature <- factor(SBS_long$Signature, levels= sigs)

cols_SBS = rev(c("cadetblue3","cadetblue2","hotpink1"))

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

#Save data
#Change name of clone
# Define mapping for replacements
clusterID_map <- c(`14` = 2, `5` = 3, `7` = 4, `23` = 5, `4` = 6, `18` = 7, 
                   `2` = 8, `19` = 9, `22` = 10, `6` = 11, `11` = 12)
# Apply transformations
SBS_long <- SBS_long %>%
  mutate(Samples = sub("Clone_","", Samples),
         Clone = ifelse( Samples %in% names(clusterID_map), clusterID_map[as.character(Samples)], Samples)
  ) %>%
   select(Clone,Signature,Count,TotalCount, Prop ) %>% 
   rename(Proportion = Prop) %>% arrange(Clone)


write.table(SBS_long, file.path(Figures_data, paste0("Figure_Sup_3B_SBS_signature_per_clone_", patient, ".txt")), sep = "\t", quote = F)