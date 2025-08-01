library(tidyverse)

Patients = c( "OO77", "OO99", "OO100")

for (patient in Patients){

#read mutation

Conifer_dir = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient)
df = read.delim(file.path(Conifer_dir,"path_specific_muts_sep_truncal/All_tree_mutations.txt"))


 #Add Annotation (Annovar + VEP) : see Mutagen script
#Annovar 
Annovar <- read.csv(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2",patient,"Somatic/MuTect2/Annovar/",paste0(patient,".filtered.PASS.variants.annovar.hg38_multianno.csv")))
Annovar <- Annovar %>% mutate(mutation_id= paste0(":",sub("chr","",Annovar$Chr),":", Annovar$Start,":", Annovar$Ref, ":", Annovar$Alt),
                             Func_Annovar = paste0(Annovar$Func.refGene,"_",Annovar$ExonicFunc.refGene), 
                             Gene_Annovar = Gene.refGene) %>%                 
                        select(mutation_id,Func_Annovar,Gene_Annovar) 
#VEP 
VEP <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2",patient,"Somatic/MuTect2",paste0(patient,".filtered.PASS.vep_clean.vcf")), header = F )

colnames(VEP) <- c("Chr","Start","End","Ref","Alt","Func_VEP","Gene_VEP")
VEP <- VEP %>% mutate(mutation_id = paste0(":", sub("chr","",VEP$Chr),":", VEP$Start,":", VEP$Ref,":", VEP$Alt) ) %>% select(mutation_id,Func_VEP,Gene_VEP)

annotation <- left_join(Annovar,VEP, by = "mutation_id")


#COmbine  Annotation to df
df <- df %>% left_join(annotation, by = "mutation_id")

intergenic = "intergenic_variant"
splicing = c("splice_acceptor_variant" ,  "splice_donor_variant" , "splice_region_variant&intron_variant" , "splice_region_variant&intron_variant&NMD_transcript_variant"       
                , "splice_region_variant&intron_variant&non_coding_transcript_variant", "splice_region_variant&synonymous_variant" )
regulatory = c("3_prime_UTR_variant" , "3_prime_UTR_variant&NMD_transcript_variant"  , 
                "5_prime_UTR_variant"  ,  "5_prime_UTR_variant&NMD_transcript_variant" ,
                 "coding_sequence_variant&5_prime_UTR_variant", "downstream_gene_variant" ,
                 "regulatory_region_variant" ,"TF_binding_site_variant" , "upstream_gene_variant", splicing )
intron_nc = c("intron_variant","intron_variant&NMD_transcript_variant" , "intron_variant&non_coding_transcript_variant" ,"non_coding_transcript_exon_variant"     )

symnonymous = c( "synonymous_variant"  , "synonymous_variant&NMD_transcript_variant" )
pot_detrimental = c( "start_lost" ,  "stop_gained" ,   "stop_gained&splice_region_variant"     , "frameshift_variant"  ,"inframe_deletion"  , "missense_variant"    ,                                              
 "missense_variant&NMD_transcript_variant"  ,         "missense_variant&splice_region_variant"  )


#Driver mutations
#Mutagene
mutagene <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/mutagene",patient,"/maf_mutagen_result.txt"))
mutagene <- mutagene %>% mutate(Mutagene = ifelse(label %in% c("Driver","Potential driver"),TRUE,FALSE))


#Chasmplus
Chasmplus <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/Chasmplus/",patient,"result.clean.tsv"))
colnames(Chasmplus) <-  Chasmplus[1,]
Chasmplus <- Chasmplus[2:nrow(Chasmplus),]
Chasmplus <- Chasmplus %>% mutate(mutation_id = paste0(":", sub("chr","",Chrom), ":", Position, ":", Ref_Base, ":", Alt_Base),
                                  Chasmplus = ifelse(Score >= 0.1, TRUE, FALSE))

#Combine
Combined_MC <- left_join(mutagene,Chasmplus,by="mutation_id")

#Candra and Cscape use hg19, need to turn back to Hg38
#Candra
Candra_df <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/CanDrA/",patient,paste0(patient,".result.txt")))
Candra_df <- Candra_df[,c(1:4,12)]
Candra_df <- Candra_df %>% mutate(Candra = ifelse(CanDrA_Category == "Driver",TRUE,FALSE),
                                   mutation_19 = paste0(":",Chrom,":",Coordinate)) %>% select(mutation_19,Candra)  

#Cscape
Cscape_df <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/Cscape/",patient,paste0(patient,".lift_over_hg19_result.txt")))
Cscape_df <- Cscape_df %>% mutate(Cscape = ifelse(Coding.Score >0 & Warning == "oncogenic (high conf.)" ,TRUE,FALSE),
                                  mutation_19 = paste0(":",X..Chromosome,":",Position))    

Combined_hg19 <- merge(Cscape_df,Candra_df, by = "mutation_19",all=T) %>% mutate(Candra =ifelse(is.na(Candra),FALSE,Candra))

#reload lift over file
lift_overed_df <- read.delim(file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/dndscv",patient,paste0(patient,"_hg38to19_liftovered.txt")))
lift_overed_df <- lift_overed_df %>% mutate( mutation_19 = paste0(":",chr_hg19,":",  pos_hg19)) %>%  select(mutation_19,chr_hg38, pos_hg38)
lift_overed_df <- lift_overed_df[lift_overed_df$mutation_19 != ":NA:NA",]
lift_overed_df <- lift_overed_df [!duplicated(lift_overed_df$mutation_19),]

#translate to hg38
Combined_hg38 <- left_join(Combined_hg19,lift_overed_df,by="mutation_19") %>% 
                 select(chr_hg38, pos_hg38,Ref..Base, Mutant.Base,Cscape ,Candra) %>%
                  mutate(mutation_id = paste0(":",sub("chr","",chr_hg38),":", pos_hg38,":",Ref..Base,":", Mutant.Base))


#Combine with Mutagen and Chasmplus
Combined_all <- left_join(Combined_MC,Combined_hg38, by="mutation_id")
Combined_all <- Combined_all %>% mutate(Cscape = ifelse(is.na(Cscape),FALSE, Cscape), Candra = ifelse(is.na(Candra),FALSE,Candra))
Combined_all <- Combined_all %>% select(mutation_id,Gene, Protein_Change ,Driver_class,Mutagene,Chasmplus, Cscape ,Candra)



#Intogen
Intogen <- read.delim("/group/poetsch_projects/poetsch_sc/Driver_predict/Intogen/2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv")
length(unique(Intogen$SYMBOL)) #619 genes
Intogen_genes <- unique(Intogen$SYMBOL)

#Cosmic
Cosmic <- read.delim("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCITE/cancer_gene_lists/Cosmic_CancerGeneCensus_Tsv_v99_GRCh38/Cosmic_CancerGeneCensus_v99_GRCh38.tsv")
length(unique(Cosmic$GENE_SYMBOL))   #743 genes
Cosmic_genes <- unique(c(Cosmic$GENE_SYMBOL,Cosmic$SYNONYMS))

#OncoKB
OncoKB <- read.delim("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCITE/cancer_gene_lists/cancerGeneList.tsv")
length(unique(OncoKB$Hugo.Symbol))  #1138
OncoKB_genes <- unique(c(OncoKB$Hugo.Symbol,OncoKB$ Gene.Aliases))

#CancerMine
CancerMine <- read.delim("/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCITE/cancer_gene_lists/cancermine_collated.tsv")
length(unique(CancerMine$gene_normalized))  #6585
CancerMine_genes <- unique(CancerMine$gene_normalized)

#Annotate 
Combined_all <- Combined_all %>% mutate(Cosmic=ifelse(Gene %in% Cosmic_genes,TRUE,FALSE),
                                        Intogen=ifelse(Gene %in% Intogen_genes,TRUE,FALSE),
                                        OncoKB= ifelse(Gene %in% OncoKB_genes,TRUE,FALSE),
                                        CancerMine = ifelse(Gene %in% CancerMine_genes,TRUE,FALSE))

df_no_branch <- df %>% select(mutation_id,treeCLUSTER,Gene_Annovar, Func_Annovar) %>% unique()
Combined_all <- left_join(Combined_all,df_no_branch,by="mutation_id") 
Combined_all$n_mut_hits = rowSums(Combined_all[,5:8])
Combined_all$n_gene_hits =  rowSums(Combined_all[,9:12])
#export the table
write.table(Combined_all ,file.path(Conifer_dir  , paste0(patient,"_All_mutation_annotated_gene_mut_hit.txt")),row.names=T,sep="\t",quote=F)


#top drivers
top_driver <- Combined_all %>% filter(n_mut_hits > 0) %>% arrange(desc(n_gene_hits),desc(n_mut_hits)) %>% filter(n_gene_hits>0 & n_mut_hits >0 ) %>% filter(Gene != "")


write.table(top_driver ,file.path(Conifer_dir  , paste0(patient,"_top_drivers.txt")),row.names=T,sep="\t",quote=F)


#Reload
top_driver = read.delim(file.path(Conifer_dir , paste0(patient,"_top_drivers.txt")))
#Gastric driver genes Totoki et al
gastric_totoki = read.delim("/group/poetsch_projects/poetsch_sc/Driver_predict/Gastric_drivers_Totoki.txt")
gastric_genes = gastric_totoki$HUGO.Symbol

#Top driver in gastric genes
top_driver = top_driver %>% mutate(Totoki = ifelse(Gene %in% gastric_genes, TRUE,FALSE)) 



#Missense mutation
Missense <- df  %>% filter(genic == "potential detrimental") %>% mutate(patient = patient)

detrimental = c( "start_lost" ,  "stop_gained" ,   "stop_gained&splice_region_variant"     , "frameshift_variant"  ,"inframe_deletion" )

write.table(Missense,file.path(Conifer_dir , paste0(patient,"_missense.txt")),row.names=T,sep="\t",quote=F)

Det_gastric_missense = Missense %>% filter( (Gene_Annovar %in% gastric_genes) & (Func_VEP %in% detrimental))
Det_gastric_missense_df <- df %>% filter(mutation_id %in% Det_gastric_missense$mutation_id )

###DRIVERS: in top_driver + Det_gastric_missense
top_drivers_mut = top_driver %>% filter(Totoki == T) %>% pull(mutation_id)
Final_drivers_mut = unique(c( top_drivers_mut ,Det_gastric_missense$mutation_id ))
Final_drivers = df %>% filter(mutation_id %in%  Final_drivers_mut )

write.table(Final_drivers,file.path(Conifer_dir , paste0(patient,"_Final_drivers.txt")),row.names=T,sep="\t",quote=F)







#Which cluster the drivers belong to
driver_cluster = top_driver %>% group_by(treeCLUSTER) %>% summarize(count=n()) %>% arrange(desc(count))
driver_cluster = driver_cluster %>% mutate(treeCLUSTER = ifelse(is.na(treeCLUSTER),"not in tree",treeCLUSTER))
driver_cluster$treeCLUSTER = factor(driver_cluster$treeCLUSTER,levels= c(setdiff(driver_cluster$treeCLUSTER,"not in tree"),"not in tree"))

driver_cluster %>% ggplot(aes(x=treeCLUSTER,y=count,fill=treeCLUSTER))+
                  geom_bar(stat="identity",width=0.7) + 
                   scale_fill_brewer(palette = "Set3", name = "Cluster")  + theme_minimal() + ylab("Number of driver mutation") +
                   theme(legend.position = "None",
                        axis.text.x.bottom=element_text(size=15),
                        axis.title.x.bottom=element_text(size=15),
                         axis.text.y.left=element_text(size=15),
                         axis.title.y.left=element_text(size=15)
                   )

ggsave(file.path(plot_dir,"Driver_cluster_distribution.pdf"),width=13,height=7) 

#plot cluster driver enrichment
clusters = intersect(names(table(df$treeCLUSTER)),names(table(top_driver$treeCLUSTER )))

n_total_cluster = c(table(df$treeCLUSTER))[names(table(df$treeCLUSTER)) %in% clusters ]

n_driver_cluster = c(table(top_driver$treeCLUSTER))[names(table(top_driver$treeCLUSTER)) %in% clusters ]

df_driver_cluster <- data.frame(cluster = clusters, n_driver_cluster = n_driver_cluster, n_total_cluster= n_total_cluster)
df_driver_cluster <- df_driver_cluster %>% mutate(driver_enrichment= 10000*n_driver_cluster/n_total_cluster )


df_driver_cluster %>% ggplot(aes(x= cluster,y= driver_enrichment,fill= cluster)) + geom_bar(stat = 'identity', width = 0.6) + theme_minimal() +
ylab("Driver Enrichment Score") +
theme(legend.position = "None",
                        axis.text.x.bottom=element_text(size=15),
                        axis.title.x.bottom=element_blank(),
                         axis.text.y.left=element_text(size=15),
                         axis.title.y.left=element_text(size=15))

ggsave(file.path(plot_dir,"Driver_cluster_enrichment.pdf"),width=9,height=5) 






#Which branchs drivers belong to

df_driver_branch <- df %>% filter(mutation_id %in% top_driver$mutation_id )

branch = c("truncal",paste0("branch",1:6))
n_drivers =  c(table(df_driver_branch$branches),"branch6"=0)[branch]
names(n_drivers)= NULL
n_total =  c(table(df$branches))[branch]
names(n_total)= NULL

df_driver_summary <- data.frame(branch = branch,
                                n_drivers = n_drivers,
                                n_total = n_total )

df_driver_summary <- df_driver_summary %>% mutate( driver_enrichment = 10000 * n_drivers/n_total)
df_driver_summary$branch <- factor(df_driver_summary$branch, levels = branch)

df_driver_summary %>% ggplot(aes(x= branch,y= driver_enrichment,fill=branch)) + geom_bar(stat = 'identity', width = 0.6) + theme_minimal() +
ylab("Driver Enrichment Score") +
theme(legend.position = "None",
                        axis.text.x.bottom=element_text(size=15),
                        axis.title.x.bottom=element_blank(),
                         axis.text.y.left=element_text(size=15),
                         axis.title.y.left=element_text(size=15))

ggsave(file.path(plot_dir,"Driver_path_enrichment.pdf"),width=9,height=5) 


#Missense mutation

Missense <- df  %>% filter(genic == "potential detrimental") %>% mutate(patient = patient)

write.table(Missense,file.path(Conifer_dir , paste0(patient,"_missense.txt")),row.names=T,sep="\t",quote=F)
}