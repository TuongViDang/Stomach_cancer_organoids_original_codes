library(tidyverse)
library(GenomicRanges)


patient = "OO77"

#SNV
Conifer_output = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/Clustering",paste0(patient,".SCoutput.CLEAN.tsv"))
snv_df = read.delim(Conifer_output) 

snv_AD = snv_df %>% select(mutation_id,SAMPLE, VAR_COUNT) %>% pivot_wider(names_from = "SAMPLE", values_from = "VAR_COUNT") %>% as.data.frame
rownames(snv_AD) = snv_AD$mutation_id
snv_AD = snv_AD[,2:4]
snv_AD = as.matrix(snv_AD)

snv_DP = snv_df %>% select(mutation_id,SAMPLE, DEPTH) %>% pivot_wider(names_from = "SAMPLE", values_from = "DEPTH") %>% as.data.frame
rownames(snv_DP) = snv_DP$mutation_id
snv_DP = snv_DP[,2:4]
snv_DP = as.matrix(snv_DP)

#CNV
tumor = paste0(patient,"_tumor")
ex = paste0(patient, "_ex")
inv = paste0(patient, "_in")

Ascat_dir= file.path("/group/poetsch_projects/poetsch_sc/AscatNGS", patient)
tumor = read.delim(file.path(Ascat_dir, tumor, "result",paste0(patient,"_",tumor,".copynumber.caveman.txt")), header=F)
colnames(tumor) <- c("CHR","START","END","M","m")
tumor$M <- tumor$M - tumor$m

tumor_gr <- GRanges(
    seqnames = tumor$CHR,
    ranges = IRanges(tumor$START, end = tumor$END),
    sample = rep("tumor",nrow(tumor)),
    M = tumor$M,
    m= tumor$m)

ex = read.delim(file.path(Ascat_dir, ex, "result", paste0(patient,"_", ex,".copynumber.caveman.txt")), header=F)
colnames(ex) <- c("CHR","START","END","M","m")
ex$M <- ex$M - ex$m

ex_gr <- GRanges(
    seqnames = ex$CHR,
    ranges = IRanges(ex$START, end = ex$END),
    sample = rep("ex",nrow(ex)),
    M = ex$M,
    m= ex$m)

inv = read.delim(file.path(Ascat_dir,inv, "result", paste0(patient,"_", inv, ".copynumber.caveman.txt")), header=F)
colnames(inv) <- c("CHR","START","END","M","m")
inv$M <- inv$M  - inv$m

inv_gr <- GRanges(
    seqnames = inv$CHR,
    ranges = IRanges(inv$START, end = inv$END),
    sample = rep("in",nrow(inv)),
    M = inv$M,
    m= inv$m)

all_gr = c(tumor_gr,ex_gr,inv_gr)

disjoint_df = as.data.frame(disjoin(all_gr)) %>% select(-strand) %>% filter((end - start) != 0) %>% mutate(cnv_id   = paste0(sub("chr","",seqnames) ,":",  start  , ":",  end ))

#Annotate CNV by ClassifyCNV
disjoint_df_bed = disjoint_df %>% mutate(Annotate = "DEL") %>% select(seqnames ,  start   ,  end ,Annotate)
colnames(disjoint_df_bed) = NULL

ClassiyCNV_dir = file.path("/group/poetsch_projects/poetsch_sc/ClassifyCNV", patient, paste0(patient,"_disjoint/"))
dir.create(ClassiyCNV_dir)

write.table(disjoint_df_bed, file.path(ClassiyCNV_dir, paste0(patient,"_disjoint.bed")), quote=FALSE, row.names = F)


#Run bash 
module load apps/anaconda
conda activate ClassifyCNV
#bash update_clingen.sh
patient=OO77
main_dir=/group/poetsch_projects/poetsch_sc/ClassifyCNV/${patient}/${patient}_disjoint/
mkdir $main_dir/output
python3 /group/poetsch_projects/poetsch_sc/ClassifyCNV/ClassifyCNV.py  --infile $main_dir/${patient}_disjoint.bed --GenomeBuild hg38 --outdir $main_dir/output

#cd $output_dir

#awk '$6 == "Pathogenic"' Scoresheet.txt  > ${patient}_disjoint.Result.Pathogenic.txt


#### R ######
#reload result of ClassifyCNV
annotated_disjoint_df = read.delim(file.path("/group/poetsch_projects/poetsch_sc/ClassifyCNV", patient,paste0(patient,"_disjoint/output/Scoresheet.txt")), header=T)
annotated_disjoint_df = annotated_disjoint_df[,c(2,3,4,44)] 
colnames(annotated_disjoint_df) = c("Chr","Start","End","Gene")
annotated_disjoint_df = annotated_disjoint_df %>% mutate(cnv_id   = paste0(sub("chr","",Chr) ,":",  Start  , ":",  End ))


#Calculate Major and minor copy number of the disjoint granges
#reload disjoinr_df

CN_for_disjoint_seq <- function(new_seq, df){
    seqnames = disjoint_df$seqnames[new_seq]
    start = disjoint_df$start[new_seq]
    end = disjoint_df$end[new_seq]
    Major <- 1
    Minor <- 1
for (old_seq in 1:nrow(df)){
if (seqnames == df$CHR[old_seq]) {
     if ((start >= df$START[old_seq]) & (end <= df$END[old_seq]) ){
          Major = df$M[old_seq]
          Minor = df$m[old_seq]
          break
          } 
          }
          }
new_df = data.frame(chr = seqnames , start = start, end = end, M = Major, m = Minor)
new_df
}

tumor_new = purrr::map_dfr(1:nrow(disjoint_df),~CN_for_disjoint_seq(.x,df=tumor))
tumor_new = tumor_new %>% mutate( cnv_id = paste0(sub("chr","",chr) ,":",  start  , ":",  end ), sample = paste0(patient,"_tumor"))
tumor_new_M = tumor_new %>% select(cnv_id,sample, M)
tumor_new_m = tumor_new %>% select(cnv_id,sample, m)

ex_new = purrr::map_dfr(1:nrow(disjoint_df),~CN_for_disjoint_seq(.x,df= ex))
ex_new = ex_new %>% mutate( cnv_id = paste0(sub("chr","",chr) ,":",  start  , ":",  end ), sample = paste0(patient,"_ex"))
ex_new_M = ex_new %>% select(cnv_id,sample, M)
ex_new_m = ex_new %>% select(cnv_id,sample, m)



inv_new =   purrr::map_dfr(1:nrow(disjoint_df),~CN_for_disjoint_seq(.x,df= inv))
inv_new = inv_new %>% mutate( cnv_id = paste0(sub("chr","",chr) ,":",  start  , ":",  end ), sample = paste0(patient,"_in"))
inv_new_M = inv_new %>% select(cnv_id,sample, M)
inv_new_m = inv_new %>% select(cnv_id,sample, m)

#COmbine into M and m
M_data <- rbind(tumor_new_M,ex_new_M,inv_new_M)
M_data <- M_data %>% pivot_wider(names_from = "sample", values_from = "M")
row.names(M_data) <- M_data$cnv_id


m_data <- rbind(tumor_new_m,ex_new_m,inv_new_m)
m_data <- m_data %>% pivot_wider(names_from = "sample", values_from = "m")
row.names(m_data) <- m_data$cnv_id

  M_long <- M_data %>%
    pivot_longer(cols = -cnv_id, names_to = "sample", values_to = "M_value")
  
  m_long <- m_data %>%
    pivot_longer(cols = -cnv_id, names_to = "sample", values_to = "m_value")

#add gene annotated 
  combined_data <- M_long %>%
    left_join(m_long, by = c("cnv_id", "sample")) %>%
    pivot_wider(names_from = "sample", values_from = c("M_value", "m_value")) %>% left_join(annotated_disjoint_df[,4:5],by="cnv_id")

#Modify to export
# Create column names dynamically based on the patient_id
  columns_to_select = c(
                        paste0("M_value_", patient, "_tumor"),
                        paste0("M_value_", patient, "_ex"),
                        paste0("M_value_", patient, "_in"),
                        paste0("m_value_", patient, "_tumor"),
                        paste0("m_value_", patient, "_ex"),
                        paste0("m_value_", patient, "_in"))

combined_data_e = disjoint_df %>% left_join(combined_data, by = "cnv_id") %>% 
                      select(seqnames , start,  end, all_of(columns_to_select)) %>% 
                      pivot_longer(cols=columns_to_select, names_to = "Feature", values_to = "Copy_number" )

combined_data_e  <- combined_data_e %>%
  separate(Feature, into = c("Allele", "Sample"), sep = "_value_") %>% pivot_wider( names_from = "Allele", values_from = "Copy_number")

write.table(combined_data_e,file.path("/group/poetsch_projects/poetsch_sc/CANOPY",patient,"All_disjoint_seg_copy_number.txt"), quote = F, sep ="\t")

#Gastric drivers Totoki
Totoki = read.delim("/group/poetsch_projects/poetsch_sc/Driver_predict/Gastric_drivers_Totoki.txt")
driver_genes = Totoki$HUGO.Symbol

Find_driver_cn = function(driver_gene, combined_data){
  # Create the pattern to match the gene
  cnv = combined_data[grepl(paste0("\\b", driver_gene, "\\b"), combined_data$Gene),]
  cnv$driver_gene = driver_gene
  # Create column names dynamically based on the patient_id
  columns_to_select = c("driver_gene",
                        "cnv_id",
                        paste0("M_value_", patient, "_tumor"),
                        paste0("M_value_", patient, "_ex"),
                        paste0("M_value_", patient, "_in"),
                        paste0("m_value_", patient, "_tumor"),
                        paste0("m_value_", patient, "_ex"),
                        paste0("m_value_", patient, "_in"))
  # Select the dynamically created columns
  cnv = cnv %>% select(all_of(columns_to_select))
  
  return(cnv)
}


drivers_cn = purrr::map_dfr(driver_genes,~Find_driver_cn(.x,combined_data = combined_data))





# Function to classify CNV categories based on M_data and m_data
#Step1: identical M and m across 3 sample (group 1), identical within each pairs of samples (group 2,3,4), the rest(group X).
#Step2: in group X: identical m or M across 3 sample (group 5), identical m within each pair of samples (group 6,7,8), the rest(group 9).
# Return the cnv_id classified in each group 1-9

# Function to classify CNV categories based on M_data and m_data
classify_cnv <- function(M_data, m_data, patient) {
  # Reshape M_data and m_data to long format and combine
  M_long <- M_data %>%
    pivot_longer(cols = -cnv_id, names_to = "sample", values_to = "M_value")
  
  m_long <- m_data %>%
    pivot_longer(cols = -cnv_id, names_to = "sample", values_to = "m_value")
  
  combined_data <- M_long %>%
    left_join(m_long, by = c("cnv_id", "sample")) %>%
    pivot_wider(names_from = "sample", values_from = c("M_value", "m_value"))
  
  # Helper functions to check pairwise similarity
  identical_all <- function(values) {
    return(length(unique(values)) == 1)
  }

  identical_pair <- function(values1, values2) {
    return(length(unique(values1)) == 1 && length(unique(values2)) == 1)
  }
  
  # Create dynamic column names based on patient_id
  M_tumor_col <- paste0("M_value_", patient, "_tumor")
  M_ex_col <- paste0("M_value_", patient, "_ex")
  M_in_col <- paste0("M_value_", patient, "_in")
  m_tumor_col <- paste0("m_value_", patient, "_tumor")
  m_ex_col <- paste0("m_value_", patient, "_ex")
  m_in_col <- paste0("m_value_", patient, "_in")

  # Apply classification
  combined_data <- combined_data %>%
    rowwise() %>%
    mutate(group = case_when(
      # Step 1
      identical_all(c_across(all_of(c(M_tumor_col, M_ex_col, M_in_col)))) &
      identical_all(c_across(all_of(c(m_tumor_col, m_ex_col, m_in_col)))) ~ "group 1",
      
      identical_pair(c_across(all_of(c(M_tumor_col, M_ex_col))), c_across(all_of(c(m_tumor_col, m_ex_col)))) ~ "group 2",
      identical_pair(c_across(all_of(c(M_tumor_col, M_in_col))), c_across(all_of(c(m_tumor_col, m_in_col)))) ~ "group 3",
      identical_pair(c_across(all_of(c(M_ex_col, M_in_col))), c_across(all_of(c(m_ex_col, m_in_col)))) ~ "group 4",
      
      TRUE ~ "group X"
    )) %>%
    rowwise() %>%
    mutate(group_2 = case_when(
      # Step 2
      (group == "group X") &
      ((identical_all(c_across(all_of(c(m_tumor_col, m_ex_col, m_in_col))))) | 
      (identical_all(c_across(all_of(c(M_tumor_col, M_ex_col, M_in_col))))))  ~ "group 5",

      (group == "group X") &
      ((identical_all(c_across(all_of(c(m_tumor_col, m_ex_col))))) |
      (identical_all(c_across(all_of(c(M_tumor_col, M_ex_col)))))) ~ "group 6",

      (group == "group X") &
      ((identical_all(c_across(all_of(c(m_tumor_col, m_in_col))))) | 
      (identical_all(c_across(all_of(c(M_tumor_col, M_in_col)))))) ~ "group 7",

      (group == "group X") &
      ((identical_all(c_across(all_of(c(m_ex_col, m_in_col))))) |
      (identical_all(c_across(all_of(c(M_ex_col, M_in_col)))))) ~ "group 8",

      (group == "group X") ~ "group 9", 
      
      TRUE ~ group
    )) %>% 
    mutate(identical = case_when(
      (group_2 == "group 5") && (identical_all(c_across(all_of(c(m_tumor_col, m_ex_col, m_in_col))))) ~ "m",
      (group_2 == "group 5") && (identical_all(c_across(all_of(c(M_tumor_col, M_ex_col, M_in_col))))) ~ "M",
      (group_2 == "group 6") &&  (identical_all(c_across(all_of(c(m_tumor_col, m_ex_col)))))  ~ "m",
      (group_2 == "group 6") && (identical_all(c_across(all_of(c(M_tumor_col, M_ex_col))))) ~ "M",
      (group_2 == "group 7") &&   (identical_all(c_across(all_of(c(m_tumor_col, m_in_col))))) ~ "m",
      (group_2 == "group 7") &&   (identical_all(c_across(all_of(c(M_tumor_col, M_in_col))))) ~ "M",
      (group_2 == "group 8") &&  (identical_all(c_across(all_of(c(m_ex_col, m_in_col))))) ~ "m",
      (group_2 == "group 8") &&  (identical_all(c_across(all_of(c(M_ex_col, M_in_col))))) ~ "M",
      group_2 %in% c("group 1", "group 2", "group 3", "group 4") ~ "mM",
      TRUE ~ "No"
    ))
  return(combined_data)
}


classified_cnv <- classify_cnv(M_data, m_data, patient = patient)
classified_cnv <- classified_cnv %>% select(-group) %>% as.data.frame

#Add classify information to drivers_cn
 columns_to_select = c("driver_gene",
                        "cnv_id",
                        paste0("M_value_", patient, "_tumor"),
                        paste0("M_value_", patient, "_ex"),
                        paste0("M_value_", patient, "_in"),
                        paste0("m_value_", patient, "_tumor"),
                        paste0("m_value_", patient, "_ex"),
                        paste0("m_value_", patient, "_in"))

drivers_cn = drivers_cn %>% left_join(classified_cnv[,c(1,8)], by = "cnv_id") %>% 
                              select(driver_gene,cnv_id,group_2, all_of(columns_to_select))


#Find deletion in any sample of Totoki driver genes

# Create dynamic column names based on patient_id
  M_tumor_col <- paste0("M_value_", patient, "_tumor")
  M_ex_col <- paste0("M_value_", patient, "_ex")
  M_in_col <- paste0("M_value_", patient, "_in")
  m_tumor_col <- paste0("m_value_", patient, "_tumor")
  m_ex_col <- paste0("m_value_", patient, "_ex")
  m_in_col <- paste0("m_value_", patient, "_in")

drivers_cn = drivers_cn %>% mutate(total_cp_tumor = .data[[M_tumor_col]] + .data[[m_tumor_col]],
                                 total_cp_ex = .data[[M_ex_col]] + .data[[m_ex_col]],
                                 total_cp_in = .data[[M_in_col]] + .data[[m_in_col]]) %>% 
  distinct(driver_gene, .keep_all = TRUE)

CANOPY_dir = file.path("/group/poetsch_projects/poetsch_sc/CANOPY", patient)
dir.create(CANOPY_dir)
write.table(drivers_cn,file.path(CANOPY_dir,"Drivers_copy_number.txt"), sep = "\t", quote=F)

Deleted_drivers <- drivers_cn %>% 
  filter((total_cp_tumor == 0) | (total_cp_ex == 0) | (total_cp_in == 0) )  %>% as.data.frame()

#Find amplification in any Allele in any sample of Totoki driver genes

Amp_drivers <- drivers_cn %>% 
  filter((total_cp_tumor > 3) | (total_cp_ex > 3) | (total_cp_in > 3) ) %>% select(driver_gene,total_cp_tumor, total_cp_ex, total_cp_in)





#driver of OO77 predicted 
#Reload
Conifer_dir = file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient)
top_driver = read.delim(file.path(Conifer_dir , paste0(patient,"_Final_drivers.txt")))

pred_drivers = top_driver %>%  select(Gene_Annovar, mutation_id)

#Find pred_drivers in drivers_cn

pred_drivers_cn = drivers_cn %>% filter(driver_gene %in% pred_drivers$Gene) %>% as.data.frame

#Find copy number of mut
All_muts = read.delim(file.path(Conifer_dir,"/output/tree/Trees/treeTable.tsv"))

pred_muts = All_muts %>% filter(mutation_id %in% pred_drivers$mutation_id) %>% left_join( pred_drivers, by = "mutation_id")

pred_muts %>% filter(Gene_Annovar == "CDH1")