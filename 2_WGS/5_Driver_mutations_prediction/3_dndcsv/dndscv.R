
library(GenomicRanges)
library(rtracklayer)

#download and unzip the file "hg38ToHg19.over.chain" from http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/
#import chain file
MAIN_DIR="/group/poetsch_projects/poetsch_sc/Driver_predict/dndscv"
chain_file_38to19=import.chain(file.path(MAIN_DIR,"hg38ToHg19.over.chain"))

#or use system file
#library(liftOver)
#path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
#ch = import.chain(path)
patient <- commandArgs(trailingOnly = TRUE)[1]

position_file  <- read.delim(file.path(MAIN_DIR,patient,paste0(patient,".mutations.tsv")), header=T)


lift_over_38to19 <- function(i){ 
    chain_file = chain_file_38to19
    chr_hg38 = position_file[i, "chr"]
    pos_hg38 = position_file[i, "pos"]

    df <- data.frame(chr_hg38 = character(0),
                     pos_hg38 = numeric(0),
                     chr_hg19 = character(0),
                     pos_hg19 = numeric(0))
    df[1,"chr_hg38"] = chr_hg38
    df[1,"pos_hg38"] = pos_hg38
    # Create a GRanges object for the current position
    ranges <- GRanges(seqnames = chr_hg38, ranges = IRanges(start = pos_hg38, end = pos_hg38))
    
    # Perform liftOver
    lifted_ranges <- liftOver(ranges, chain_file)
    
    # If liftOver was successful and produced non-empty result
    if (length(lifted_ranges[[1]]) == 1 ) {
        # Extract the end position from the lifted result
         hg19_df<-as.data.frame(lifted_ranges)
         df[1,"chr_hg19"] <- sub("chr","",as.character(hg19_df$seqnames))
         df[1, "pos_hg19"] <- hg19_df$end
    } else {
        # If liftOver failed, set NA for hg19 position
        df[1,"chr_hg19"] <- NA
        df[1, "pos_hg19"] <- NA
    }
    df
    }

lift_overed_df<- purrr::map_dfr(1:nrow(position_file),lift_over_38to19)
write.table(lift_overed_df, file.path(MAIN_DIR,patient,paste0(patient,"_hg38to19_liftovered.txt")),sep="\t",quote=F)


#reload file
lift_overed_df <- read.delim(file.path(MAIN_DIR,patient,paste0(patient,"_hg38to19_liftovered.txt")))
mutations <- position_file
mutations$chr <- lift_overed_df$chr_hg19
mutations$pos <- lift_overed_df$pos_hg19
mutations <- mutations[!(is.na(mutations$pos)),]
mut_1= stringr::str_length(mutations$mut)==1
ref_1= stringr::str_length(mutations$ref)==1
mutations <- mutations[mut_1 & ref_1,]

library(dndscv)
dndsout <- dndscv(mutations)

sel_cv <- dndsout$sel_cv
significant = (sel_cv$pmis_cv < 0.05 | sel_cv$ptrunc_cv < 0.05 | sel_cv$pallsubs_cv < 0.05 | sel_cv$qmis_cv < 0.05 | sel_cv$qtrunc_cv < 0.05 | sel_cv$qallsubs_cv < 0.05)
sig_results <- sel_cv[significant,]

write.table(sig_results,file.path(MAIN_DIR,patient,"sig_result.txt"),sep="\t",quote=F)

#export mutatations to CScape input
mutations["note"] = "."
mutations=mutations[,c(2,3,6,4,5)]
colnames(mutations) <- NULL
rownames(mutations) <- NULL
Cscape_dir=file.path("/group/poetsch_projects/poetsch_sc/Driver_predict/Cscape",patient)
dir.create(Cscape_dir)
write.table(mutations,file.path(Cscape_dir,paste0(patient,".lift_over_hg19.vcf")), sep = "\t", quote = F, row.names=F)