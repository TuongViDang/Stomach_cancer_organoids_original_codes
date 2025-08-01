
#MUTAGENE
patient=OO100
vcf=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vcf

#transform vcf to maf
MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/mutagene
mkdir -p $MAIN_DIR/$patient
cd $MAIN_DIR



#install vcftomaf, install vep
conda activate vep

#mkdir -p $MAIN_DIR/.vep/homo_sapiens/102_GRCh38/
#rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz $MAIN_DIR/.vep/
#tar -zxf $MAIN_DIR/.vep/homo_sapiens_vep_102_GRCh38.tar.gz -C $MAIN_DIR/.vep/
#rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/homo_sapiens/dna_index/ $MAIN_DIR/.vep/homo_sapiens/102_GRCh38/



#create maf from vcf
perl vcf2maf-1.6.21/vcf2maf.pl --input-vcf $vcf --output-maf $MAIN_DIR/$patient/$patient.vep.maf --ncbi-build GRCh38 \
 --vcf-tumor-id ${patient}_${patient}_tumor \
 --vcf-tumor-id ${patient}_${patient}_ex \
 --vcf-tumor-id ${patient}_${patient}_in \
 --vcf-normal-id ${patient}_${patient}_normal \
  --vep-overwrite --vep-data $MAIN_DIR/.vep/ \
  --ref-fasta /group/poetsch_projects/poetsch_sc/nf-core_run/reference/Homo_sapiens_assembly38.fasta \
  --vep-path /home/biotec_poetsch/thda453f/.conda/envs/vep/bin/
 #--vcf-tumor-id ${patient}_${patient}_tumor2  

#ANNOTATE with VEP and ANNOVAR
#The above code will also create this vep annotated file

vcf_vep=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vep.vcf 

#extract function and gene name

grep -v "^#" $vcf_vep | awk -F"\t" -v OFS="\t" '{split($8, info, "|"); print $1, $2, $2, $4, $5, info[3], info[5]}' > /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vep_clean.vcf 

#Annotate with Annovar

ANNOVAR=/group/poetsch_users/poetsch_vi/annovar 
hummandb=$ANNOVAR/humandb
out_dir=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/Annovar
mkdir -p $out_dir

zless $vcf | grep -v '#' |  tr '\t' '-' | awk -F'-' -v OFS='\t' '{print $1,$2,$2,$4,$5,$0}' \
> $out_dir/${patient}.filtered.PASS.avinput


perl $ANNOVAR/table_annovar.pl $out_dir/${patient}.filtered.PASS.avinput \
                $hummandb -buildver hg38 \
                -out $out_dir/${patient}.filtered.PASS.variants.annovar \
                -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
                -nastring . -csvout -polish --otherinfo

less $out_dir/${patient}.filtered.PASS.variants.annovar.hg38_multianno.csv
awk -F, 'NR != 1 && $6 =="\"exonic\"" && $9 != "\"synonymous SNV\"" {print}' $out_dir/${patient}.filtered.PASS.variants.annovar.hg38_multianno.csv | wc -l #OO99: 154 missense, OO100:112, OO7OO7: 443




#run mutagene
conda activate mutagene
mutagene -v rank -i $MAIN_DIR/$patient/$patient.vep.maf -g hg38 -o $MAIN_DIR/$patient/$patient.ranking.txt -c pancancer 
mutagene -v rank -i $MAIN_DIR/$patient/$patient.vep.maf -g hg38 -o $MAIN_DIR/$patient/$patient.ranking_gastric.txt -c gastric_cancer 

#Combine Mutagen result with maf file 
### R
library(tidyverse)
patient="OO100"
MAIN_DIR="/group/poetsch_projects/poetsch_sc/Driver_predict/mutagene"

setwd(MAIN_DIR)
maf = read.delim(file.path(patient,paste0(patient,".vep.maf")),comment.char = "#")

cols <- c("Hugo_Symbol","Chromosome","Start_Position" ,"Reference_Allele" ,"Tumor_Seq_Allele2" ,"HGVSp_Short")
maf = maf[,cols]
maf$HGVSp_Short <- sub("p.","",maf$HGVSp_Short) 
maf$mutation_id <- paste0(":",sub("chr","",maf$Chromosome),":",maf$Start_Position,":",maf$Reference_Allele,":",maf$Tumor_Seq_Allele2)

mutagen <- read.delim(file.path(patient,paste0(patient,".ranking.txt")))

maf_mutagen <- left_join(maf,mutagen,join_by( HGVSp_Short == mutation, Hugo_Symbol == gene)) 
maf_mutagen<- maf_mutagen[,c("mutation_id","Hugo_Symbol","label")]
write.table(maf_mutagen,file.path(MAIN_DIR,patient,"maf_mutagen_result.txt"),sep="\t",quote=F)


mutagen_gastric <- read.delim(file.path(patient,paste0(patient,".ranking_gastric.txt")))

maf_mutagen_gastric <- left_join(maf,mutagen_gastric,join_by( HGVSp_Short == mutation, Hugo_Symbol == gene)) 
maf_mutagen_gastric <- maf_mutagen_gastric[,c("mutation_id","Hugo_Symbol","label")]
write.table(maf_mutagen_gastric,file.path(MAIN_DIR,patient,"maf_mutagen_gastric_result.txt"),sep="\t",quote=F)
