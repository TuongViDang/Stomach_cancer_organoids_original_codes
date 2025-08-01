#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=10:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 40  # how many processors
#SBATCH --mem=100g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J SCOMATIC_OO77_all  # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCOMATIC/log    # folder for job error+logs.
#SBATCH -o SCOMATIC_OO77_all.log               # job log name
#SBATCH -e SCOMATIC_OO77_all.err               # job error output name (error output just means output without log flag - not all of it is errors)


module load apps/anaconda
eval "$(conda shell.bash hook)"  #to use in batch mode 
conda activate SComatic 


#Detect Mutect2
#Prepare Mutect2 file
patient=$1
main_dir=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCOMATIC/$patient
bulk_var_dir=$main_dir/WGS_somatic_result
mkdir -p $bulk_var_dir
bulk_var_file=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vcf 


awk -F'\t' -v OFS='\t' '{print $1,$2,$2,$4,$5,".","Cell_type",".",".",".",".",".",".","."}' $bulk_var_file  > $bulk_var_dir/${patient}.Mutect.snv.indel.tsv


#Detect in scRNA-seq data

organoid_list=(${patient}-no ${patient}-vitro ${patient}-vivo )

for organoid in ${organoid_list[@]}; do

SCOMATIC=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCOMATIC/${organoid}

cd $SCOMATIC

sample=no_cluster

output_dir=$SCOMATIC/$sample/results
output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

META=$SCOMATIC/$sample/${sample}_Cell_type.tsv
REF=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa

#take BAM file
config=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/txt_files/gastric_meta.txt
library=$(awk -v organoid="$organoid" '$5 == organoid {print $1}' $config)
BAM=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/cellranger_count_outputs/$library/outs/possorted_genome_bam.bam

python $SCOMATIC/../scripts/SplitBam/SplitBamCellTypes.py --bam $BAM\
        --meta $META \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1


#assign mutation
bulk_var=${patient}_all
VAR_DIR=$main_dir/Mutect2_assign_to_cell/$organoid
mkdir -p $VAR_DIR

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
    
    temp=$VAR_DIR/temp_${bulk_var}_${cell_type}
    mkdir -p $temp

    python $SCOMATIC/../scripts/SingleCellGenotype/SingleCellGenotype.py --bam $bam  \
        --infile $bulk_var_dir/${patient}.Mutect.snv.indel.tsv \
        --nprocs 40  \
        --meta $META   \
        --outfile $VAR_DIR/${cell_type}.$bulk_var.single_cell_genotype.tsv  \
        --tmp_dir $temp  \
        --ref $REF

    rm -rf $temp
done

mv $VAR_DIR/${cell_type}.$bulk_var.single_cell_genotype.tsv  $VAR_DIR/$bulk_var.single_cell_genotype.tsv


#SNV
awk 'NR==1 || length($4) == length($5) && $5 == $10  ' $VAR_DIR/$bulk_var.single_cell_genotype.tsv\
 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$10}' > $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.tsv

#Indel
awk 'NR == 1 || ($10 == "D" && length($4) > length($5)) || ($10 == "I" && length($4) < length($5)) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$10}'  $VAR_DIR/$bulk_var.single_cell_genotype.tsv \
>  $VAR_DIR/$bulk_var.kept.single_cell_genotype.indel.tsv

#Combine
(cat $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.tsv ; grep -v "^#" $VAR_DIR/$bulk_var.kept.single_cell_genotype.indel.tsv ) > $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.tsv 

done 


#to count the unique mutation
cat $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.tsv  | grep -v "^#" | \
awk '{ key = $1"\t"$2"\t"$3"\t"$4"\t"$5; count[key]++} END { for (k in count) print k, count[k] }' | sort -k6nr | awk '{print $1"\t"$2"\t"$3}' > $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.unique_mut.tsv

cat $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.tsv | grep -v "^#" | \
awk '{ key = $1"\t"$2"\t"$3"\t"$4"\t"$5; count[key]++} END { for (k in count) print k, count[k] }' | wc -l   #equal number of unique detected mutations 

#to count the number of cell detected 
cat $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.tsv | grep -v "^#" | \
awk '{ key = $6; count[key]++} END { for (k in count) print k, count[k] }' | sort -k2nr | less

cat $VAR_DIR/$bulk_var.kept.single_cell_genotype.snv.indel.tsv | grep -v "^#" | \
awk '{ key = $6; count[key]++} END { for (k in count) print k, count[k] }' | sort -k2nr | wc -l    