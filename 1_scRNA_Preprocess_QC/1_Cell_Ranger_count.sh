#!/bin/bash
# =============================================================================
# Run cellranger count on 24 samples in cellranger/fastq folder
# =============================================================================
#
#SBATCH -p all.q             
#SBATCH --time=100:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -J run_cellranger_count_24_organoids_OO100-no
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --mem=64G
#SBATCH -D /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/bash_log       # folder for job error+logs.
#SBATCH -o run_cellranger_count_24_organoids_OO100-no.log
#SBATCH -e run_cellranger_count_24_organoids_OO100-no.err



cd /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/cellranger  #otherwise it looks for script in -D direction

module load apps/cellranger/7.1.0 

# Get list of runids
RUNID=($(ls /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/fastq/*_S1_L001_R1_001.fastq.gz | sed -n 's|.*/\(.*\)_.*_S1_L001_R1_001.fastq.gz|\1|p'))
unique_RUNID=($(printf "%s\n" "${RUNID[@]}" | sort -u))


fastq="/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/fastq"
reference="/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/reference/refdata-gex-GRCh38-2020-A"



for runid in "${unique_RUNID[0]}"; do
   sample=$(ls "$fastq" | grep -oE "$runid"_[0-9]+ | sort -u | paste -sd ",")
   echo $sample
   logfile="/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/log/$runid.cellranger_count.log"
   (rm -rf $runid
   cellranger count --id=$runid --fastqs=$fastq --description=$runid --transcriptome=$reference \
    --sample=$sample --include-introns=true --nosecondary --jobmode=local --localcores=8 --localmem=57
   find "$runid" -type f -user "$(whoami)" -exec chmod 660 {} \;
   find "$runid" -type d -user "$(whoami)" -exec chmod 770 {} \;
   mv "$runid" "/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/cellranger_count_outputs/") 2> "$logfile"
   echo $logfile

done
