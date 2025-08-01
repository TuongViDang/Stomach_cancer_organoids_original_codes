#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=24:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 30 # how many processors
#SBATCH --mem=200g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J dndscv   # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/Driver_predict/dndscv/log     # folder for job error+logs.
#SBATCH -o dndscv.log               # job log name
#SBATCH -e dndscv.err               # job error output name (error output just means output without log flag - not all of it is errors)



#Load R
module load apps/R/4.3.1




patient=$1

MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/dndscv
cd $MAIN_DIR
mkdir -p $MAIN_DIR/$patient

vcf=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vcf 

(echo -e SampleID"\t"chr"\t"pos"\t"ref"\t"mut ; grep -v "^#" $vcf | awk -v patient=$patient '{print patient"\t"$1"\t"$2"\t"$4"\t"$5 }' ) \
 > $MAIN_DIR/$patient/${patient}.mutations.tsv

cd $MAIN_DIR/Vi_scripts

Rscript --vanilla dndscv.R $patient


