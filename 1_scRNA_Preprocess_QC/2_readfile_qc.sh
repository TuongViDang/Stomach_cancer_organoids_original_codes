#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=24:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 1  # how many processors
#SBATCH --mem=200g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J Readfile_qc   # name of the job   
#SBATCH -D /group/poetsch_projects/Projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/log     # folder for job error+logs.
#SBATCH -o Readfile_qc.log               # job log name
#SBATCH -e Readfile_qc.err               # job error output name (error output just means output without log flag - not all of it is errors)
 
module load apps/R  # load a module
 
cd /group/poetsch_projects/Projects/poetsch_sc/scrna_gastric_cancer_230613/scripts/downstream_analysis

Rscript --vanilla readfile_qc.R
