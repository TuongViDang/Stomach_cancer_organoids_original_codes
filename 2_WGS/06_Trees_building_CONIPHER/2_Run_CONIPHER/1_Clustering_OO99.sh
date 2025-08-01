#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=20-00:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 120 # how many processors
#SBATCH --mem=200g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J OO99_all # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/Conifer/log    # folder for job error+logs.
#SBATCH -o OO99_all.log               # job log name
#SBATCH -e OO99_all.err               # job error output name (error output just means output without log flag - not all of it is errors)
 
module load apps/anaconda
eval "$(conda shell.bash hook)"
conda activate pyclone


cd /group/poetsch_projects/poetsch_sc/Conifer/script/Vi_scripts/


Rscript --vanilla Conipher_OO99_all.R