#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=20-00:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 100 # how many processors
#SBATCH --mem=200g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J OO77_OO100 # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/Conifer/log    # folder for job error+logs.
#SBATCH -o OO77_OO100.log               # job log name
#SBATCH -e OO77_OO100.err               # job error output name (error output just means output without log flag - not all of it is errors)
#SBATCH --array=1-2



module load apps/anaconda
eval "$(conda shell.bash hook)"
conda activate pyclone


cd /group/poetsch_projects/poetsch_sc/Conifer/script/Vi_scripts/

patients=(OO77 OO100)
patient=${patients[$SLURM_ARRAY_TASK_ID-1]}

Rscript --vanilla Conipher_OO77_OO100.R $patient