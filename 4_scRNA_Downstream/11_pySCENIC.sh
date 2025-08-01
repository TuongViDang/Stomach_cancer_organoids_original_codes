#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=600:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 50  # how many processors
#SBATCH --mem=200g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J pySCENIC_cont # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCENIC/log    # folder for job error+logs.
#SBATCH -o pySCENIC_cont.log               # job log name
#SBATCH -e  pySCENIC_cont.err               # job error output name (error output just means output without log flag - not all of it is errors)
 
module load apps/singularity
cd /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/SCENIC
# mkdir -p results

#singularity shell aertslab-pyscenic-0.12.1.sif


singularity exec aertslab-pyscenic-0.12.1.sif \
  pyscenic grn \
    expression_matrix.tsv \
    allTFs_hg38.txt \
    -o results/expr_mat.adjacencies.tsv  \
    --num_workers 50


# singularity exec aertslab-pyscenic-0.12.1.sif \
   pyscenic ctx \
   results/expr_mat.adjacencies.tsv \
   hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
   hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
  --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname expression_matrix.tsv\
  --output results/regulons.csv \
  --num_workers 50


singularity exec aertslab-pyscenic-0.12.1.sif \
  pyscenic aucell \
  expression_matrix.tsv \
  results/regulons.csv \
  --output results/auc_mtx.csv \
  --num_workers 50
