#pull image
module load apps/singularity
#singularity pull docker://quay.io/wtsicgp/ascatngs:4.5.0


#run 

patient=OO77
samples=(${patient}_ex ${patient}_in ${patient}_tumor )

for sample in ${samples[@]}; do
OUT_DIR=/group/poetsch_projects/poetsch_sc/AscatNGS/$patient/$sample/result
mkdir -p $OUT_DIR
singularity exec \
/group/poetsch_projects/poetsch_sc/AscatNGS/image/ascatngs_4.5.0.sif \
ascat.pl \
   -r /group/poetsch_projects/poetsch_sc/nf-core_run/reference/Homo_sapiens_assembly38.fasta \
   -t /group/poetsch_projects/poetsch_sc/nf-core_run/results/preprocessing/recalibrated/$sample/$sample.recal.bam \
   -n /group/poetsch_projects/poetsch_sc/nf-core_run/results/preprocessing/recalibrated/${patient}_normal/${patient}_normal.recal.bam\
   -snp_gc /group/poetsch_projects/poetsch_sc/AscatNGS/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv \
   -pr WGS \
   -g XX \
   -gc chrY \
   -rs "Homo sapiens" \
   -ra GRCh38 \
   -pl ILLUMINA \
   -q 20 \
   -c 45 \
   -o $OUT_DIR
done
