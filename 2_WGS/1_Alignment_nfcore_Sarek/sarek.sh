#at login node

module load apps/nf-core
module load apps/singularity
module load apps/java/16.0.1

#nextflow pull nf-core/sarek -revision 3.2.1


cd /group/poetsch_projects/poetsch_sc/nf-core_run
mkdir -p results
nextflow run nf-core/sarek --input samplesheet_new.csv --outdir results --genome GATK.GRCh38 -profile singularity -revision 3.2.1 \
--tools mutect2,vep -c config_ori.txt  --save_output_as_bam -resume

