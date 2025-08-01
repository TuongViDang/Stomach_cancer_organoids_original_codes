#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=60:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 40  # how many processors
#SBATCH --mem=100g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J Join_Mutect2 # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/logs    # folder for job error+logs.
#SBATCH -o Join_Mutect2.log               # job log name
#SBATCH -e Join_Mutect2.err               # job error output name (error output just means output without log flag - not all of it is errors)



#Load module
module load apps/gatk
module load apps/samtools

patient=OO99


MAIN_DIR=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic
mkdir -p $MAIN_DIR
cd $MAIN_DIR
[ ! -d bulk_align ] && mkdir -p bulk_align
[ ! -d bulk_align ] && mkdir -p scRNA_align


align_dir=$MAIN_DIR/bulk_align/$patient
mkdir -p $align_dir
cd $align_dir
samples=(${patient}_ex ${patient}_tumor ${patient}_tumor2 ${patient}_in  ${patient}_normal)

for sample in ${samples[@]};do
echo $align_dir/$sample.recal.bam
ln -s /group/poetsch_projects/poetsch_sc/nf-core_run/results/preprocessing/recalibrated/$sample/$sample.recal.bam  $align_dir/.
ln -s /group/poetsch_projects/poetsch_sc/nf-core_run/results/preprocessing/recalibrated/$sample/$sample.recal.bam.bai $align_dir/.
samtools view -bh $align_dir/$sample.recal.bam \
  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
  chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 \
  chrM chrX chrY > $align_dir/$sample.filter.recal.bam
done


for bam in $align_dir/*.filter.recal.bam ;do
samtools index $bam
done


MuTect2_dir=$MAIN_DIR/MuTect2
mkdir -p $MuTect2_dir


# http://www.bio-info-trainee.com/3386.html

tumor=${patient}_tumor
tumor2=${patient}_tumor2
ex=${patient}_ex
in=${patient}_in
normal=${patient}_normal
REF=/group/poetsch_projects/poetsch_sc/nf-core_run/reference/Homo_sapiens_assembly38.fasta

gatk Mutect2 -R $REF -I $align_dir/$tumor.filter.recal.bam  -I $align_dir/$tumor2.filter.recal.bam -I $align_dir/$ex.filter.recal.bam -I $align_dir/$in.filter.recal.bam -I $align_dir/$normal.filter.recal.bam -normal ${patient}_${patient}_normal -O $MuTect2_dir/${patient}.somatic.vcf.gz
	

#-------------------------------------------------------------
# Step 2c: FilterMutectCalls on joint calling vcfs
# FilterMutectCalls applies filters to the raw output of 
# Mutect2
#-------------------------------------------------------------


gatk FilterMutectCalls \
	-R $REF \
	-V $MuTect2_dir/${patient}.somatic.vcf.gz \
	-O $MuTect2_dir/${patient}.filtered.vcf.gz
	
gunzip -c $MuTect2_dir/${patient}.filtered.vcf.gz > $MuTect2_dir/${patient}.filtered.vcf
(grep "^#" $MuTect2_dir/${patient}.filtered.vcf; awk '$7 == "PASS"' $MuTect2_dir/${patient}.filtered.vcf)  > $MuTect2_dir/${patient}.filtered.PASS.vcf
grep -v "^#" $MuTect2_dir/${patient}.filtered.PASS.vcf | wc -l
