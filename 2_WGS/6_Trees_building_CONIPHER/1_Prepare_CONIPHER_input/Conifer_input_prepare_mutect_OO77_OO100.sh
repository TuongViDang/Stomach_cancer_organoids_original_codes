#!/bin/bash
 
#SBATCH -p all.q             
#SBATCH --time=5-00:00:00 ## Time until auto-shutoff                         
#SBATCH --mail-user=thi_tuong_vi.dang@tu-dresden.de ## Change the email address
#SBATCH --no-requeue   ## You can remove it - with it on, if it fails, it can't restart
#SBATCH -n 4  # how many processors
#SBATCH --mem=100g  # how much memory
#SBATCH -N 1   # how many nodes 
#SBATCH --mail-type=ALL   # how much info do you want to get by mail (ie - start time, end time+endstate)
#SBATCH -J OO77_OO100_prepare  # name of the job   
#SBATCH -D /group/poetsch_projects/poetsch_sc/Conifer/log    # folder for job error+logs.
#SBATCH -o OO77_OO100_prepare.log               # job log name
#SBATCH -e OO77_OO100_prepare.err               # job error output name (error output just means output without log flag - not all of it is errors)
#SBATCH --array=1-2

patients=(OO77 OO100)

patient=${patients[$SLURM_ARRAY_TASK_ID-1]}

Mutect_dir=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2
Mutect_file=$Mutect_dir/${patient}.filtered.PASS.vcf

Main_dir=/group/poetsch_projects/poetsch_sc/Conifer/$patient
Input_dir=$Main_dir/data
mkdir -p $Input_dir
cd $Main_dir

grep -v "^#" $Mutect_file | awk -F'\t' -v patient=$patient  '{split($10, LF, ":"); split(LF[2],AD,","); print patient"\t"patient"_ex""\t"$1"\t"$2"\t"$4"\t"$5"\t"AD[1]"\t"AD[2]"\t"LF[4]}' \
  > $Input_dir/${patient}_ex.mutect.tsv

grep -v "^#" $Mutect_file | awk -F'\t'  -v patient=$patient  '{split($11, LF, ":"); split(LF[2],AD,","); print patient"\t"patient"_in""\t"$1"\t"$2"\t"$4"\t"$5"\t"AD[1]"\t"AD[2]"\t"LF[4]}' \
  > $Input_dir/${patient}_in.mutect.tsv

grep -v "^#" $Mutect_file | awk -F'\t'  -v patient=$patient '{split($13, LF, ":"); split(LF[2],AD,","); print patient"\t"patient"_tumor""\t"$1"\t"$2"\t"$4"\t"$5"\t"AD[1]"\t"AD[2]"\t"LF[4]}' \
  > $Input_dir/${patient}_tumor.mutect.tsv


samples=(${patient}_ex ${patient}_in ${patient}_tumor )


for sample in ${samples[@]}; do

mutect_var=$Input_dir/$sample.mutect.tsv
ascat_dir=/group/poetsch_projects/poetsch_sc/AscatNGS/$patient/$sample/result

#Ascat segmentation copy number
less  $ascat_dir/${patient}_${sample}.copynumber.caveman.csv | awk -F"," '{print $2"\t"$3"\t"$4"\t"$7"\t"$8}' >  $ascat_dir/${patient}_${sample}.copynumber.caveman.txt

ascat_seg=$ascat_dir/${patient}_${sample}.copynumber.caveman.txt

purity=$(echo "1 - $(grep "NormalContamination" "$ascat_dir/${patient}_${sample}.samplestatistics.txt" | awk -F' ' '{print $2}')" | bc)
ploidy=$(grep "Ploidy" "$ascat_dir/${patient}_${sample}.samplestatistics.txt" | awk -F' ' '{print $2}')

#generate the input.tsv file
python /group/poetsch_projects/poetsch_sc/Conifer/script/Vi_scripts/add_ascat_to_variant.py  $ascat_seg $mutect_var $purity $ploidy
mv modified_mutect_file_ascat.txt $Input_dir/$sample.mutect.ASCAT.tsv


#add the header
echo -e "CASE_ID\tSAMPLE\tCHR\tPOS\tREF\tALT\tREF_COUNT\tVAR_COUNT\tDEPTH\tCOPY_NUMBER_A\tCOPY_NUMBER_B\tACF\tPLOIDY" > header.txt

cat header.txt $Input_dir/$sample.mutect.ASCAT.tsv \
 > $Input_dir/$sample.mutect.ASCAT.header.tsv
rm header.txt
rm $Input_dir/$sample.mutect.ASCAT.tsv


#change CHR from string to integer 
(awk 'NR == 1 {print $0}' $Input_dir/$sample.mutect.ASCAT.header.tsv ;
 awk 'BEGIN{OFS="\t"} NR>1{$3=substr($3,4); if ($3 == "X") $3 = 23; else if ($3 == "Y") $3 = 24; else $3 = int($3); print }' $Input_dir/$sample.mutect.ASCAT.header.tsv ) \
    > ${Input_dir}/${sample}.mutect.ASCAT.header.conipher_input.tsv

#generate the input_seg.tsv

less "$ascat_seg" | awk -F'\t' "{printf \"%s\t%s\n\", \"$sample\", \$0}" > "$ascat_dir/${patient}_${sample}.copynumber.caveman.seg.tsv"


echo -e "SAMPLE\tCHR\tSTARTPOS\tENDPOS\tCOPY_NUMBER_A\tCOPY_NUMBER_B" > header.txt

cat header.txt $ascat_dir/${patient}_${sample}.copynumber.caveman.seg.tsv  | awk 'BEGIN{OFS="\t"} NR>1{$2=substr($2,4)}1' |\
awk '$2 !~ /^(X|Y)$/'  \
> $Input_dir/${patient}_${sample}.copynumber.caveman.seg.header.tsv

rm header.txt

done


#Combined inout file and seg file
awk 'FNR == 1 && NR > 1 {next} {print}' $Input_dir/*.mutect.ASCAT.header.conipher_input.tsv  | \
awk 'NR==1; NR>1 {print $0 | "sort -k3,3n -k4,4n"}' > $Input_dir/$patient.All.mutect.ASCAT.header.conipher_input.tsv


awk 'FNR == 1 && NR > 1 {next} {print}' $Input_dir/*.copynumber.caveman.seg.header.tsv | \
awk 'NR==1; NR>1 {print $0 | "sort -k2,2n -k3,3n"}' > $Input_dir/$patient.All.copynumber.caveman.seg.header.tsv 
