
patients=("OO77" "OO99" "OO100")

for patient in ${patients[@]}; do

Main_dir=/group/poetsch_projects/poetsch_sc/Conifer/$patient/tree30/each_clone
cd $Main_dir
mkdir $Main_dir/vcf

SigProfiler_dir=/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/$patient/tree_30/vcf_SigProfiler/
rm -r $SigProfiler_dir
mkdir -p $SigProfiler_dir

clones=($(ls $Main_dir/mutation_files| grep -oP  '(?<=Mutations_).+(?=\.txt)'))
vcf=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/$patient.filtered.PASS.vcf 
zless $vcf | \
grep -v "^#"  > /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/$patient.filtered.PASS.no_header.txt

file_a=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/$patient.filtered.PASS.no_header.txt

for clone in ${clones[@]}; do
file_b=$Main_dir/mutation_files/Mutations_$clone.txt

(zless $vcf| grep "^#" ; \
 awk 'NR == FNR {seen[$1 "\t" $2]++; next} ($1 "\t" $2 in seen)' $file_b $file_a ) >  $Main_dir/vcf/${clone}_specific.mutation.vcf


grep -v "^#" $Main_dir/vcf/${clone}_specific.mutation.vcf | awk -v file=${clone} '{print $1"\t"$2"\t"file"\t"$4"\t"$5}' \
> $SigProfiler_dir/${clone}_specific.mutation.vcf

done

done
