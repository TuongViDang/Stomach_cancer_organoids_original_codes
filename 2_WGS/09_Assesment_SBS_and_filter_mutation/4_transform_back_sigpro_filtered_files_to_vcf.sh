patient=OO77
Main_dir=/group/poetsch_projects/poetsch_sc/Conifer/$patient/tree30/each_clone
cd $Main_dir
rm -r $Main_dir/Filtered/vcf_filtered
mkdir -p $Main_dir/Filtered/vcf_filtered

vcf=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/$patient.filtered.PASS.vcf 

file_a=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/$patient.filtered.PASS.no_header.txt

Sigpro_dir=/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/$patient/tree_30/Filtered/vcf_filtered

clones=($(ls $Sigpro_dir/| grep "vcf" | sed 's/_filtered.vcf//'))

for clone in ${clones[@]}; do
file_b=$Sigpro_dir/${clone}_filtered.vcf

(zless $vcf| grep "^#" ; \
 awk 'NR == FNR {seen[$1 "\t" $2]++; next} ($1 "\t" $2 in seen)' $file_b $file_a ) >  $Main_dir/Filtered/vcf_filtered/${clone}_filtered.vcf

done