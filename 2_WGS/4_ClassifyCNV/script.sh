conda activate ClassifyCNV

patients=("OO100" "OO77" "OO99")
treatments=("tumor" "ex" "in")

for patient in ${patients[@]}; do
for treatment in ${treatments[@]}; do

sample=${patient}_${treatment}


Main_dir=/group/poetsch_projects/poetsch_sc/ClassifyCNV/$patient/$sample
input_dir=$Main_dir/input
output_dir=$Main_dir/output
mkdir -p $input_dir
mkdir -p $output_dir

less /group/poetsch_projects/poetsch_sc/AscatNGS/$patient/$sample/result/${patient}_${sample}.copynumber.caveman.txt | \
awk '
{
  if (int($4) > 2) {
    category = "DUP"
  } else if (int($4) < 2) {
    category = "DEL"
  } else if (int($5) == 0) {
    category = "LOH"
  }
  if (category == "DUP" || category == "DEL") {
    print $1"\t"$2"\t"$3"\t"category
  }
}' | awk '$3 - $2 != 0' > $input_dir/$sample.copynumber.input.bed



bash update_clingen.sh


python3 /group/poetsch_projects/poetsch_sc/ClassifyCNV/ClassifyCNV.py  --infile $input_dir/$sample.copynumber.input.bed --GenomeBuild hg38 --outdir $output_dir

cd $output_dir

awk '$6 == "Pathogenic"' Scoresheet.txt  > $sample.Result.Pathogenic.txt
done
done
