patient=OO100
MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/CanDrA
mkdir -p $MAIN_DIR/$patient

#use hg19 lift over vcf


awk -F"\t" '{print $1"\t"$2"\t"$4"\t"$5"\t""+"}' /group/poetsch_projects/poetsch_sc/Driver_predict/Cscape/${patient}/${patient}.lift_over_hg19.vcf > $MAIN_DIR/$patient/$patient.vep.CanDra.input.txt
input_file=$MAIN_DIR/$patient/$patient.vep.CanDra.input.txt

cd $MAIN_DIR

perl CanDrA.v1.0/open_candra.pl CRC  $input_file > $MAIN_DIR/$patient/$patient.result.txt

#awk -F"\t" '$12 == "Driver"' $MAIN_DIR/$patient/$patient.result.txt | less