
patient=OO100
MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/Cscape
cd $MAIN_DIR
mkdir -p $MAIN_DIR/$patient


#input taken from dndscv
vcf=/group/poetsch_projects/poetsch_sc/Driver_predict/Cscape/${patient}/${patient}.lift_over_hg19.vcf

#submit online: http://cscape.biocompute.org.uk/

mv 7d0adbc1-6347-4fa3-8f96-c954f75e9541.tab ${patient}.lift_over_hg19_result.txt
less ${patient}.lift_over_hg19_result.txt

awk -F"\t"  '$7 == "oncogenic (high conf.)"' ${patient}.lift_over_hg19_result.txt | wc -l  
awk -F"\t" '$5 > 0.9'  ${patient}.lift_over_hg19_result.txt | less    | wc -l   #OO77: 86, OO99: 26, OO100: 18
awk -F"\t" '$6 > 0.9'  ${patient}.lift_over_hg19_result.txt | less    