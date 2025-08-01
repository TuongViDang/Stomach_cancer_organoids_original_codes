
patient=OO100
vcf=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vcf 
#submit vcf file online https://run.opencravat.org/submit/nocache/index.html


cd /group/poetsch_projects/poetsch_sc/Driver_predict/Chasmplus/
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$12"\t"$14"\t"$15 }' 240507-095735_export_variant.tsv > result.clean.tsv

awk -F"\t" '$8 > 0.1 {print $5"\t"$7"\t"$8}' result.clean.tsv | less