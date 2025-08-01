


#Prepare Mutect2 file
patients=(OO100 OO77)

for patient in ${patients[@]};do
MAIN_DIR=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scReadCounts/$patient
WGS_DIR=$MAIN_DIR/WGS_somatic_result
mkdir -p $WGS_DIR

bulk_var_file=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/CANOPY2/$patient/Somatic/MuTect2/${patient}.filtered.PASS.vcf 
less  $bulk_var_file| grep -v "^#" | awk ' {print $1" "$2" "$4" "$5}' >   $WGS_DIR/${patient}.filtered.PASS.txt

organoid_list=(${patient}-no ${patient}-vitro ${patient}-vivo )

for organoid in ${organoid_list[@]}; do
scRNA_DIR=$MAIN_DIR/$organoid

#move to the folder where barcodes.tsv file already available from Make_CB_file_from_seurat.R
cd $scRNA_DIR

#take the bam file
config=/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/Downstream_analysis/txt_files/gastric_meta.txt
library=$(awk -v org="$organoid" '$5 == org {print $1}' $config)

OUT_DIR=$scRNA_DIR/result
mkdir -p $OUT_DIR    

/group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/scReadCounts/SCReadCounts-1.3.2.Linux-x86_64/bin/scReadCounts \
                 -s  $WGS_DIR/${patient}.filtered.PASS.txt \
                 -r  /group/poetsch_projects/poetsch_sc/scrna_gastric_cancer_230613/cellranger_rnaseq/cellranger_count_outputs/$library/outs/possorted_genome_bam.bam \
                 -C 'CellRanger'\
                 -U 'CellRanger' \
                 -t 30 \
                 -o $OUT_DIR/$organoid.scReadCounts.tsv

awk 'NR ==1 || $10 > 0' $OUT_DIR/$organoid.scReadCounts.tsv > $OUT_DIR/$organoid.scReadCounts.PASS.tsv 

#count unique mut
cat  $OUT_DIR/$organoid.scReadCounts.PASS.tsv | grep -v "^#" | awk '{ key = $1"\t"$2"\t"$3"\t"$4; count[key]++} END { for (k in count) print k, count[k] }' | wc -l   #1053

#count unique cell
cat $OUT_DIR/$organoid.scReadCounts.PASS.tsv  | grep -v "^#" | awk '{ key = $5; count[key]++} END { for (k in count) print k, count[k] }' | sort -k2nr | wc -l   #4544
done

done