README
1_scRNA_Preprocess_QC :  preprocess of scRNA-seq data

2_WGS : 
  - Aligment with nf-core sarek
  - SNV calling with Mutect2
  - CNV calling with Ascat
  - Tree-building with Conipher
  - Mutation signature extraction with SigProfiler
  - Driver prediction 

3_map_WGS_to_scRNA
  - Search mutation from bulk tree in scRNA-seq data using scReadCounts and SComatic
  - Clone assignment based on the mutation mapping

4_scRNA_Downstream
  - Downstream analysis that leads to plots

functions.R : a collection of functions used in some codes.