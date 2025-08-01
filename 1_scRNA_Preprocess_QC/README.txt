README
Original directory: 
1_scRNA_Preprocess_QC:
    1_Cell_Ranger_count.sh: 
        - Fastqc files are processed by Cell Ranger 

    2_readfile_qc.R, 2_readfile_qc.sh
        - Output of Cell Ranger was imported to Seurat Object and QC and filtering are applied

    3_doublet_detect_remove.R, 3_doublet_detect_remove.sh
        - Doublet remove after QC using DoubletFinder