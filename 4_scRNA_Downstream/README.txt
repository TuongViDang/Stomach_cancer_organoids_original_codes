README
1_integration_3_patient_after_db_rm.R: integrate all samples  

2_after_QC.R: For Figure Supplementary 1A

3_Stomach_celltypes_markers.R: For Figure Supplementary 1B

4_by_patient_separate_in_ex_new.R: Plot all samples before integration and linear models for Figure 1

5_Transcription_dissimilarity_by_patient.R: For Figure 2A

6_Differential_vitro_vivo_each_patient_GO_LF0.5_p0.01.R: For Figure 2BC

7_pathway_score_NFKb_KS.R : Scoring HALLMARK pathways activities in single cell

8_Clustering_after_integration.R: Perform unsupervised clustering

9_cluster_marker.R: Find cluster markers

10_GO_enrichment_clusters.R: Gene Ontology form clusters'markers

11_pySCENIC.sh: Run 11_pySCENIC

12_pySCENIC_Downstream.R: Downtream analysis of pySCENIC results

13_NFKB_regulon_transcriptomic_clusters.R, 14_NFKB_regulon_expression_clusters.R: NFKB1 regulon activities in transcriptional clusters

15_Liana_run_each_sample.R: Run Liana in each sample and save as a list

16_Liana_post_analysis.R: Downstream analysis of Liana and plot Figure 7

17_NFKB_Liana.R: Plot a network of clusters and color with NFKB activities

18_ITH_genetic.R: Calculate genetic intra tumor heterogeneity (ITH)

19_ITH_transcriptomic.R: Calculate transcriptomic intra tumor heterogeneity (ITH)

20_Clone_Transcriptomic_VariantPartition.R: Linear model of gene expression with variate clone labels