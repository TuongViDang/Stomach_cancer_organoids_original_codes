README
1_Alignment_nfcore_Sarek
- sample sheet, configuration file and the code to run sarek on biocluster 
  (note: run interactively on log-in node not compute node)
  
2_Mutect2
- Jointly call SNV from all samples of each patient
  
3_Ascat
- Call CNV calling with Ascat

4_ClassifyCNV
- Annotate the CNV

5_Driver_mutations_prediction
-  Mutagene code also give VEP and Annovar annotation
-  3_dndcsv is just to do the liftover hg38 to hg19 because Candra and Cscape use hg19
-  Cscape and Chasmplus: submit on their website and download the results as in the scripts
-  Run Driver_mutation_prediction.R to summarize the results

6_Trees_building_CONIPHER
- Use 1_Prepare_CONIPHER_input to turn Mutect2 and Ascat result into input of Conipher
- 2_Run_CONIPHER:
    1. Do the 1st Step: Clusterization => ".SCoutput.CLEAN.tsv" file
    2. Adapt some changes to get the ".SCoutput.CLEAN_final.tsv" before moving to 2nd step
    3. 2nd step: Tree building : use the ".SCoutput.CLEAN_final.tsv" file as input

7_Prepare_SBS_signature_extraction
   1. Use 1_split_each_clone.R to split mutations into clone
   2. Use 2_Subset_clone_specific_vcf.sh to extract the mutations in step 1 from the original vcf files. These new vcf files are input of SigProfiler.

8_SBS_signature_extraction
   1. Run SigProfiler Assignment witht the vcf files above with 1_SigProfiler.sh
   2. Go to 9_Assesment_SBS_and_filter_mutation:
      - run 1_Signature_Clone_OO77.R to assess and filter mutation likely to be artifact
      - Check context using 2_Mutational_Pattern_OO77.R
      - Do the filtering in the mutation list in the Conipher tree: 3_Filter_mutation_in_Conipher_clusters_OO77.R
      - prepare filtered vcf files to do SigprofierExtractor: transform_back_sigpro_filtered_files_to_vcf.sh
   3. Run SigProfilerExtractor using the filtered vcf files
  
9_Assesment_SBS_and_filter_mutation
    See 8
    Run 5_Signature_Clone_OO77_Filtered.R to have the final result
    
10_CloneMap_FishPlot
   To plot the clone composition changes 

11_Driver_CNV
   Use ClassifyCNV and Ascat result to obtain copy number of driver genes : 1_CANOPY_OO77.R
   Summarize the result and visualize with Copy_number.R 