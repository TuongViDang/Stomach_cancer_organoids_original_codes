
#install Sigprofiler
patient=OO77
MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/$patient/tree_30/Filtered
cd $MAIN_DIR


mkdir -p result_SBS


python3
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh38')

from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples = "vcf_filtered", 
                   output = "result_SBS", 
                   input_type="vcf", 
                   context_type="96",
                   collapse_to_SBS96=True, 
                   cosmic_version=3.4,
                   exome=False,
                   genome_build="GRCh38" ,
                   export_probabilities_per_mutation = True
                  #,exclude_signature_subgroups = ['Artifact_signatures']
                   )
python3
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor(input_type = "vcf",
                         output = "result_Extractor",
                         input_data = "vcf_filtered",
                         reference_genome = "GRCh38" ,
                         opportunity_genome = "GRCh38" ,
                        context_type = "96",
                        cosmic_version = "3.4")