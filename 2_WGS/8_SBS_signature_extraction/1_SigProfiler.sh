
#install Sigprofiler
patient=OO77
MAIN_DIR=/group/poetsch_projects/poetsch_sc/Driver_predict/SigProfilerExtractor/$patient/tree_30/
cd $MAIN_DIR


mkdir -p result_SBS
#mkdir -p result_DBS
#mkdir -p result_ID

python3
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh38')

from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples = "vcf_SigProfiler", 
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


Analyze.cosmic_fit(samples = "vcf_SigProfiler", 
                   output = "result_DBS", 
                   input_type="vcf", 
                   context_type="DINUC",
                   collapse_to_SBS96= False, 
                   cosmic_version=3.4,
                   exome=False,
                   genome_build="GRCh38",
                   export_probabilities_per_mutation = True
                   #,exclude_signature_subgroups = ['Artifact_signatures']
                   )


Analyze.cosmic_fit(samples = "vcf_SigProfiler", 
                   output = "result_ID", 
                   input_type="vcf", 
                   context_type="ID",
                   collapse_to_SBS96= False, 
                   cosmic_version=3.4,
                   exome=False,
                   genome_build="GRCh38",
                   export_probabilities_per_mutation = True
                   #,exclude_signature_subgroups = ['Artifact_signatures']
                   )


#Use SigProfilerExtractor
python3
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor(input_type = "vcf",
                         output = "result_Extractor",
                         input_data = "vcf_SigProfiler",
                         reference_genome = "GRCh38" ,
                         opportunity_genome = "GRCh38" ,
                        context_type = "ID",
                        cosmic_version = "3.4")


sig.sigProfilerExtractor(input_type = "vcf",
                         output = "result_Extractor_DBS",
                         input_data = "vcf_SigProfiler",
                         reference_genome = "GRCh38" ,
                         opportunity_genome = "GRCh38" ,
                        context_type = "DINUC",
                        cosmic_version = "3.4")

sig.sigProfilerExtractor(input_type = "vcf",
                         output = "result_Extractor_SBS",
                         input_data = "vcf_SigProfiler",
                         reference_genome = "GRCh38" ,
                         opportunity_genome = "GRCh38" ,
                        context_type = "96",
                        cosmic_version = "3.4")