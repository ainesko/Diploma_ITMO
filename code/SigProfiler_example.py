import SigProfilerExtractor
import pandas as pd
import numpy as np
import os 

# File should be in tsv format
#For SBS data
path_to_skin_table = '/Users/ksenia/SigProfilerExtractor/SigProfilerExtractor/data/TextInput/simulated_example.Skin.Melanoma.X.txt'
sig.sigProfilerExtractor("matrix", "SKCM_musical_final_results_SigProfiler", path_to_skin_table, opportunity_genome="GRCh38", minimum_signatures=1, maximum_signatures=15)
#For ID data
path_to_table = '/Users/ksenia/Desktop/ITMO/indels/Arm_A_Post.ID83.exome'
sig.sigProfilerExtractor("matrix", "/Users/ksenia/Desktop/ITMO/indels/new_SigProfiler_results/Arm_A_Post", path_to_table, opportunity_genome="GRCh38", context_type='ID', minimum_signatures=1, maximum_signatures=4, nmf_init="nndsvd")
