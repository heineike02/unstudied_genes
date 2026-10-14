import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil
import numpy as np

base_dir = os.path.normpath('/home/heineikeb/alphafold')
aln_dir = base_dir + os.sep + os.path.normpath('msas/structural/tm_align') 

align_files = os.listdir(aln_dir + os.sep + os.path.normpath('fasta_renamed'))
selected_alignments = [fname.split('.')[0] for fname in align_files]

mafft --auto bioinformatics/panther/SPBC115_03/concatenated_families.fasta > bioinformatics/panther/SPBC115_03/concatenated_families_aln.fasta

clipkit bioinformatics/panther/SPBC115_03/concatenated_families_aln.fasta -m gappy -g 0.1 -l

iqtree -s bioinformatics/panther/SPBC115_03/concatenated_families_aln.fasta.clipkit -nt 7 -bb 1000 -alrt 1000
#"AUTO"  automatically determines number of threads but 7 was performing well
#"-m", 'LG+I+G4',  #'MF', #only runs model finder 
                          #"-o", 'Spom_AF-Q10208-F1-model_v2'  #Outgroup for rooting should be pombe  for now using default. 
                         ]

trim_msa_thresh = 0.25  # Threshold to remove clusters that have poor alignments.  If the strict trimming MSA length is less than .25 * median sequence length, the cluster is removed. 

# As of 10 Jan 2023 that resulted in three OGs being removed. 
# Orthogroups filtered out when running 20221206_struct_align_dnds_msas.py because strict trimming of alignment was below trim_msa_thresh=0.25 * average sequence length threshold
# OG2147_REF_Scer_AF-P39692-F1-model_v2
# OG1306_REF_Scer_AF-P38298-F1-model_v2
# OG1746_REF_Scer_AF-P32642-F1-model_v2


ogs_filtered = []

#os.path.normpath('/home/heineike_wsl2/Crick_LMS/projects/diverse_yeasts/alphafold')
# output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')

#selected_ogs = ['OG2645']#['OG4150', 'OG2603', 'OG3677', 'OG2845']

#selected_og_refs = ['OG2645_REF_Scer_AF-P05375-F1-model_v2']#['OG4150_REF_Scer_AF-P07256-F1-model_v2', 'OG2603_REF_Scer_AF-P50076-F1-model_v2', 'OG2845_REF_Scer_AF-P43577-F1-model_v2', 'OG3677_REF_Scer_AF-P47125-F1-model_v2', 'OG1299_REF_Scer_AF-P00549-F1-model_v2']


## Get this to handle non REF ones or filter them out
#OG1111_alloascoidea_hylecoeti__OG1111__0_3867
#og_ref = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'

for alignment in selected_alignments: 
    #     print(alignment)
    #     if alignment.split('_')[1]=='REF':
    #         og,ref = alignment.split('_REF_')
    #     else: 
    og = alignment.split('_')[0]
    
    og_pep_msa_fname = aln_dir + os.sep + os.path.normpath('fasta_renamed/' + alignment + '.tm.fasta')
    
    #strict trimming for codon alignments    
    print('Strict trimming')
    clipkit_cmd = ['clipkit', og_pep_msa_fname, '-m','gappy' ,'-g','0.1', '-l']
    subprocess.run(clipkit_cmd)

    #Move clipkit log and output
    suffixes = ['.clipkit','.clipkit.log' ]

    for suffix in suffixes: 
        #Move clipkit files to new folder
        shutil.move(og_pep_msa_fname + suffix, 
                    aln_dir + os.sep + os.path.normpath('trim_strict/' + alignment + '.tm.fasta' + suffix)
                   )