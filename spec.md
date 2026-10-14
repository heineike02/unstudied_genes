Now I want to extract the pLDDT data from the .cif files for all s. pombe genes and in the dataframe that contains that information, list whether genes are unstudied and also whether they are priority unstudied genes. 

For all S. pombe genes as listed in the file gene_IDs_names_products.tsv, please make a new ipython notebook called plddt_all_genes_vs_PUGs.ipynb that extracts gene_systematic_id, gene_name, and uniprot_id, and puts that into a dataframe.  Then for each uniprot_id, you should then extract the alphafoldDB id and the pLDDT score from the compressed .cif.gz file for the alphafolddb structure saved in G:\My Drive\Crick_LMS\projects\a_unstudied_genes\pombe_alphafold_structures\UP000002485_284812_SCHPO_v6\ with a filename whose structure is AF-<uniprot_id>-F1-model_v6.cif.gz.  

The alphafoldDB id is designated as _entry.id at the top of that file and the pLDDT score is a numeric value between 0 and 100 designated as _ma_qa_metric_global.metric_value near a line that reads 1 global pLDDT 1 pLDDT.  

Next, make a boolean column in the dataframe to label all unknown genes as defined in the file all_unstudied_genes.tsv, and another boolean column to label all priority unstudied genes as defined in the file priority_unstudied_genes.tsv

Finally, make a barplot comparing the average pLDDT score between all S. pombe genes, the priority unstudied genes, and the unstudied genes that are not priority unstudied genes.


In the cell labeled "Load unstudied gene classifications", 
when defining priority unstudied genes and unstudied genes the systematic id column is "Systematic ID" so there is no need to use different options 

I made a mistake in the instructions for the structure of the .cif.gz files.  The filename structure should be AF-<uniprot_id>-F1-model_v6.cif.gz.  