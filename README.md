# sciRNAseq3_process_pipeline

sci3_primer_sequences_plate.xls: primer sequences for reverse transcription, ligation, and PCR.

sci3_main.sh: main processing script for sci-RNA-seq3.

script_folder: folder for sub-scripts called by sci3_main.sh.

gene_count_processing_sciRNAseq.R: R script for processing the gene count data - the function of “sciRNAseq_gene_count_summary” accepts the gene count folder and then return a cell annotation data frame, a gene annotation data frame and a gene count sparse matrix. The output can be used as input to commonly used single cell RNA-seq analysis packages. 
