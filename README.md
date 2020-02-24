# sciRNAseq3_process_pipeline

## Pipeline steps

0. [Outside of pipeline] Demultiplex with bcl2fastq. After demultiplexing, there will be a pair of FASTQ files for each unique i5/i7 index (typically 96).
1. Extract cell barcode and UMI from read 1, add to read 2 FASTQ. Validate barcodes against list of RT and ligation barcodes used. By default allows for 1 mismatch _(I don't think this is actually implemented, although the variable is set to 1)_. Output in `UMI_attach` folder. 
2. Trim A+ sequence from FASTQ files using trim galore. Output in `trimmed_fastq` folder.
3. Align with STAR. Output in `STAR_alignment`.
4. Filter SAM files, hardcoded to retain reads with mapping quality >=30 and multimappers (`-q 30 -F 4`). Output in `filtered_sam`.
5. Remove duplicates—reads with the same UMI, barcode and genomic position. Output in `rmdup_sam`.
6. Repeat duplicate removal–this time include with up to 1 mismatch in barcode+UMI sequence. Output in `rmdup_sam_2`.
7. Splits the deduplicated SAM file into one file per valid barcode, at a given read cut-off. Output in `sam_splitted`.
8. Generate gene/exon counts using HTSeq. Output in `report/human_mouse_gene_count`.
9. [Outside of pipeline] Combine exonic and intronic reads and produce a count matrix across all cell barcodes and genes.

---
## Original Readme
sci3_primer_sequences_plate.xls: primer sequences for reverse transcription, ligation, and PCR.

sci3_main.sh: main processing script for sci-RNA-seq3.

script_folder: folder for sub-scripts called by sci3_main.sh.

gene_count_processing_sciRNAseq.R: R script for processing the gene count data - the function of “sciRNAseq_gene_count_summary” accepts the gene count folder and then return a cell annotation data frame, a gene annotation data frame and a gene count sparse matrix. The output can be used as input to commonly used single cell RNA-seq analysis packages. 
