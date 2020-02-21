
#!/bin/bash
# the scRNA-seq pipeline accept a input folder, and then use the default parameter for sequencing processing and generating the gene count matrix for downstream analysis

# The script is for processing sci-RNA-seq3 sequencing reads in UW genomic science cluster. For running in other environment, some modules, R and python packages are needed to be installed.

# define the fastq folder including all fastq files
# fastq_folder="/mnt/scratch/20191212-Pool-1058_sci-rna-seq/deml"
fastq_folder="/mnt/scratch/20191212-Pool-1058_sci-rna-seq/results_demux"

# define the PCR group sample id for each fastq file
sample_ID="/mnt/cbbi-data00/internal_data/20191212-Pool-1058_sci-rna-seq/samples.txt"

# define the output folder
# all_output_folder="/mnt/cbbi-data00/internal_data/20191212-Pool-1058_sci-rna-seq/sci3"
all_output_folder="/mnt/scratch/20191212-Pool-1058_sci-rna-seq/sci3"

# define the core number for parallele processing
core=15 # for most steps
samtools_core=4 # for reads filtering and sorting - this number is normally lower than the core number used in other script

# define the number of UMI cutoff for splitting single cell; cells with UMIs less than this number will be discarded
cutoff=200

# define the location of index files for reads alignment with STAR
index="/mnt/cbbi-data00/references/star/mm10"

# define the gtf file for gene counting
gtf_file="/mnt/cbbi-data00/references/gencode/gencode.vM15.annotation.gtf"

#define the mismatch rate for removing duplicates:
mismatch=1

# Define the location of the sub script folder
script_folder="/mnt/home/pjonsson/git/sci-RNA-seq3_pipeline/script_folder"

#define the bin of python (python V2.7)
python_path="/mnt/home/pjonsson/venv/scirnaseq/bin"

# load required modules from UW GS cluster
# module load modules modules-init modules-gs
# module load samtools/1.4
# module load STAR/2.5.2b
# module load python/2.7.3
# module load cutadapt/1.8.3
# module load trim_galore/0.4.1

# define the location of the ligation barcodes (they are in the script folder)
ligation_barcode=$script_folder/lig_384_bc.pickle2
# define the location of the RT barcodes
RT_barcode=$script_folder//RT_384_bc.pickle2
# define the location of the combined RT and ligation barcodes
barcodes=$script_folder//combined_384_bc.txt
# define the location of the R script for multi-core processing
R_script=$script_folder/sci3_bash_input_ID_output_core.R
script_path=$script_folder


now=$(date)
echo "Current time : $now"

############ UMI attach
# the script take an input folder, a sample ID list, an output folder, the RT barcode list, the ligation barcode list and core number. Then it extract the RT and ligation barcode from read1, correct them to the nearest RT and ligation barcode (with edit distance <= 1), and attach the RT and ligation barcode and UMI sequence to the read name of read2. Reads with unmatched RT or ligation barcodes are discarded.

input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
script=$script_path/UMI_barcode_attach_gzipped_with_dic.py
echo "Changing the name of the fastq files..."
# for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*r1.fq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*r2.fq.gz $input_folder/$sample.R2.fastq.gz; done

echo "Attaching barcode and UMI...."
mkdir -p $output_folder
$python_path/python $script $input_folder $sample_ID $output_folder $ligation_barcode $RT_barcode $core
echo "Barcode transformed and UMI attached."

################# Trimming the read2
echo
echo "Start trimming the read2 file..."
echo $(date)

trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach
bash_script=$script_path/sci3_trim.sh

Rscript $R_script $bash_script $UMI_attached_R2 $sample_ID $trimmed_fastq $core

############align the reads with STAR, filter the reads based on q > 30, and remove duplicates based on UMI sequence and tagmentation site

#define the output folder for mapping
input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

#align read2 to the index file using STAR
echo "Start alignment using STAR..."
echo input folder: $input_folder
echo sample ID file: $sample_ID
echo index file: $index
echo output_folder: $STAR_output_folder
#make the output folder
mkdir -p $STAR_output_folder
#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
#start the alignment
for sample in $(cat $sample_ID); do echo Aligning $sample;STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*gz --outFileNamePrefix $STAR_output_folder/$sample --genomeLoad LoadAndKeep; done
#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
echo "All alignment done."

#Filter and sort the sam file

echo
echo "Start filter and sort the sam files..."
echo input folder: $STAR_output_folder
echo output folder: $filtered_sam_folder
bash_script=$script_path/sci3_filter.sh
Rscript $R_script $bash_script $STAR_output_folder $sample_ID $filtered_sam_folder $samtools_core

# make a folder for rmdup_sam_folder, 
# Then for each filtered sam file, remove the duplicates based on UMI and barcode, chromatin number and position

# Remove duplicates based on UMI sequence (exact match) and tagmentation site

echo
echo "Start removing duplicates..."
echo input folder: $filtered_sam_folder
echo output folder: $rmdup_sam_folder
mkdir -p $rmdup_sam_folder
# module unload python

bash_script=$script_path/sci3_rmdup_nomismatch.sh # for removing duplicates only considering exact match
## bash_script=$script_path/sci3_rmdup.sh
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch

# repeat the rmdup process to remove duplicates based on edit distance of UMI sequence
echo
echo "Start removing duplicates..."
echo input folder: $all_output_folder/rmdup_sam
echo output folder: $all_output_folder/rmdup_sam_2
mkdir -p $all_output_folder/rmdup_sam_2
# module unload python

# bash_script=$script_path/sci3_rmdup_nomismatch.sh # for removing duplicates only considering exact match
bash_script=$script_path/sci3_rmdup.sh
filtered_sam_folder=$all_output_folder/rmdup_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam_2
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch

################# split the sam file based on the barcode, and mv the result to the report folder
sam_folder=$all_output_folder/rmdup_sam_2
output_folder=$all_output_folder/sam_splitted

echo
echo "Start splitting the sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_ID
echo ouput folder: $output_folder
echo barcode file: $barcodes
echo cutoff value: $cutoff
# module unload python

bash_script=$script_path/sci3_split.sh
Rscript $R_script $bash_script $sam_folder $sample_ID $output_folder $core $barcodes $cutoff


cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $output_folder/../barcode_samples.txt
# output the report the report/barcode_read_distribution folder
mkdir -p $output_folder/../report/barcode_read_distribution
mv $output_folder/*.txt $output_folder/../report/barcode_read_distribution/
mv $output_folder/*.png $output_folder/../report/barcode_read_distribution/
echo
echo "All sam file splitted."

################# gene count
# count reads mapping to genes
output_folder=$all_output_folder/report/human_mouse_gene_count/
input_folder=$all_output_folder/sam_splitted
script=$script_path/sciRNAseq_count.py
sample_ID=$all_output_folder/barcode_samples.txt
echo "Start the gene count...."
python $script $gtf_file $input_folder $sample_ID $core

echo "Make the output folder and transfer the files..."
mkdir -p $output_folder
cat $input_folder/*.count > $output_folder/count.MM
rm $input_folder/*.count
cat $input_folder/*.report > $output_folder/report.MM
rm $input_folder/*.report
mv $input_folder/*_annotate.txt $output_folder/
echo "All output files are transferred~"

now=$(date)
echo "Current time : $now"
