
#!/bin/bash
# this scRNA-seq pipeline accept a input folder, and then use the default parameter for the data processing and analysis

# define the fastq folder including all fastq files
fastq_folder="/net/shendure/vol10/projects/scRNA/nobackup/171221_sci3_mouse_embryo/fastq/"
# define the output folder
all_output_folder="/net/shendure/vol10/projects/scRNA/nobackup/171221_sci3_mouse_embryo/output"
# define the PCR group id after demultiplexing
sample_ID="./sample_ID.txt"
# define the core number for parallele processing
core=15
# define the number of unique reads cutoff for splitting single cell
cutoff=200 # the number of unique reads cutoff for splitting single cell

# define the location for the index files used in STAR
index="/net/shendure/vol1/home/cao1025/../../../vol10/projects/scRNA/reference/index/STAR/STAR_mm10_RNAseq/"
# define the gtf file for gene counting
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/mm10/gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"

#define the mismatch rate for removing duplicates:
mismatch=1


# Define the location of the script folder
script_folder="/net/shendure/vol1/home/cao1025/analysis_script/sci3/"
# define the location of the ligation barcodes
ligation_barcode=$script_folder/lig_384_bc.pickle2
# define the location of the RT barcodes
RT_barcode=$script_folder//RT_384_bc.pickle2
# define the location of the combined RT and ligation barcodes
barcodes=$script_folder//combined_384_bc.txt
# define the location of the R script for multi-core processing
R_script=$script_path/sci3_bash_input_ID_output_core.R


#define the bin of python (python V2.7)
python_path="/net/shendure/vol1/home/cao1025/anaconda2/bin/"
#define the location of script:
script_path="/net/shendure/vol1/home/cao1025/analysis_script/sci3/"

now=$(date)
echo "Current time : $now"
module load samtools/1.3
module load bedtools/2.24.0

############ UMI attach
# this script take in a input folder, a sample ID, a output folder, a oligo-dT barcode file, a corresponding N5 barcode file, and
# it pass the factors to the python script
input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
script=$script_path/UMI_barcode_attach_gzipped_with_dic.py
echo "changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*r1.fq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*r2.fq.gz $input_folder/$sample.R2.fastq.gz; done

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

############align the reads with STAR, filter the reads based on q > 30, and remove duplicates based on exactly UMI sequence and tagmentation site
#define the output folder
input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

#align read2 to the index file using STAR with default setting
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

#make the filter sam folder, and filter and sort the sam file 
#make the flltered sam folder
echo
echo "Start filter and sort the sam files..."
echo input folder: $STAR_output_folder
echo output folder: $filtered_sam_folder
module load samtools/1.4
core=4
bash_script=$script_path/sci3_filter.sh
Rscript $R_script $bash_script $STAR_output_folder $sample_ID $filtered_sam_folder $core

core=20
# make a folder for rmdup_sam_folder, 
# Then for each filtered sam file, remove the duplicates based on UMI and barcode, chromatin number and position
echo
echo "Start removing duplicates..."
echo input folder: $filtered_sam_folder
echo output folder: $rmdup_sam_folder
mkdir -p $rmdup_sam_folder
module unload python

# bash_script=$script_path/sci3_rmdup_nomismatch.sh # for removing duplicates only considering exact match
bash_script=$script_path/sci3_rmdup.sh
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch

#mv the reported files to the report/duplicate_read/ folder
mkdir -p $input_folder/../report/duplicate_read
mv $rmdup_sam_folder/*.csv $input_folder/../report/duplicate_read/
echo "removing duplicates completed.."
echo
echo "Alignment and sam file preprocessing are done."  

################# split the sam file based on the barcode, and mv the result to the report folder
sam_folder=$all_output_folder/rmdup_sam
output_folder=$all_output_folder/sam_splitted

echo
echo "Start splitting the sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_ID
echo ouput folder: $output_folder
echo barcode file: $barcodes
echo cutoff value: $cutoff
module unload python

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
$python_path/python $script $gtf_file $input_folder $sample_ID $core

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