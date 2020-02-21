
input_folder=$1
sample=$2
output_folder=$3

# module load python/2.7.3
# module load cutadapt/1.8.3
# module load trim_galore/0.4.1
echo Trimming sample: $sample
trim_galore -j 8 $input_folder/$sample*.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $output_folder
# module unload python/2.7.3
echo Trimming $sample done.