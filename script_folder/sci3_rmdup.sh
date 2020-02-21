
input_folder=$1
sample=$2
output_folder=$3
mismatch=$4

# python="/net/shendure/vol1/home/cao1025/anaconda2/bin/python2.7"
python_script="/mnt/home/pjonsson/git/sci-RNA-seq3_pipeline/script_folder/rm_dup_barcode_UMI.py"

echo Filtering sample: $sample

python $python_script $input_folder/$sample.sam $output_folder/$sample.sam $mismatch

echo Filtering $sample done.