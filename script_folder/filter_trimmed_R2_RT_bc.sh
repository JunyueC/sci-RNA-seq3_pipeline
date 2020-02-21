input_folder=$1
sample_ID=$2
output_folder=$3
RT_barcode_file=$4
core=$5

python="/net/shendure/vol1/home/cao1025/anaconda2/bin/python2"
script="/net/shendure/vol1/home/cao1025/analysis_script/sci3/filter_trimmed_R2_RT_bc.py"
mkdir -p $output_folder
$python $script $input_folder $sample_ID $output_folder $RT_barcode_file $core
echo All analysis done.