"""
Created on Tue Jan  16 10:58:54 2018

@author: Junyue
"""

import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial

'''
    this script accept a read2 file and a RT barcode list and then extract the reads matching the RT barcodes
'''
    
def filter_RT(sample, input_folder, output_folder, RT_barcode_list):
    #open the read1, read2, and output file
    Read1 = input_folder + "/" + sample + ".R2_trimmed.fq.gz"
    output_file = output_folder + "/" + sample + ".R2.fq.gz"

    f1 = gzip.open(Read1)
    f3 = gzip.open(output_file, 'wb')
    
    line1 = f1.readline()
    total_line = 0
    filtered_line = 0

    while (line1):
        total_line += 1
        line1_bc = line1.split(",")[0]
        
        find = False
        for barcode in RT_barcode_list:
            if (barcode in line1_bc):
                find = True
                break
        if (find == True):
            filtered_line += 1
            f3.write(line1)
            
            line1 = f1.readline()
            f3.write(line1)
            
            line1 = f1.readline()
            f3.write(line1)
            
            line1 = f1.readline()
            f3.write(line1)
            
            line1 = f1.readline()

        else:
            line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()
    f1.close()
    f3.close()
    
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f" 
          %(sample, total_line, filtered_line, float(filtered_line) / float(total_line)))
        

# this function accept an input folder and a output folder and then generate the output file with the index
def filter_RT_files(input_folder, sampleID, output_folder, RT_barcode_file, core):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    RT barcode file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, RT_barcode_file)
    
    print(init_message)
    
    # generate the RT barcode list:
    RT_barcode_list = []
    barcodes = open(RT_barcode_file)
    for barcode in barcodes:
        RT_barcode_list.append(barcode.strip())
    barcodes.close()
    
    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core))
    #print("Processing core number: ", core_number)
    func = partial(filter_RT, input_folder = input_folder, output_folder=output_folder, RT_barcode_list=RT_barcode_list)
    #sciRNAseq_count(sample, input_folder, exons, genes, gene_end)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    RT_barcode_file = sys.argv[4]
    core=sys.argv[5]
    filter_RT_files(input_folder, sampleID, output_folder, RT_barcode_file, core)