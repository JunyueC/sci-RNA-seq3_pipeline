
'''
This script accept a input sorted sam file, a output sam file, and a mismatch rate, then it will remove
duplicates based on the barcode + UMI (edit distance <= 1), and chromatin and start site, at the same
time, it will output the duplication number for each read, and generate the histogram plot for the read
per duplication number
'''
from Levenshtein import distance
import sys
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
sys.setrecursionlimit(1000000)
def rm_dup_samfile(samfile, output_file, mismatch):
    f1 = open(samfile)
    f2 = open(output_file, 'w')
    f3 = open(output_file+'.csv', 'w')
    pre_barcode = []
    pre_line = []
    unique_id = []
    pre_chrom = 0
    pre_site = 0
    dup = False
    mismatch=int(mismatch)
    
    pre_dup_num = 0
    cur_dup_num = 0
    
    counter = 0
    for line in f1:
        counter += 1
        if counter%100==0:
            print(counter)
        if (line[0] == '@'):
            f2.write(line)
        else:
            name = (((line.split('\t'))[0]).split(','))
            barcode_UMI = name[0] + name[1]
            chrom_num = (line.split('\t'))[2]
            start_site = (line.split('\t'))[3]
            
            if ((start_site == pre_site) and (chrom_num == pre_chrom)):
                if (barcode_UMI in pre_barcode):
                    cur_dup_num += 1
                else:
                    pre_barcode.append(barcode_UMI)
                    pre_line.append(line)
                '''
                dup = False
                for each_barcode in pre_barcode:
                    edit_dis = distance(each_barcode, barcode_UMI)
                    if edit_dis <= mismatch:
                        dup = True
                        break
                
                if dup == False:
                    pre_dup_num = cur_dup_num
                    cur_dup_num = 1
                    f2.write(line)
                    pre_barcode.add(barcode_UMI)
                    f3.write('%d' %(pre_dup_num))
                    f3.write('\n')
                else:
                    cur_dup_num += 1
                    if edit_dis == 0:
                        continue
                    else:
                        pre_barcode.add(barcode_UMI)
    '''
            else:
                if (pre_barcode != []):
                    unique_id = index_unique(pre_barcode, mismatch)
                    for i in unique_id:
                        f2.write(pre_line[i])
                    cur_dup_num = cur_dup_num + len(pre_barcode) - len(unique_id)
                
                pre_dup_num = cur_dup_num
                cur_dup_num = 1
                unique_id = []
                pre_barcode = []
                pre_line = []
                pre_chrom = chrom_num
                pre_site = start_site
                
                pre_barcode.append(barcode_UMI)
                pre_line.append(line)
                if (pre_dup_num != 0):
                    f3.write("%d" % (pre_dup_num))
                    f3.write('\n')
    
    # also count for the last reads
    if (pre_barcode != []):
        unique_id = index_unique(pre_barcode, mismatch)
        for i in unique_id:
            f2.write(pre_line[i])
        cur_dup_num = cur_dup_num + len(pre_barcode) - len(unique_id)
    
    f1.close()
    f2.close()
    f3.close()
    
    '''
    #plot the histogram for the read duplication number
    dups = (pd.read_csv(output_file+'.csv', header=None))[0]
    fig = plt.figure()
    plt.hist(dups, bins=100)
    plt.xlabel("Duplication number")
    plt.ylabel("Read number")
    fig.savefig(output_file + '.png')
    '''

def index_unique(barcodes_list, cutoff_value):
    # first I am going to create a matrix that include all the pairs of barcodes
    list_length = len(barcodes_list)
    distance_matrix = np.arange(list_length * list_length).reshape(list_length, list_length)
    
    for i in range(list_length):
        for j in range(list_length):
            if (j < i):
                distance_matrix[i][j] = distance_matrix[j][i]
            elif (j == i):
                distance_matrix[i][j] = 0
            else:
                distance_matrix[i][j] = distance(barcodes_list[i], barcodes_list[j])
    #print distance_matrix
    # Then I am going to create a list of indexs
    unique_index = []
    non_visited_indexes = range(list_length)
    # for each element in the index list, I am going remove its duplicates based on the cutoff_value
    for i in range(list_length):
        if i in non_visited_indexes:
            unique_index.append(i)
            rm_dup(non_visited_indexes, distance_matrix, i, cutoff_value)
    
    return unique_index

def rm_dup(indexes, distance_matrix, n, cutoff_value):
    # this script accept a index list, a distance matrix, and a index n; and then it remove the duplicates of index value n
    # and cutoff_value
    if n in indexes:
        indexes.remove(n)
    
    all_indexes = indexes[:]
    for i in range(len(all_indexes)):
        #print all_indexes, indexes, i, n, all_indexes[i]
        if all_indexes[i] in indexes:
            if (distance_matrix[all_indexes[i], n] <= cutoff_value):
                #print i, n, distance_matrix[all_indexes[i], n], indexes
                indexes = rm_dup(indexes, distance_matrix, all_indexes[i], cutoff_value)
    return indexes

if __name__ == "__main__":
    samfile = sys.argv[1]
    output_file = sys.argv[2]
    mismatch = sys.argv[3]
    rm_dup_samfile(samfile, output_file, mismatch)