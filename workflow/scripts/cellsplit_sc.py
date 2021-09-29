#! /usr/bin/env python

import sys
import os
import pysam
import argparse
from collections import defaultdict, Counter

############################################
# Check if running as snakemake script
try:
    snakemake
    snakemake_flag = True
except:
    snakemake_flag = False
############################################
def split_bam_by_cell_barcode(bamfile, selected_barcodes_file, dest, log, barcode_tag, umicode_tag):

    import pandas as pd

    bam_in = pysam.AlignmentFile(bamfile) #creates AlignmentFile object
    bam_header = str(bam_in.header).strip() #get the header for the large bamfile
    file_handles = {} #dictionaries of filehandles
    reads_per_umis = defaultdict(set) #umis counter per cell
    reads_per_barcode = Counter() #count how many barcodes
    data_name,extension = os.path.basename(bamfile).split('.')
    selected_barcodes = set(pd.read_csv(selected_barcodes_file, header=None)[0].tolist()) #read the selected barcodes and turn them into a list



    if not os.path.exists(dest): #Make directories if they dont exists
        os.makedirs(dest) #if not, create corresponding directories

    for read in bam_in.fetch(until_eof=True):#For each read in bamfile
        if(read.has_tag(barcode_tag) and read.has_tag(umicode_tag)): # if the read has the selected barcode
            cbc = read.get_tag(barcode_tag) #get the barcode
            umi_code = read.get_tag(umicode_tag) #get the umicode

            if(cbc.split('-')[0] in selected_barcodes): #if the read is in the filtered barcodes file
                reads_per_barcode[cbc] += 1 #counter for reads per cellbarcode
                reads_per_umis[cbc].add(umi_code) #store the umi for the cellbarcode

                if(cbc not in file_handles):# if its not already created, create file handle
                    file_handle = os.path.join(dest, 'cbc_%s.%s.sam' % (cbc, data_name)) #create the file handle
                    file_handles[cbc] = file_handle #add the filehandle to dictionary of filehandles
                    #open the file and write the header
                    with open(file_handle, 'w') as f: #create and open the file
                        print(bam_header, file=f) #append the header
                        print(read.to_string(), file=f) #append the read
                        f.close() #close the file
                else:
                    with open(file_handles[cbc], 'a') as f: #if it already exists
                        print(read.to_string(), file=f) #print the aignment
                        f.close() #close the file


    ##########
    # Print report to the log to either file/sterr
    if log is None:
        outh = sys.stderr #print log to stderr
    else:
        outh = open(log, 'w') #if log exists, open it
    # Print log header
    print('barcode\treads_per_barcode\tumis_per_cell\tfile_path', file=outh)
    for barcode in file_handles.keys(): # log(3) : cell_barcode    number_of_reads    path_to_bam
        # print(barcode)
        n_umis_per_cell = len(set(reads_per_umis[barcode]))
        print(barcode+'\t'+str(reads_per_barcode[barcode])+'\t'+str(n_umis_per_cell)+'\t'+str(file_handles[barcode]),file=outh)
    #     # print('%s\t%d\t%d' % (barcode, float(n_reads_per_cell[barcode]), paths[barcode]), file=outh)
    # #Close the log
    if log is not None:
        outh.close()
###########################################
def count_barcodes(clustering_csv,bamfile):
    import pandas as pd
    import pysam

    df = pd.read_csv(clustering_csv)
    barcodes = df['Barcode']
    unique_barcodes = set(barcodes)
    print("Barcodes in clustering: "+str(len(barcodes)))
    print('Unique Barcodes in clustering: '+str(len(unique_barcodes)))

    barcodes_bam = []
    bam_in = pysam.AlignmentFile(bamfile, "rb") #reads bam file
    for read in bam_in.fetch(until_eof=True): #uses pysam to read bamfile line by line
        if(read.has_tag(barcode_tag)):
            cbc = read.get_tag(barcode_tag)
            barcodes_bam.append(cbc)

    unique_barcodes_bam = set(barcodes_bam)
    print('Barcodes in BAM: '+str(len(barcodes_bam)))
    print('Unique Barcodes in BAM: '+str(len(unique_barcodes_bam)))
############################################
def create_smaller_dataset(bamfile, nreads):

    print("#######################################")
    print("Creating smaller data set with "+str(nreads)+" reads")
    print("#######################################")
    bam_in = pysam.AlignmentFile(bamfile, "rb") #reads bam file
    bam_path = os.path.dirname(bamfile)+'/test_data.bam'
    bam_out = pysam.AlignmentFile(bam_path, 'wb', template=bam_in) #create the bamfile with the file handle
    counter = 0
    for read in bam_in.fetch(until_eof=True): #uses pysam to read bamfile line by line
        _ret = bam_out.write(read)
        counter += 1
        if(counter >= nreads):
            bam_out.close()
            break
    print('wrote '+str(counter)+' lines to: '+bam_path)

    i = 0
    print("---------------------------------")
    print("final file:")
    print(bam_path)
    print("---------------------------------")
    bam_in = pysam.AlignmentFile(bam_path, "rb") #reads bam file
    for read in bam_in.fetch(until_eof=True):
        i = i+1
        # print("READ: "+str(i))
        # print(read)
        # print("---------")
    print("File Ended with "+str(i)+'/'+str(counter)+" reads!")
############################################
def failed_bam_appending(bamfile, nreads):

    bam_in = pysam.AlignmentFile(bamfile, "rb") #reads bam file
    bam_path = os.path.dirname(bamfile)+'/test_data.bam'
    # bam_out = pysam.AlignmentFile(bam_path, 'wb', template=bam_in) #create the bamfile with the file handle
    counter = 0
    for read in bam_in.fetch(until_eof=True): #uses pysam to read bamfile line by line
        bam_out = pysam.AlignmentFile(bam_path, 'wb', template=bam_in) #create the bamfile with the file handle
        _ret = bam_out.write(read)
        counter += 1
        bam_out.close()
        if(counter > nreads):
            # bam_out.close()
            break

    print('---------------------------------')
    print('tried to write '+str(counter)+' times to: '+bam_path)
    # #Now, read the damn thing
    i = 0
    print("---------------------------------")
    print("final file:")
    print(bam_path)
    print("----------")
    bam_in = pysam.AlignmentFile(bam_path, "rb") #reads bam file
    for read in bam_in.fetch(until_eof=True):
        i = i+1
        print("READ: "+str(i))
        print(read)
        print("---------")
    print("File Ended with "+str(i)+'/'+str(counter)+" reads!")
############################################
def test_bam_appending(bamfile, nreads):

    bam_in = pysam.AlignmentFile(bamfile, "rb") #reads bam file
    bam_path = os.path.dirname(bamfile)+'/test_data.bam'
    counter = 0
    for read in bam_in.fetch(until_eof=True): #uses pysam to read bamfile line by line

        bam_out = pysam.AlignmentFile(bam_path) #create the bamfile with the file handle
        _ret = bam_out.write(read)
        counter += 1
        bam_out.close()
        if(counter > nreads):
            # bam_out.close()
            break

    print('---------------------------------')
    print('tried to write '+str(counter)+' times to: '+bam_path)
    # #Now, read the damn thing
    i = 0
    print("---------------------------------")
    print("final file:")
    print(bam_path)
    print("----------")
    bam_in = pysam.AlignmentFile(bam_path, "rb") #reads bam file
    for read in bam_in.fetch(until_eof=True):
        i = i+1
        print("READ: "+str(i))
        print(read)
        print("---------")
    print("File Ended with "+str(i)+'/'+str(counter)+" reads!")
############################################
def create_selected_barcodesfile(bamfile,n_barcodes):


    dest_folder = os.path.dirname(bamfile)
    barcodes_file = dest_folder+'/selected_barcodes.tsv'
    print(barcodes_file)
    bam_in = pysam.AlignmentFile(bamfile) #reads bam file
    i = 0
    barcodes = []
    barcode_tag = 'CB'
    f = open(barcodes_file, "w")
    for read in bam_in.fetch(until_eof=True): #uses pysam to read bamfile line by line
        if(read.has_tag(barcode_tag)):
            cbc = read.get_tag(barcode_tag).split('-')[0]
            if(cbc not in barcodes):
                barcodes.append(cbc)
                i+=1
            if(i >= n_barcodes):
                break
    #Write the barcode file
    for bc in barcodes:
        print(bc,file=f)
        # f.write(bc+'\n')
    f.close()
#########################################################################
# big_bam = '/home/santiago/github/data/bam/500_PBMC_3p_LT_Chromium_Controller_possorted_genome_bam.bam'
# test_bam = '/home/santiago/github/data/bam/test_data.bam'
# path_to_clustering = '/home/santiago/github/data/bam/clusters.csv'
# results_folder = '/home/santiago/github/data/bam/splitted_bam/'
# log = '/home/santiago/github/data/bam/splitted_bam/cellsplit_log.txt'
# selected_barcodes = '/home/santiago/github/data/bam/selected_barcodes.tsv'
# barcode_tag = 'CB' #cell barcode tag
# umicode_tag = 'UB'
# split_bam_by_cell_barcode(test_bam, selected_barcodes,results_folder, log, barcode_tag,umicode_tag)
# count_barcodes(path_to_clustering,path_to_bam)
# create_smaller_dataset(big_bam,100)
# test_bam_appending(test_bam,99)
# create_selected_barcodesfile(big_bam,100)
#########################################################################


# If the snakemake flag is true / run with snakemake
if snakemake_flag:

    if snakemake.params and 'barcode_tag' in snakemake.params:
        barcode_tag = snakemake.params.barcode_tag
    if snakemake.params and 'umicode_tag' in snakemake.params:
        umicode_tag = snakemake.params.umicode_tag
    #bamfile, selected_barcodes_file, dest, log, barcode_tag, umicode_tag
    split_bam_by_cell_barcode(snakemake.input.bam, snakemake.input.tsv, snakemake.output[0], snakemake.output[1],barcode_tag,umicode_tag)



# If not runned by snakemake

#If not snakemake get parameters from command line
elif __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Split a BAM file by cell barcode. It is assmued that the cell
                       barcode has been added to the read name'''
    )

    parser.add_argument("--bam",
        help="Path to BAM"
    )

    parser.add_argument("--brcds",
        help="Path to filtered barcodes from starsolo"
    )
    parser.add_argument("--dest",
        help='''Destination path. A BAM file for each cell barcode is created within the
                destination path (i.e. cbcXXXXXX). Default is in the same directory as
                the source BAM.'''
    )

    parser.add_argument("--log",
        help='Filename to write barcode count information. Default is STDERR'
    )

    parser.add_argument("--barcode_tag", default='CB',
        help='Delimiter to split cel_barcode . Default is "CB"'
    )

    parser.add_argument("--umicode_tag", default='UB',
        help='Delimiter to get the umi_code. Default is "UB"'
    )

    args = parser.parse_args()
    split_bam_by_cell_barcode(args.bam, args.brcds,args.dest, args.log,args.barcode_tag,args.umicode_tag)
