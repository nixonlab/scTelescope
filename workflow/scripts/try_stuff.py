#! /usr/bin/env python

import os


# bamfile = "Esta.Chingadera.Tiene.Unchingo.depuntos.bam"
# data_name = '.'.join(os.path.basename(bamfile).split('.')[0:-1])
# extension = os.path.basename(bamfile).split('.')[-1]

dest = '/home/santiago/github/nixonlab/scTelescope/workflow/scripts/results'
parent_dir = '/'.join(dest.split('/')[0:-1])
print(parent_dir)


#
# if not os.path.exists(dest): #Make directories if they dont exists
#     os.makedirs(dest) #if not, create corresponding directories
#
# fofn_filepath = os.path.pardir(os.path.abspath(dest)) #create file of filenames
# print(fofn_filepath)

# print(data_name)
# print(extension)
# exit()
