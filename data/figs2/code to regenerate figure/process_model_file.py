import sys, os, re

from Bio import SeqIO
import csv 
import h5py
import numpy as np


def process_model_file(model_file):
    model_kmers = dict()
    with open(model_file, 'rb') as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        d = list(reader)
        for r in range(0, len(d)):
            #print r
            kmer = d[r][0]
            #print kmer
            mean = d[r][1] 
            #print type(mean)
            try:
                if (float(mean) <= 5): 
                    print "Looks like you have a poorly formatted model file. These aren't the means you are looking for.\n"
                    print "The value supplied for "+kmer+" was "+str(mean)
                    exit()
            except Exception,err:
                #print "Problem with means - but it isn't terminal - we assume this is the header line!"
		1
            model_kmers[kmer]=mean
    return     model_kmers
