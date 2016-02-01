#!/usr/bin/python
# -*- coding: utf-8 -*-
# --------------------------------------------------
# File Name: mlpy_dtw_alignReal.py
# Purpose:
# Creation Date: 12-01-2016
# Last Modified: Tue Jan 19 14:01:33 2016
# Author(s): The DeepSEQ Team, University of Nottingham UK
# Copyright 2016 The Author(s) All Rights Reserved
# Credits: 
# --------------------------------------------------

import sys 
import numpy as np
from time import time as tm

# Unbuffered IO
import os
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

from os.path import isfile, join
def getFiles(path):
	return [f for f in os.listdir(path) if isfile(join(path, f))]

#------------------------------------------------------------------------------
# Get args ....

fasta_file = sys.argv[1]
template_model_file = sys.argv[2] 
complement_model_file = sys.argv[3] 
readsFolder = sys.argv[4] 
winSz = int(sys.argv[5])
frm = int(sys.argv[6])
to = int(sys.argv[7])
offset = int(sys.argv[8])
norm = sys.argv[9]

trg = frm+(to-frm)/2
ampSz = to-frm

file_leader = 'llssbzms2p35x_20151004_readuntiludududududu_RU21_lambdaPCR_2922_1_'
amplicon = "3" # (readsFolder.split('/')[-2])[9:]

#------------------------------------------------------------------------------
# Z-score Noramlisation ....

'''

def scale(a): 
    mu = np.mean(a, None)
    sigma = np.std(a)
    return (a - mu) / sigma
'''

#------------------------------------------------------------------------------
# Parse HDF5 Read File ....

import h5py

def process_hdf5(read_file):
        hdf = h5py.File(read_file, 'r')

	reads = 'Analyses/EventDetection_000/Reads/'
        for read in hdf[reads]:
                events = reads+ '/' + read + '/' + 'Events' 
                event_collection=[] 

                for event in hdf[events]:
                        event_collection.append(float(event['mean']))


                if hdf[ reads + read ].attrs['hairpin_found'] == 1:
			hairpinPos = hdf[ reads + read ].attrs['hairpin_event_index']
		else: hairpinPos = len(event_collection)

                read_t = event_collection[:hairpinPos]
                read_c = event_collection[hairpinPos+50:]

                return read_t, read_c

#------------------------------------------------------------------------------
# Parse and Process Model Files ....

# Get Models
from process_model_file import process_model_file

t_model = process_model_file(template_model_file)
c_model = process_model_file(complement_model_file)

#------------------------------------------------------------------------------
# Parse and Process Reference Files ....

# Get kmer means
from processRefFasta import process_ref_fasta_raw, scale

#print "t_kmers:"
t_kmer_means = process_ref_fasta_raw(fasta_file, t_model)
#print "c_kmers:"
c_kmer_means = process_ref_fasta_raw(fasta_file, c_model)

# Get Ref Forward Arrays ....
if norm is 'y': flag = 'prime'
else: 		flag = ''

for ks in t_kmer_means:
        ref_Ft = t_kmer_means[ks]['F' + flag]
for ks in t_kmer_means:
        ref_Fc = c_kmer_means[ks]['F' + flag]

# Get Ref Reversed Arrays ....
for ks in c_kmer_means:
        ref_Rt = t_kmer_means[ks]['R' + flag]
for ks in c_kmer_means:
        ref_Rc = c_kmer_means[ks]['R' + flag]

ref_arrays =  \
	[ (ref_Ft, ref_Rt) # Template strand: F and R models
	, (ref_Fc, ref_Rc) # Complement strand: F and R models
	]

refLen = len(ref_Ft)
#frm_ = refLen - frm
#to_ = refLen - to

#------------------------------------------------------------------------------
# DTW ...

import mlpy

def dtw(qry, (ref_F, ref_R)):

	if len(qry) is 0 : return '-1', -1

	dist_F, _, path_F = mlpy.dtw_subsequence(qry, ref_F)
	dist_R, _, path_R = mlpy.dtw_subsequence(qry, ref_R)

	if dist_F < dist_R: 	return 'F', path_F[1][0] - offset
	else: 			return 'R', refLen - path_R[1][0] - ampSz  + offset

#------------------------------------------------------------------------------
# Generate Synthectic Read...

def generateRead():
	sz =  ampSz
        trg = np.random.randint(0, refLen-sz, 1)[0]

	# Using t_model ...
	i,j = trg, trg+sz # NB ascending numerically
	#print trg, i, j
        read_t = ref_Ft[i:j] 

	# Using c_model 
	i,j = refLen - (trg+sz), refLen - trg  # NB decending numerically
	#print trg, i, j
        read_c = ref_Rc[i:j] 

        return trg, (read_t, read_c)

#------------------------------------------------------------------------------
# Process reads ....

if readsFolder == 'synthetic': 
	files = xrange(100) # Use Synthetic reads 
else:
	files = getFiles(readsFolder) # Use Real reads

for read_file in files: 


    if readsFolder == 'synthetic': 
    	fileId = "0, 0"
	trg, reads = generateRead() 
	frm = trg
	to = trg+ampSz
	#frm_ = refLen - frm
	#to_ = refLen - to
	#print "bounds:", frm,to, frm_,to_

    else: 
    	ch,fl = read_file[len(file_leader):].split('_')[:2]
    	ch_ = ch[2:]
    	fl_ = fl[4:]
    	fileId = ', '.join([ch_, fl_])
	reads = process_hdf5( readsFolder + '/' + read_file )

    res = []
    for label, (read, ref_F_R_pair) in \
			zip(['t','c'], zip(reads, ref_arrays)):

	qry = read[offset : offset+winSz]
	if len(qry) > 0 and norm == 'y' : qry = scale(qry) 

	tic=tm()
	dr, pos = dtw(qry, ref_F_R_pair) # Do DTW ...
	toc=tm()
	dt = round(toc - tic, 3)


	rangeTol = 100 # Given leader can affect results ...
	is_ok = frm - rangeTol <= pos <= to +rangeTol


	outBy =  abs(trg - pos) 

	hdr = [readsFolder, norm, fileId, amplicon ] \
		+ map(str, [winSz, len(read), len(qry), offset])
		

	res.append(
		( outBy
		, [label,dr] + map(str, [dt, pos , trg, outBy, is_ok ]) 
		) )

    (i,a),(j,b) = res
	
    tol = 50 # 1000 # winSz 
    dim = '0d'
    if i<tol or j<tol : dim = '1d'
    if i<tol and j<tol : dim = '2d'
    
    quasi2d = a[-1] == b[-1] == 'True' \
	      or (a[-1]=='True' and b[1]=='-1') # no nohairpin 

    print ', '.join(hdr +  a + b + [dim, str(quasi2d)])
    #print ', '.join(a + b + [dim, str(quasi2d)])
    sys.stdout.flush



