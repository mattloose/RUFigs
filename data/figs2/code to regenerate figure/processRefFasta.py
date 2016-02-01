#!/usr/bin/python
# -*- coding: utf-8 -*-
# --------------------------------------------------
# File Name: processRefFasta.py
# Purpose:
# Creation Date: 2014 - 2015
# Last Modified: Mon Jan 18 15:33:55 2016
# Author(s): The DeepSEQ Team, University of Nottingham UK
# Copyright 2015 The Author(s) All Rights Reserved
# Credits:
# --------------------------------------------------

import os
import sys
import re
from Bio import SeqIO
import numpy as np
import subprocess


# ---------------------------------------------------------------------------

# Naive z-score
def my_scale(a):  # MS
    mu = np.mean(a) # , None)
    sigma = np.std(a)
    return (a - mu) / sigma

import sklearn.preprocessing

def scale(a): 
	'''
	return my_scale(a)
	'''
	return sklearn.preprocessing.scale(a
		, axis=0, with_mean=True, with_std=False, copy=True) 
# ---------------------------------------------------------------------------

def process_ref_fasta(
    args,
    valid_ref_dir,
    bwa_index_dir,
    last_index_dir,
    ref_fasta,
    ref_fasta_hash,
    ):
    #print 'processing the reference fasta.'

    refdict = dict()
    refdict['seq_len'] = dict()
    refdict['refid'] = dict()
    refdict['seq_file'] = dict()
    refdict['seq_file_len'] = dict()

    # refdict["sequence"]=dict()

    refdict['kmer'] = dict()

    files = ref_fasta.split(',')
    ref_basename = ''

    if len(files) == 1:
        ref_basename = os.path.splitext(os.path.basename(files[0]))[0]

    if 1 < len(files):
        b = os.path.splitext(os.path.basename(files[0]))[0]
        ref_basename = '%s_plus_%s_more_seqs' \
            % (os.path.splitext(os.path.basename(files[0]))[0],
               str(len(files) - 1))

    validated_ref = os.path.join(os.path.sep, valid_ref_dir,
                                 ref_basename + '_valid.fasta')

    refdict['big_name'] = validated_ref
    refdict['big_len'] = 0

    # refdict["prefix"]=ref_basename

    refdict['path'] = os.path.dirname(files[0])

    if os.path.isfile(validated_ref) is False \
        or os.stat(validated_ref).st_size == 0:
        valid_fasta_handle = open(validated_ref, 'w')
        for fasta_file in files:
            #print 'FASTA file:', fasta_file
            fasta_records = list(SeqIO.parse(fasta_file, 'fasta'))
            if len(fasta_records) == 0:
                os.remove(validated_ref)
                err_string = \
                    "Error with your reference sequence FASTA file: %s: It's an empty file" \
                    % fasta_file
                #print >> sys.stderr, err_string
                sys.exit(1)

            try:
                for record in fasta_records:

                    # ref_fasta_hash["seq_file"][record.id]=os.path.splitext(os.path.basename(fasta_file))[0]

                    if len(record.seq) == 0:
                        os.remove(validated_ref)
                        err_string = \
                            'Error with your reference sequence FASTA file: %s: SEQID %s' \
                            % (fasta_file, record.id)
                        #print >> sys.stderr, err_string
                        sys.exit(1)
                    else:
                        seq = record.seq.upper()
                        record.seq = seq
                        record.description = \
                            os.path.basename(fasta_file)
                        SeqIO.write([record], valid_fasta_handle,
                                    'fasta')
            except Exception, err:

                os.remove(validated_ref)
                err_string = \
                    'Error with your reference sequence FASTA file: %s: %s' \
                    % (fasta_file, err)
                #print >> sys.stderr, err_string
                sys.exit(1)
        valid_fasta_handle.close()

    for record in SeqIO.parse(validated_ref, 'fasta'):
        if args.verbose is True:
            print 'processing seq: ', record.id

        refdict['seq_len'][record.id] = len(record.seq)
        disc = re.split('\s+', record.description)[1]
        refdict['seq_file'][record.id] = disc
        refdict['big_len'] += len(record.seq)

        # ####
        # refdict["sequence"][record.id]=record.seq
        # ####

        if disc in refdict['seq_file_len']:
            refdict['seq_file_len'][disc] += len(record.seq)
        else:
            refdict['seq_file_len'][disc] = len(record.seq)

        # ## do kmer

        if args.telem is True:
            recomp = record.seq.reverse_complement()
            km = kmer_count_fasta(record.seq, recomp, 5)
            refdict['kmer'][record.id] = km

    if args.last_align is True:
        last_index = os.path.join(os.path.sep, last_index_dir,
                                  ref_basename + '.last.index')
        last_index_file_part = os.path.join(os.path.sep,
                last_index_dir, ref_basename + '.last.index.bck')
        if os.path.isfile(last_index_file_part) is False \
            or os.stat(last_index_file_part).st_size == 0:
            #print 'Building LAST index for reference fasta...'
            cmd = 'lastdb -Q 0 %s %s' % (last_index, validated_ref)  # format the database
            if args.verbose is True:
                print cmd
            proc = subprocess.Popen(cmd, shell=True)
            status = proc.wait()
        refdict['last_index'] = last_index

    if args.bwa_align is True:
        bwa_index = os.path.join(os.path.sep, bwa_index_dir,
                                 ref_basename + '.bwa.index')
        bwa_index_file_part = os.path.join(os.path.sep, bwa_index_dir,
                ref_basename + '.bwa.index.bwt')
        if os.path.isfile(bwa_index_file_part) is False \
            or os.stat(bwa_index_file_part).st_size == 0:
            #print 'Building BWA index for reference fasta...'
            cmd = 'bwa index -p %s %s' % (bwa_index, validated_ref)  # format the database
            if args.verbose is True:
                print cmd
            proc = subprocess.Popen(cmd, shell=True)
            status = proc.wait()
        refdict['bwa_index'] = bwa_index

    #print 'finished processing reference fasta.'
    ref_fasta_hash[ref_basename] = refdict


# ---------------------------------------------------------------------------
'''
def process_ref_fasta_raw(ref_fasta, model_kmer_means):
    #print 'processing the reference fasta.'
    kmer_len = 7
    kmer_means = dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id] = dict()
        kmer_means[record.id]['F'] = list()
        kmer_means[record.id]['R'] = list()
        kmer_means[record.id]['Fprime'] = list()
        kmer_means[record.id]['Rprime'] = list()
        #print 'ID', record.id
        #print 'length', len(record.seq)
        #print 'FORWARD STRAND'

        seq = record.seq
        for x in range(len(seq) + 1 - kmer_len):
            kmer = str(seq[x:x + kmer_len])
            kmer_means[record.id]['F'
                                  ].append(float(model_kmer_means[kmer]))

            # if model_kmer_means[kmer]:
                # #print x, kmer, model_kmer_means[kmer]

        #print 'REVERSE STRAND'
        seq = revcomp = record.seq.reverse_complement()
        for x in range(len(seq) + 1 - kmer_len):
            kmer = str(seq[x:x + kmer_len])
            kmer_means[record.id]['R'
                                  ].append(float(model_kmer_means[kmer]))

        # @MS kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=False, copy=True)

        kmer_means[record.id]['Fprime'] = \
            scale(kmer_means[record.id]['F'])  # , axis=0, with_mean=True, with_std=False, copy=True)

        # @MS kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=False, copy=True)

        kmer_means[record.id]['Rprime'] = \
            scale(kmer_means[record.id]['R'])  # , axis=0, with_mean=True, with_std=False, copy=True)
    return kmer_means
'''

def process_ref_fasta_raw(ref_fasta, model_kmer_means):
    #print 'processing the reference fasta.'

    kmer_len = len(model_kmer_means.keys()[0]) # 7 # MS
    kmer_means = dict()
    for record in SeqIO.parse(ref_fasta, 'fasta'):
        kmer_means[record.id] = dict()
        kmer_means[record.id]['F'] = list()
        kmer_means[record.id]['R'] = list()
        kmer_means[record.id]['Fprime'] = list()
        kmer_means[record.id]['Rprime'] = list()
        #print 'ID', record.id
        #print 'length', len(record.seq)
        #print 'FORWARD STRAND'

        seq = record.seq
	'''
	print "For 5'->3':", seq[:20], seq[-20:]
	print "For 5'->3' 100:", seq[100:120]
	'''
        for x in range(len(seq) + 1 - kmer_len): 
            kmer = str(seq[x:x + kmer_len])
            kmer_means[record.id]['F'
                                  ].append(float(model_kmer_means[kmer]))

            # if model_kmer_means[kmer]:
                # #print x, kmer, model_kmer_means[kmer]

        #print 'REVERSE STRAND'
        seq = revcomp = record.seq.reverse_complement()
	refLen = len(seq)
	'''
	print "Rev 5'->3':" , seq[:20], seq[-20:]
	print "Rev 5'->3' 100:", seq[100:120]
	print "Rev 5'->3' 100:", seq[refLen-120:refLen-100]
	print "Rev 3'->5' 100:", seq[::-1][100:120]
	#print "===="
	'''
        for x in range(len(seq) + 1 - kmer_len):
            kmer = str(seq[x:x + kmer_len])
            kmer_means[record.id]['R'
                                  ].append(float(model_kmer_means[kmer]))

        kmer_means[record.id]['Fprime'] = \
            scale(kmer_means[record.id]['F']) 

        kmer_means[record.id]['Rprime'] = \
            scale(kmer_means[record.id]['R']) 
    return kmer_means



#---------------------------------------------------------------------------

def kmer_count_fasta(seq, revcompseq, kmer_len):
                kmerhash =dict()
                seqs = [seq, revcompseq]
                for x in range(len(seq)+1-kmer_len):
                                for s in seqs:
                                                kmer = str(s[x:x+kmer_len])
                                                if kmer in kmerhash:
                                                                kmerhash[kmer]+=1
                                                else:
                                                                kmerhash[kmer]=1
                return kmerhash
#---------------------------------------------------------------------------
