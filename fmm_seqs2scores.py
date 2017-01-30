#!/usr/bin/env python2

"""
Eilon Sharon, FMM 

"""

from __future__ import division

import sys
import os
import getopt
import argparse

import numpy as np

import xmltodict
import numpy as np
from Bio import SeqIO

class feature:
    
    def __init__(self, feature_dict):
        self.weight = float(feature_dict['@Weight'])

        self.pos = []
        self.letter_int = []
        
        if int(feature_dict['SequenceFeature']['@PositionsNum']) == 1:
            self.pos.append(int(feature_dict['SequenceFeature']['LettersAtPosition']['@Position']))
            self.letter_int.append(int(feature_dict['SequenceFeature']['LettersAtPosition']['@Letters']))
        else:
            for letteratpos in feature_dict['SequenceFeature']['LettersAtPosition']:
                self.pos.append(int(letteratpos['@Position']))
                self.letter_int.append(int(letteratpos['@Letters']))
    
    def cal_weight(self, dna_int, offset):
        is_feature = True
        for i,l in zip(self.pos,self.letter_int):
            if dna_int[offset+i] != l:
                is_feature = False
        
        if is_feature:
            return(self.weight)
        else:
            return(0.0)


class fmm:
    """
    fmm class, init from gxw(xml) file
    """
    
    def __init__(self, gxw_filename):
        with open(gxw_filename) as fd:
            doc = xmltodict.parse(fd.read())
        
        self.model = doc['WeightMatrices']['WeightMatrix']

        self.Zfunc = float(self.model['@LogPartitionFuncZ'])
        self.alphabet = self.model['@Alphabet']
        self.alphabet_size = int(self.model['@EffectiveAlphabetSize'])
        self.weights_num = int(self.model['@WeightsNum'])
        self.positions_num = int(self.model['@PositionsNum'])

        self.features = []

        for feature_dict in self.model['WeightedSequenceFeatureVecs']['SequenceFeatureVec']:
            self.features.append(feature(feature_dict))

    def dna2int(self,dna_seq):
        for idx,nt in enumerate(list(self.alphabet)):
            dna_seq = dna_seq.replace(nt, str(idx)) 
        dna_int = map(int,list(dna_seq))
        return(dna_int)
    
    def seq2scores(self,dna_seq):
        dna_int = self.dna2int(dna_seq)
        
        scores = []
        for offset in range(len(dna_seq)-self.positions_num):
            cur_score = -self.Zfunc
            for feature in self.features:
                cur_score += feature.cal_weight(dna_int, offset)
            scores.append(cur_score)
        
        return(scores)
            


def main():
    
    parser = argparse.ArgumentParser("Score sequences by FMM model. output is vectors of socres for each position")

    # input files
    parser.add_argument("input_gxw_filename", help= "FMM GXW (xml) file name, only firs model will be used")
    parser.add_argument("input_fa_filename", help= "sequences fasta files (only A,C,G,T upper case are allowed)")
    parser.add_argument("output_filename", help= "output file name")
    
    args = parser.parse_args()
    
    cur_fmm = fmm(args.input_gxw_filename)      
            
    fasta_sequences = SeqIO.parse(open(args.input_fa_filename),'fasta')
    with open(args.output_filename,'w') as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            print("scoring %s" %(name))
            out_file.write("%s = %s\n" % (name, str(cur_fmm.seq2scores(sequence))))



if __name__ == '__main__':
    main()
