from re import X
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq 
import sanger_seq
import reference_seq
import pickle
import os 
import sys
import pandas as pd 

WINDOW = 10
QUAL_SCORE = 20 
DATA_DIR = "data/"

#TODO:
# Flags:
# -phred
# -window 
# -tide or clone 
#------------------------------------------------------------------- 
# - find all targets that overlap with the alignment                |
# - determine how to convey off-target or somatic alteration        | 
# -                                                                 |
# -                                                                 |
# -------------------------------------------------------------------

RefSeq = reference_seq.RefSeq


def smith_waterman(): 
    

def process_ab1(pathname):
    sanger_data = []
    for _, _, files in os.walk(pathname): 
        for name in files:
            if ".ab1" in name:
                sanger_data.append(sanger_seq.SangerObject(os.path.join(pathname, name)))
    return sanger_data


def align(seqs):
    #load in all of the reference sequences 
    f = open("reference.pickle", 'rb') 
    ref = pickle.load(f)


    #load in the lookup table
    df = pd.read_csv('lookup.csv')
    conv_id = dict(zip(df['sample_id'], df['sample_description']))

    alignments = []
    for seq in seqs:
        gene_target = conv_id[seq.id.split("_")[0]]
        window = seq.find_alignable_window(WINDOW, QUAL_SCORE)
        if not window:
            continue
        if window['max_window'] is None: 
            continue
        start, stop = window['max_window']
        seq1 = seq.data[start:stop].upper()
        seq1_rc = Seq(seq1).reverse_complement()
        seq2 = ref[gene_target].seq.upper()
        a_f = pairwise2.align.localms(seq1, seq2, 2, -1, -.5, -.1, one_alignment_only=True)[0]  # read docs to parse object
        a_r = pairwise2.align.localms(seq1, seq2, 2, -1, -.5, -.1, one_alignment_only=True)[0]
        if a_f.score > a_r.score: 
            print(a_f.score)
            alignments.append([seq1, str(seq2[a_f.start:a_f.end])])
        else:
            print(a_r.score) 
            alignments.append([seq1_rc, str(seq2[a_r.start:a_r.end])])
    print(alignments)
    
def main(pathname):
    seqs = process_ab1(pathname)
    aigned_seq = align(seqs)


if __name__ == "__main__": 
    main(DATA_DIR)

