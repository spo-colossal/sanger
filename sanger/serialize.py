from Bio import SeqIO
import pickle 
import sys


#>TRPM8:[5180017,5192158]:[5168377,5168380]::55:5121806-5212478
class RefSeq:
    def __init__(self, header, seq):
        self.h_targets = header[0]
        self.l_targets = header[1] 
        self.scaffold = header[3] 
        self.range = header[4].split("-")
        self.seq  = seq



def main(pathname): 
    seqs = {}
    fastaFile = SeqIO.parse(open(pathname), 'fasta')
    for fasta in fastaFile: 
        header = fasta.id.split(":")
        data = RefSeq(header[1:], fasta.seq)
        seqs[header[0]] = data 
    with open('reference.pickle', 'wb') as handle: 
        pickle.dump(seqs, handle)
        


if __name__ == "__main__": 
    main("../mammoth.fasta")
