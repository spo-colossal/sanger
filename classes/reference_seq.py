class RefSeq:
    def __init__(self, header, seq):
        self.h_targets = header[0]
        self.l_targets = header[1] 
        self.scaffold = header[3] 
        self.range = header[4].split("-")
        self.seq  = seq
