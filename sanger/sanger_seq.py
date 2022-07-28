import os
import numpy as np
from Bio import SeqIO, Align

class SangerObject:

    def __init__(self, pathname):
        self.id = os.path.basename(pathname)
        self.rev = (self.id.split("_")[1][-1] == "R")
        print(self.id.split("_")[1][-1])
        record = SeqIO.read(pathname, 'abi')
        
        traces_ = record.annotations['abif_raw']
        traces  = {}
        for k, v in traces_.items():
            if isinstance(v, bytes):
                v = v.decode()
            else:
                v = v
            traces[k] = v
        self.traces = traces

        self.data = record.seq

        phred_scores = []
        for c in traces['PCON2']:
            phred_scores.append(ord(c))
        self.phred_scores = phred_scores
    
    def find_alignable_window(self, window_size=30, QUAL_CUTOFF=20):
        '''
        This code finds the optimal window for aligning the edited trace to the reference sequence
        1) It looks at the phred scores and calculates a moving average
        2) It then selects the largest window where the moving average is greater than QUAL_CUTOFF
        3) This segment and other data is then returned in a dictionary
        all coordinates are absolute with respect to the sanger object
        '''

        phreds = np.asarray(self.phred_scores)
        # find a running mean
        window_size = window_size
        windowed_phred = running_mean(phreds, window_size)

        if len(windowed_phred) is 0:
            return {}
        
        # find regions with windowed mean > QUAL_CUTOFF
        y = find_regions_greater_than(windowed_phred, QUAL_CUTOFF)
        

        max_window_size = 0
        max_window = None
        for region in y:
            w_size = region[1] - region[0]
            if w_size > max_window_size:
                max_window_size = w_size
                max_window = (region[0], region[1])

        return {
            'regions': y,
            'max_window': max_window,
            'windowed_phred': windowed_phred
        }




def find_regions_greater_than(array, comparison):
    condition = array > comparison
    return contiguous_regions(condition)

 
def align(self, refSeq):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(refSeq, self.seq)

def running_mean(x, n):
    """
    https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    """

    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[n:] - cumsum[:-n]) / n

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.

    # From
    # https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array

    """

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx