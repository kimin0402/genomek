import os
import sys
import psutil
import itertools
import pysam
from collections import defaultdict
from pandas.api.types import CategoricalDtype

chromosomes_37 = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]
chromosomes_38 = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y', 'M']]
chrom_cat_type_37 = CategoricalDtype(categories=chromosomes_37, ordered=True)
chrom_cat_type_38 = CategoricalDtype(categories=chromosomes_38, ordered=True)
chrom_sort_dict_37 = {k:v for k,v in zip(chromosomes_37, range(len(chromosomes_37)))}
satellite_path_37 = "/home/users/kimin/projects/00_Reference/repeatmasker/GRCh37.rmsk.satellite.bed"
mask_path_37 = "/home/users/kimin/projects/00_Reference/repeatmasker/GRCh37.kimin.mask.bed"
## mask_path_37 was created from this notebook -> /home/users/kimin/projects/07_Gastric_Cancer/20221008_pon_satellite_study.ipynb 
fasta_37 = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")



revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','[':'[',']':']','>':'>'}[B] for B in x][::-1])


def memory_usage():
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse1 = py.memory_info()[0]/2.**30  # memory use in GB...I think
    print('\n************** Reported Current Memory Use: '+ str(round(memoryUse1,2))+" GB *****************\n") 
    #print('\n************** Reported Current Memory Use: '+ str(round(memoryUse2,2))+" GB *****************\n")


def perm(n, seq):
    '''
    From Sig profiler
    
    Generates a list of all available permutations of n-mers.
    Parameters:
               n  -> length of the desired permutation string
             seq  -> list of all possible string valuesvcf_pathvcf_path
    Returns:
          permus  -> list of all available permutations
    '''
    permus = []
    for p in itertools.product(seq, repeat=n):
        permus.append("".join(p))
    return(permus)


def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, flush=True, **kwargs)
    return



def get_vaf(a, b):
    if a+b == 0:
        return None
    return a / (a+b) 
