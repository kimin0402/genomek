import os
import sys
import psutil
import itertools
from collections import defaultdict


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


