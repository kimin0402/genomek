import sys

from . import tools as tl
from . import vcf as vcf
from . import signature as sig
from . import coverage as cov
from . import annotation as an
from . import visualization as viz
from . import bam as bam
from . import botseq as bot
from . import sv as sv
from . import qc as qc
from . import snv as snv

__author__ = 'Kyoung Il Min'


def hello():
    print("scripts folder is loaded!")
  
if __name__ == 'scripts':
    hello()


