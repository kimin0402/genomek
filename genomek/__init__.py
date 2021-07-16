import sys

from . import tools as tl
from . import vcf as vcf
from . import signature as sig
from . import snv as snv
from . import coverage as cov
from . import annotation as an
from . import visualization as viz
from . import bam as bam
from . import botseq as bot

__author__ = 'KyoungIl Min'


def hello():
    print("scripts folder is loaded!")
  
if __name__ == 'scripts':
    hello()


from pandas.api.types import CategoricalDtype
chromosomes_37 = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]
chromosomes_38 = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y', 'M']]
chrom_cat_type_37 = CategoricalDtype(categories=chromosomes_37, ordered=True)
chrom_cat_type_38 = CategoricalDtype(categories=chromosomes_38, ordered=True)