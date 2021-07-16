import os
import sys
import io
import gzip
import tempfile
import subprocess

import re
import pysam
import cyvcf2
import pandas as pd
from pandas.api.types import CategoricalDtype

mosdepth_path = "/home/users/kimin/bin/mosdepth"

def mito_gc_cov(bam_path: str, bed_path: str) -> pd.DataFrame:

    with tempfile.TemporaryDirectory() as tempdir:
        output_prefix = f"{tempdir}/mito_gc"
        output_path = f"{tempdir}/mito_gc.regions.bed.gz"
        cmd = f"{mosdepth_path} -n -t 4 -c MT -b {bed_path} {output_prefix} {bam_path}"
        subprocess.run(cmd, shell=True)

        df = pd.read_csv(output_path, compression='gzip', names=['CHROM', 'start', 'end', 'GC', 'cov'], index_col=None, sep='\t') #, quotechar='"', error_bad_lines=False)

    return df

