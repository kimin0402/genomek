import os
import tempfile
import subprocess
import pandas as pd
from ..tools._gadgets import chrom_cat_type_37 as chrom_cat_type

mosdepth_path = "/home/users/kimin/bin/mosdepth"

def mito_gc_cov(bam_path: str, bed_path: str) -> pd.DataFrame:

    with tempfile.TemporaryDirectory() as tempdir:
        output_prefix = f"{tempdir}/mito_gc"
        output_path = f"{tempdir}/mito_gc.regions.bed.gz"
        cmd = f"{mosdepth_path} -n -t 4 -c MT -b {bed_path} {output_prefix} {bam_path}"
        subprocess.run(cmd, shell=True)

        df = pd.read_csv(output_path, compression='gzip', names=['CHROM', 'start', 'end', 'GC', 'cov'], index_col=None, sep='\t') #, quotechar='"', error_bad_lines=False)

    return df


def mosdepth_summary_to_info(summary_path):
    if os.path.exists(summary_path) and os.path.getsize(summary_path) > 0:
        mosdepth_summary = pd.read_csv(summary_path, header=0, index_col=0, sep='\t')
        total_length, total_bases = mosdepth_summary.filter(like='region', axis=0).iloc[:-2,0:2].sum(axis=0)
        avg = total_bases/total_length
        try:
            mt_depth = float(mosdepth_summary.loc['MT', 'mean'])
        except:
            mt_depth = None
        return avg, mt_depth
    else:
        return None, None


def bam_to_cov(bam_path: str, bed_path: str, quality: int=30, threads: int=4) -> pd.DataFrame:
    '''
    input bam_path, bed_path and returns dataframe of mosdepth output
    Mosdepth Options:
      -t --threads <threads>     number of BAM decompression threads. (use 4 or fewer) [default: 0]
      -c --chrom <chrom>         chromosome to restrict depth calculation.
      -b --by <bed|window>       optional BED file or (integer) window-sizes.
      -n --no-per-base           dont output per-base depth. skipping this output will speed execution
                                 substantially. prefer quantized or thresholded values if possible.
      -f --fasta <fasta>         fasta file for use with CRAM files.
      --d4                       output per-base depth in d4 format. This is much faster.
    '''

    with tempfile.TemporaryDirectory() as tempdir:
        output_prefix = f"{tempdir}/mosdepth_temp"
        output_path = f"{tempdir}/mosdepth_temp.regions.bed.gz"
        summary_path = f"{tempdir}/mosdepth_temp.mosdepth.summary.txt"
        cmd = f"{mosdepth_path} -n -t {threads} -Q {quality} -b {bed_path} {output_prefix} {bam_path}"
        subprocess.run(cmd, shell=True)
        df = pd.read_csv(output_path, sep='\t', compression='gzip', header=None, index_col=None, names=['Chromosome', 'Start', 'End', 'GC', 'cov'],
                     dtype={'Chromosome': chrom_cat_type, 'Start': int, 'End': int, 'cov': float}) #, quotechar='"', error_bad_lines=False)
        avg_depth_mosdepth, MT_depth_mosdepth = mosdepth_summary_to_info(summary_path)

    return df, avg_depth_mosdepth, MT_depth_mosdepth



def get_foldback_n(bam, chrom, start, end):
    '''
    calculate number of foldback reads in the given segment
    '''
    n = 0
    for read in bam.fetch(contig=chrom, start=start, end=end):
        if read.is_duplicate or read.is_qcfail or read.mate_is_unmapped or read.is_unmapped: continue
        if read.is_reverse == read.mate_is_reverse and read.reference_name == read.next_reference_name and read.template_length <= 5000:
            n += 1
    return n