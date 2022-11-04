import pyranges as pr
from math import factorial as f
from math import log
import numpy as np


def log10(x):
    if x != 0:
        return log(x, 10)
    else:
        return log(0.000001, 10)

def nCr(n,r):
    return log10(f(n))-log10(f(r))-log10((f(n-r)))

def df_to_p(df, chrom_key='CHROM', pos_key='POS', ref_key='REF'):
    p = pr.PyRanges(chromosomes=df[chrom_key], starts=df[pos_key]-1, ends=df[pos_key]-1+df[ref_key].map(len))
    return p


def seg_to_p(seg, female=True):
    '''
    seg: pandas dataframe
    '''
    if not 'Y' in seg.iloc[:, 0].values:
        if female:
            seg.loc[len(seg), seg.columns] = 'Y', 0, 59373566, 0, 0, 0
        else:
            seg.loc[len(seg), seg.columns] = 'Y', 0, 59373566, 1, 1, 0
    p = pr.PyRanges(seg.set_axis(['Chromosome', 'Start', 'End'] + seg.columns[3:].to_list(), axis=1))
    return p


def df_annotate_cnt(df, bed, female=True, nb_cpu=1):
    '''
    df: pandas dataframe 
    bed: pandas dataframe
    return: numpy array of cnt, a, b, and distance
    ** nb_cpu raises error. pyranges problem.
    '''
    p = df_to_p(df)
    bed = seg_to_p(bed, female=female)
    result = p.nearest(bed, suffix='segment', nb_cpu=nb_cpu)
    cnt, a, b, d = getattr(result, 'CNt').to_numpy(int), getattr(result, 'A').to_numpy(int), getattr(result, 'B').to_numpy(int), getattr(result, 'Distance').to_numpy()
    return cnt, a, b, d


def get_ccf(purity, n_read, v_read, t_cn):
    if t_cn == 0 or n_read == 0 or v_read == 0:
        return 0, 0
    n_cn = 2
    t_fraction = t_cn * purity/(t_cn*purity + n_cn*(1-purity))
    max_binomP = 0
    max_cn = 1
    for m_cn in range(1, max(t_cn, 1)+1):
        current_binomP = sum(pow(10, nCr(n_read, i) + i*log10(t_fraction) + (n_read-i)*log10(1-t_fraction) + nCr(i, v_read) + v_read*log10((m_cn/t_cn)) + (i-v_read)*log10((1-m_cn/t_cn))) for i in range(v_read, n_read+1))
        if current_binomP >= max_binomP:
            max_binomP = current_binomP
            max_cn = m_cn
    ccf = v_read/n_read * t_cn/t_fraction/max_cn
    return ccf, max_cn


def variant_in_bed(df, bed):
    '''
    df: pandas dataframe 
    bed: pyranges pr object
    '''
    P_df = pr.PyRanges(chromosomes=df['CHROM'], starts=df['POS']-1, ends=df['POS']-1+df['REF'].map(len)).overlap(bed, how='first').as_df()
    if len(P_df):
        P_df = P_df.drop('End', axis=1).set_axis(['CHROM', 'POS'], axis=1).drop_duplicates(['CHROM', 'POS'])
        P_df['POS'] = P_df['POS'] + 1
        P_df['mask'] = True
    else:
        df['mask'] = False
        return df['mask'].to_numpy()
    if 'mask' in df.columns:
        df.drop(columns='mask', inplace=True)
    df = df.merge(P_df, on=['CHROM', 'POS'], how='left')
    df['mask'] = df['mask'].fillna(False)
    return df['mask'].to_numpy()