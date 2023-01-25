import re
import numpy as np
import pandas as pd
from ..tools._gadgets import chrom_sort_dict_37 as chrom_sort_dict
from ..tools._gadgets import chrom_cat_type_37 as chrom_cat_type
from ._sv import get_svtype
basic_columns_sv = ['CHR1', 'POS1', 'CHR2', 'POS2', 'SVTYPE', 'CT']

def row_func_sort_chr1_pos1_chr2_pos2(row):
    if row['CHR1'] == row['CHR2']:
        if row['POS1'] <= row['POS2']:
            return row['CHR1'], row['POS1'], row['CHR2'], row['POS2'], row['CT']
        else:
            return row['CHR2'], row['POS2'], row['CHR1'], row['POS1'], row['CT'][::-1].replace('ot', 'to')
    elif chrom_sort_dict[row['CHR1']] > chrom_sort_dict[row['CHR2']]:
        return row['CHR2'], row['POS2'], row['CHR1'], row['POS1'], row['CT'][::-1].replace('ot', 'to')
    else:
        return row['CHR1'], row['POS1'], row['CHR2'], row['POS2'], row['CT']

def process_delly(df):
    subset_bool = df['CHR2'].isna() | df['POS2'].isna()
    df.loc[subset_bool, 'CHR2'] = df.loc[subset_bool, 'CHROM']
    df.loc[subset_bool, 'POS2'] = df.loc[subset_bool, 'END']
    df = df.rename({'CHROM':'CHR1', 'POS':'POS1'}, axis=1)
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int, 'QUAL':float}).dropna(subset=['CHR1', 'CHR2'], axis=0)
    df[['CHR1', 'POS1', 'CHR2', 'POS2', 'CT']] = df.apply(row_func_sort_chr1_pos1_chr2_pos2, axis=1, result_type='expand')
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int}).dropna(subset=['CHR1', 'CHR2'], axis=0)
    df['SVTYPE'].replace('BND', 'TRA', inplace=True)
    df = df.sort_values(['CHR1', 'CHR2', 'POS1', 'POS2'])
    df = df[basic_columns_sv + [x for x in df.columns if x not in basic_columns_sv]]
    df.set_axis(basic_columns_sv + [f"DELLY_{x}" for x in df.columns if x not in basic_columns_sv], axis=1, inplace=True)
    return df

def row_func_ct(row):
    if row[0]:
        low_end_ct = '3'
    elif row[5]:
        low_end_ct = '5'
    else:
        raise ValueError('row parsing error')
    assert row[1] == row[4]
    if row[1] == ']':
        high_end_ct = '3'
    else:
        high_end_ct = '5'
    return low_end_ct + "to" + high_end_ct

def process_svaba(df, prefix):
    original_colnames = df.filter(like=':::').columns.tolist()
    if '01_bam/' in original_colnames[0]:
        fixed_colnames = [re.sub("/.*", "", re.sub("(:::).*01_bam/", r"\1", x)) for x in original_colnames]
        rename_dict = {k:v for k,v in zip(original_colnames, fixed_colnames)}
        df.rename(rename_dict, axis=1, inplace=True)
    df_parsed = df['ALT'].str.extract(r"([ACTGN]*)([\[\]])(.*):(\d+)([\[\]])([ACTGN]*)")
    df['CHR2'] = df_parsed[2]
    df['POS2'] = df_parsed[3]
    df = df.rename({'CHROM':'CHR1', 'POS':'POS1'}, axis=1)
    df['CT'] = df_parsed.apply(row_func_ct, axis=1)
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int, 'QUAL':float}).dropna(subset=['CHR1', 'CHR2'], axis=0)
    df[['CHR1', 'POS1', 'CHR2', 'POS2', 'CT']] = df.apply(row_func_sort_chr1_pos1_chr2_pos2, axis=1, result_type='expand')
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int}).dropna(subset=['CHR1', 'CHR2'])
    df = df.sort_values(['CHR1', 'CHR2', 'POS1', 'POS2'])
    df['SVTYPE'] = df.apply(lambda x: get_svtype(x['CHR1'], x['CHR2'], x['CT']), axis=1)
    df = df[basic_columns_sv + [x for x in df.columns if x not in basic_columns_sv]]
    df.set_axis(basic_columns_sv + [f"{prefix}_{x}" for x in df.columns if x not in basic_columns_sv], axis=1, inplace=True)
    dedup_bool = np.all(df[['CHR2', 'POS2']].astype(str).to_numpy() == df[f'{prefix}_ALT'].str.extract(r"[\[\]](.*):(\d+)[\[\]]").to_numpy(), axis=1)
    df = df[dedup_bool].reset_index(drop=True)
    return df

def process_smoove_non_bnd(df, prefix):
    df['CHR2'] = df['CHROM']
    df['POS2'] = df['END']
    df = df.rename({'CHROM':'CHR1', 'POS':'POS1'}, axis=1)
    ct = df['STRANDS'].str.split(":").str[0].str.replace('-','5').str.replace('+','3')
    df['CT'] = ct.str[0] + "to" + ct.str[1]
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int, 'QUAL':float}).dropna(subset=['CHR1', 'CHR2'], axis=0)
    df[['CHR1', 'POS1', 'CHR2', 'POS2', 'CT']] = df.apply(row_func_sort_chr1_pos1_chr2_pos2, axis=1, result_type='expand')
    df = df.astype({'CHR1':chrom_cat_type, 'POS1':int, 'CHR2':chrom_cat_type, 'POS2':int}).dropna(subset=['CHR1', 'CHR2'], axis=0)
    df = df.sort_values(['CHR1', 'CHR2', 'POS1', 'POS2'])
    df = df[basic_columns_sv + [x for x in df.columns if x not in basic_columns_sv]]
    df.set_axis(basic_columns_sv + [f"{prefix}_{x}" for x in df.columns if x not in basic_columns_sv], axis=1, inplace=True)
    return df

def process_smoove(df, prefix):
    df = pd.concat([process_smoove_non_bnd(df.query("SVTYPE != 'BND'"), prefix='SMOOVE'),
                    process_svaba(df.query("SVTYPE == 'BND'"), prefix='SMOOVE')], axis=0).sort_values(['CHR1', 'CHR2', 'POS1', 'POS2']).reset_index(drop=True)
    return df


def merge_svs(df, padding=500):
    df = df.reset_index(drop=True)
    fillna_bool = df.filter(regex="^(DELLY|SVABA|SMOOVE)_([TXO]{1}|CALL)").columns
    for i in range(len(df)-1):
        # if (df.loc[i, ['DELLY_CALL', 'SVABA_CALL', 'SMOOVE_CALL']] == df.loc[i+1, ['DELLY_CALL', 'SVABA_CALL', 'SMOOVE_CALL']]).all():
        #     continue
        if df.loc[i, 'CHR1'] == df.loc[i+1, 'CHR1'] and \
           df.loc[i, 'CHR2'] == df.loc[i+1, 'CHR2'] and \
           df.loc[i, 'SVTYPE'] == df.loc[i+1, 'SVTYPE'] and \
           df.loc[i, 'CT'] == df.loc[i+1, 'CT'] and \
           abs(df.loc[i, 'POS1'] - df.loc[i+1, 'POS1']) < padding and \
           abs(df.loc[i, 'POS2'] - df.loc[i+1, 'POS2']) < padding:
            df.loc[i+1, fillna_bool] = df.loc[i:i+1, fillna_bool].any(axis=0)
            df.loc[i, fillna_bool] = False    
    return df
