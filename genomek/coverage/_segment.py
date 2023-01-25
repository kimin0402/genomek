import re
import numpy as np

def df_filter_coord(df, coords, chrom_key='CHROM', pos_key='POS', return_drop=True):
    '''
    input: pandas df
    output: pandas df
    '''
    result_bool = np.repeat(False, len(df))
    for coord in coords:
        chrom, pos1, pos2 = re.findall(r'([\dXY]+):(\d+)-(\d+)', coord)[0]
        pos1, pos2 = int(pos1), int(pos2)
        filter_bool = (df[chrom_key] == chrom) & (df[pos_key] >= pos1) & (df[pos_key] <= pos2)
        result_bool = result_bool | filter_bool
    if return_drop:
        return df[result_bool]
    else:
        return df[~result_bool]
