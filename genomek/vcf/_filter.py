from pandas.core.frame import DataFrame
import pyranges as pr

class df_filter(DataFrame):
    @classmethod
    def new(cls, df, original_keys, edit_keys, new_keys):
        assert len(original_keys) == len(edit_keys)
        df.__class__ = df_filter
        df.edit = df.copy()
        df.edit[edit_keys] = df[original_keys]
        for k in new_keys:
            df.edit[k] = None
        return df

    def update_value(self, i, edit_key):
        updated_value = input(f"{edit_key}: ")
        self.edit.loc[i, edit_key] = updated_value
            
    def mark_true(self, i, edit_key):
        self.edit.loc[i, edit_key] = True


def variant_in_bed(df, bed, chrom_key='CHROM', pos_key='POS', ref_key='REF'):
    '''
    df: pandas dataframe 
    bed: pyranges pr object
    '''
    P_df = pr.PyRanges(chromosomes=df[chrom_key], starts=df[pos_key]-1, ends=df[pos_key]-1+df[ref_key].map(len)).overlap(bed, how='first').as_df()
    if len(P_df):
        P_df = P_df.drop('End', axis=1).set_axis([chrom_key, pos_key], axis=1).drop_duplicates([chrom_key, pos_key])
        P_df[pos_key] = P_df[pos_key] + 1
        P_df['mask'] = True
    else:
        df['mask'] = False
        return df['mask'].to_numpy()
    if 'mask' in df.columns:
        df.drop(columns='mask', inplace=True)
    df = df.merge(P_df, on=[chrom_key, pos_key], how='left')
    df['mask'] = df['mask'].fillna(False)
    return df['mask'].to_numpy()