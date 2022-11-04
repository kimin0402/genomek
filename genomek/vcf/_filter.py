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