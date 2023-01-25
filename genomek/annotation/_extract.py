import pandas as pd

def get_df_part(db_dir, group, format_key, chrom, pos, ref, alt):
    db_path = f"{db_dir}/{sample_name}.{format_key}.{chrom}.anraw.ftr"
    df_part = pd.read_feather(db_path).query(f"POS == {pos} and REF == {ref} and ALT == {alt}")
    return df_part

def df_part_to_f1r2_ratio(df_part):
    index = pd.MultiIndex.from_product([[True, False], [True, False]], names=['reverse_flag', 'read1_flag'])
    result = np.nan
    if (filter_bool := df_part['type'] == 1).sum() > 0:
        if (strand_filter_bool := filter_bool & ~df_part['sv_flag']).sum() > 0:
            strand_bias_series = df_part[strand_filter_bool][['reverse_flag', 'read1_flag']].value_counts().reindex(index).fillna(0)
            result = (strand_bias_series[(True, False)] + strand_bias_series[(False, True)])/strand_bias_series.sum()
    return result