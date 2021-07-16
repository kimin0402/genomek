import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


list1 = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
list2 = ['A', 'C', 'G', 'T']
context_96 = []
for i in list1:
    i1 = "[" + i + "]"
    for j in list2:
        i2 = j + i1
        for k in list2:
            i3 = i2 + k
            context_96.append(i3)

context_96_for_sigprofiler = []
for i in list2:
    str1 = i + "["
    for j in list1:
        str2 = str1 + j + "]"
        for k in list2:
            str3 = str2 + k
            context_96_for_sigprofiler.append(str3)


def data_generation_96(df: pd.DataFrame, column_name: str, sample_id: str) -> pd.DataFrame:
    '''
    Generate sigprofiler input

    df: df object constructed from vcf
    column_name: Name of the column containing context information. 
    sample_id: String used for the column name of the 96 data matrix.
    '''

    df_96 = df[column_name].value_counts()
    df_96 = pd.DataFrame(index=context_96_for_sigprofiler).join(df_96).fillna(0).rename(columns={column_name:sample_id})
    df_96.index.rename("MutationType", inplace=True)
    
    return df_96




def sig_plot_12(df, title):
    mtx = pd.crosstab(df['substitution'], df['strand'])
    mtx = pd.DataFrame(index=list1).join(mtx).fillna(0).reset_index().rename(columns={'index':'context'})
    mtx = pd.melt(mtx, id_vars='context', var_name='strand', value_name='count')
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(list1, inplace=True, ordered=True)
    mtx = mtx.sort_values(['context', 'strand']).reset_index().drop(columns='index')
    
    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(24,6), sharey=True)
    for i, ax in enumerate(axes):
        temp = mtx.loc[i*2:(i+1)*2 - 1].astype({'context':'str'})
        sns.barplot(data=temp, x='context', y='count', hue='strand', ax=ax)
        if i != 0: ax.yaxis.set_visible(False)
        if i != 5: ax.get_legend().remove()
        ax.set_xticklabels(temp['context'].unique(), rotation=90)
        ax.set_xlabel(list1[i], fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.tight_layout(w_pad=0.5)
    plt.suptitle(title, fontsize=24)
    fig.show()


def sig_plot_96(df, title):
    mtx = df['context'].value_counts()
    mtx = pd.DataFrame(index=context_96).join(mtx).fillna(0).reset_index().rename(columns={'context':'count', 'index':'context'})
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(context_96, inplace=True, ordered=True)
    color_dict = {0: '#03BCEE',
                 1:'#010101',
                 2:'#E32926',
                 3:'#CAC9C9',
                 4:'#A1CE63',
                 5:'#EBC5C3'}

    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(24,6), sharey=True)
    for i, ax in enumerate(axes):
        temp = mtx.loc[i*16:(i+1)*16 - 1].astype({'context':'str'})
        sns.barplot(data=temp, x='context', y='count', color=color_dict[i], ax=ax)
        if i != 0: ax.yaxis.set_visible(False)
        ax.set_xticklabels(temp['context'].unique(), rotation=90)
        ax.set_xlabel(list1[i], fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.tight_layout(w_pad=0.5)
    plt.suptitle(title, fontsize=24)
    fig.show()


def sig_plot_192(df, title):
    mtx = pd.crosstab(df['context'], df['strand'])
    mtx = pd.DataFrame(index=context_96).join(mtx).fillna(0).reset_index().rename(columns={'index':'context'})
    mtx = pd.melt(mtx, id_vars='context', var_name='strand', value_name='count')
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(context_96, inplace=True, ordered=True)
    mtx = mtx.sort_values(['context', 'strand']).reset_index().drop(columns='index')
    
    fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(24,6), sharey=True)
    for i, ax in enumerate(axes):
        temp = mtx.loc[i*32:(i+1)*32 - 1].astype({'context':'str'})
        sns.barplot(data=temp, x='context', y='count', hue='strand', ax=ax)
        if i != 0: ax.yaxis.set_visible(False)
        if i != 5: ax.get_legend().remove()
        ax.set_xticklabels(temp['context'].unique(), rotation=90)
        ax.set_xlabel(list1[i], fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.tight_layout(w_pad=0.5)
    plt.suptitle(title, fontsize=24)
    fig.show()
