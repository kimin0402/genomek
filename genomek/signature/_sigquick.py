import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import permutations
from ._context import revcompl



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

# context_192 = []
# for i in list2:
#     for mut1 in [f"[{i}>{x}]" for x in list2 if x != i]:
#         for j in list2:
#             mut2 = j+mut1
#             for k in list2:
#                 mut3 = mut2+k
#                 context_192.append(mut3)

context_192 = [f"{pre}[{ref}>{alt}]{post}" for ref, alt in permutations("ACGT",2) for pre in 'ACGT' for post in 'ACGT']


color_dict = {0: '#03BCEE',
             1:'#010101',
             2:'#E32926',
             3:'#CAC9C9',
             4:'#A1CE63',
             5:'#EBC5C3'}


norm_12 = np.array([0.04363834, 0.10423708, 0.04363834, 0.10423708, 0.04363834,
       0.10423708, 0.10309029, 0.08236762, 0.10309029, 0.08236762,
       0.10309029, 0.08236762])
norm_96 = np.array([0.01102527, 0.01197087, 0.00396346, 0.01152825, 0.01132706,
       0.0140029 , 0.00444632, 0.01440528, 0.0066393 , 0.00849026,
       0.00420489, 0.00927491, 0.01217206, 0.00971753, 0.00488894,
       0.00975777, 0.01102527, 0.01197087, 0.00396346, 0.01152825,
       0.01132706, 0.0140029 , 0.00444632, 0.01440528, 0.0066393 ,
       0.00849026, 0.00420489, 0.00927491, 0.01217206, 0.00971753,
       0.00488894, 0.00975777, 0.01102527, 0.01197087, 0.00396346,
       0.01152825, 0.01132706, 0.0140029 , 0.00444632, 0.01440528,
       0.0066393 , 0.00849026, 0.00420489, 0.00927491, 0.01217206,
       0.00971753, 0.00488894, 0.00975777, 0.0139023 , 0.00975777,
       0.01162884, 0.01420409, 0.01571302, 0.01102527, 0.00762514,
       0.01060277, 0.01068324, 0.00553275, 0.01024062, 0.01205134,
       0.0149485 , 0.01024062, 0.0116892 , 0.01559231, 0.0139023 ,
       0.00975777, 0.01162884, 0.01420409, 0.01571302, 0.01102527,
       0.00762514, 0.01060277, 0.01068324, 0.00553275, 0.01024062,
       0.01205134, 0.0149485 , 0.01024062, 0.0116892 , 0.01559231,
       0.0139023 , 0.00975777, 0.01162884, 0.01420409, 0.01571302,
       0.01102527, 0.00762514, 0.01060277, 0.01068324, 0.00553275,
       0.01024062, 0.01205134, 0.0149485 , 0.01024062, 0.0116892 ,
       0.01559231])
norm_192 = np.array([0.00201191, 0.00901336, 0.00160953, 0.01036134, 0.00156929,
       0.00239417, 0.00323918, 0.00828907, 0.00199179, 0.00933526,
       0.00144858, 0.01255432, 0.00160953, 0.00283679, 0.00350072,
       0.01090455, 0.00247465, 0.00416465, 0.00303798, 0.00545228,
       0.00311846, 0.00108643, 0.00567359, 0.00360132, 0.00382263,
       0.00834943, 0.00245453, 0.007263  , 0.00245453, 0.00243441,
       0.0035812 , 0.00617657, 0.00201191, 0.00901336, 0.00160953,
       0.01036134, 0.00156929, 0.00239417, 0.00323918, 0.00828907,
       0.00199179, 0.00933526, 0.00144858, 0.01255432, 0.00160953,
       0.00283679, 0.00350072, 0.01090455, 0.00247465, 0.00416465,
       0.00303798, 0.00545228, 0.00311846, 0.00108643, 0.00567359,
       0.00360132, 0.00382263, 0.00834943, 0.00245453, 0.007263  ,
       0.00245453, 0.00243441, 0.0035812 , 0.00617657, 0.00201191,
       0.00901336, 0.00160953, 0.01036134, 0.00156929, 0.00239417,
       0.00323918, 0.00828907, 0.00199179, 0.00933526, 0.00144858,
       0.01255432, 0.00160953, 0.00283679, 0.00350072, 0.01090455,
       0.00247465, 0.00416465, 0.00303798, 0.00545228, 0.00311846,
       0.00108643, 0.00567359, 0.00360132, 0.00382263, 0.00834943,
       0.00245453, 0.007263  , 0.00245453, 0.00243441, 0.0035812 ,
       0.00617657, 0.00651859, 0.00738371, 0.00229358, 0.00746419,
       0.00836955, 0.0032593 , 0.00756478, 0.0066393 , 0.00519073,
       0.01052229, 0.00259536, 0.00842991, 0.0040037 , 0.00362144,
       0.00420489, 0.00639788, 0.0075849 , 0.00309834, 0.00340013,
       0.00213263, 0.00913407, 0.00110655, 0.00995896, 0.00209239,
       0.00832931, 0.00661919, 0.00404394, 0.00619668, 0.00935538,
       0.00233382, 0.01054241, 0.0050499 , 0.00651859, 0.00738371,
       0.00229358, 0.00746419, 0.00836955, 0.0032593 , 0.00756478,
       0.0066393 , 0.00519073, 0.01052229, 0.00259536, 0.00842991,
       0.0040037 , 0.00362144, 0.00420489, 0.00639788, 0.0075849 ,
       0.00309834, 0.00340013, 0.00213263, 0.00913407, 0.00110655,
       0.00995896, 0.00209239, 0.00832931, 0.00661919, 0.00404394,
       0.00619668, 0.00935538, 0.00233382, 0.01054241, 0.0050499 ,
       0.00651859, 0.00738371, 0.00229358, 0.00746419, 0.00836955,
       0.0032593 , 0.00756478, 0.0066393 , 0.00519073, 0.01052229,
       0.00259536, 0.00842991, 0.0040037 , 0.00362144, 0.00420489,
       0.00639788, 0.0075849 , 0.00309834, 0.00340013, 0.00213263,
       0.00913407, 0.00110655, 0.00995896, 0.00209239, 0.00832931,
       0.00661919, 0.00404394, 0.00619668, 0.00935538, 0.00233382,
       0.01054241, 0.0050499 ])


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




def sig_plot_12(df:pd.DataFrame, title:str, normalize_MT=True):
    mtx = pd.crosstab(df['substitution'], df['strand']).reindex(index=list1, columns=['Heavy', 'Light']).fillna(0).reset_index().rename(columns={'substitution':'context'})
    mtx = pd.melt(mtx, id_vars='context', var_name='strand', value_name='count')
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(list1, inplace=True, ordered=True)
    mtx = mtx.sort_values(['context', 'strand']).reset_index().drop(columns='index')
    if normalize_MT:
        S = mtx['count'].sum()
        mtx['count'] = mtx['count'] / S / norm_12
    
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
    # fig.show()


def sig_plot_96(df:pd.DataFrame, title:str, normalize_MT=True):
    mtx = df['context'].value_counts()
    mtx = pd.DataFrame(index=context_96).join(mtx).fillna(0).reset_index().rename(columns={'context':'count', 'index':'context'})
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(context_96, inplace=True, ordered=True)
    if normalize_MT:
        S = mtx['count'].sum()
        mtx['count'] = mtx['count'] / S / norm_96
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
    # fig.show()


def sig_plot_192(df:pd.DataFrame, title:str, normalize_MT=True):
    mtx = pd.crosstab(df['context'], df['strand']).reindex(index=context_96, columns=['Heavy', 'Light']).fillna(0).reset_index().rename(columns={'substitution':'context'})
    mtx = pd.melt(mtx, id_vars='context', var_name='strand', value_name='count')
    mtx['context'] = mtx['context'].astype('category')
    mtx['context'].cat.set_categories(context_96, inplace=True, ordered=True)
    mtx = mtx.sort_values(['context', 'strand']).reset_index().drop(columns='index')
    if normalize_MT:
        S = mtx['count'].sum()
        mtx['count'] = mtx['count'] / S / norm_192
    
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
    # fig.show()



def sig_plot_64(df, title):
    color_dict = {0: "#2E8534",
                  1: "#181DB9",
                  2: "#CC8245",
                  3: "#D01027"}
    base = ['A', 'C', 'G', 'T']

    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16,4), sharey=True)
    for i, ax in enumerate(axes):
        _df = df.loc[(i)*16:(i+1)*16-1]
        sns.barplot(data=_df, x='context', y='count', color=color_dict[i], ax=ax)
        if i != 0: ax.yaxis.set_visible(False)
        ax.set_xticklabels(_df['context'].unique(), rotation=90)
        ax.set_xlabel(base[i], fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)
    plt.tight_layout(w_pad=0.5)
    plt.suptitle(title, fontsize=24)
    fig.show()


def sig_plot_32(df, title):
    color_dict = {0: "#181DB9",
                  1: "#D01027"}
    base = ['C', 'T']

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,6), sharey=True)
    for i, ax in enumerate(axes):
        _df = df.loc[(i)*16:(i+1)*16-1]
        sns.barplot(data=_df, x='context', y='count', color=color_dict[i], ax=ax)
        if i != 0: ax.yaxis.set_visible(False)
        ax.set_xticklabels(_df['context'].unique(), rotation=90)
        ax.set_xlabel(base[i], fontsize=16)
        for spine in ax.spines.values():
            spine.set_visible(False)
    plt.tight_layout(w_pad=1.5)
    plt.suptitle(title, fontsize=24)
    fig.show()











# MT_6 = Counter(fasta.fetch('MT'))
# norm_6 = np.array([[MT_6['G'], MT_6['C']], 
#           [MT_6['G'], MT_6['C']],
#           [MT_6['G'], MT_6['C']],
#           [MT_6['A'], MT_6['T']],
#           [MT_6['A'], MT_6['T']],
#           [MT_6['A'], MT_6['T']]]) / len(fasta.fetch('MT')) / 3

# MT_4 = {'G': 2169, 'A': 5124, 'T': 4094, 'C': 5181, 'N': 1}
# norm_12 = np.array([MT_4['G'], MT_4['C'], MT_4['G'], MT_4['C'], MT_4['G'], MT_4['C'],
#                    MT_4['A'], MT_4['T'], MT_4['A'], MT_4['T'], MT_4['A'], MT_4['T']]) / 16568 / 3
# reference_path = "/home/users/kimin/projects/12_MT/00_reference/MT_hg19.fasta"
# fasta = pysam.FastaFile(reference_path)
# seq = fasta.fetch('MT')
# mt_64_context_dict = dict(sorted(Counter([x for x in map("".join, zip(seq, seq[1:], seq[2:])) if not 'N' in x]).items(), key=lambda x: (x[0][1], x[0][0], x[0][2])))
# mt_32_context_dict = defaultdict(int)
# for context, count in mt_64_context_dict.items():
#     if context[1] == "A" or context[1] == "G":
#         mt_32_context_dict[revcompl(context)] += count
#     else:
#         mt_32_context_dict[context] += count
# mt_32_context_dict = dict(sorted(mt_32_context_dict.items(), key=lambda x: (x[0][1], x[0][0], x[0][2])))
# norm_32 = list(mt_32_context_dict.values())

# result = {}
# for i, (context, count) in enumerate(mt_64_context_dict.items()):
#     if context[1] == 'A' or context[1] == 'G':
#         strand = 'Heavy'
#         context_adj = revcompl(context)
#     else:
#         strand = 'Light'
#         context_adj = context
#     result[i] = {'context': context_adj, 'strand': strand, 'count': count}
    
# norm_64 = pd.DataFrame.from_dict(dict(sorted(result.items(), key=lambda x: (x[1]['context'][1], x[1]['context'][0], x[1]['context'][2], x[1]['strand']))), orient='index')['count'].to_list()
# norm_192 = np.array(norm_64[0:32] + norm_64[0:32] + norm_64[0:32] + norm_64[32:64] + norm_64[32:64] + norm_64[32:64]) / 16568 / 3