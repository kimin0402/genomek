import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def kataegis(df, ax, hue_column, y_column='logDIFF'):
    '''
    fig, ax = scripts.viz.genome_figure()
    fig.suptitle("Title")
    for i, chrom in enumerate(chromosomes[:-1]):
        df_tmp = df.query("CHROM == @chrom")[['CHROM', 'POS', 'AF', 'substitution']]
        scripts.viz.kataegis(df_tmp, ax[chrom], y_column='AF', hue_column='substitution')
        if i != 0:
            try: 
                ax[chrom].get_legend().remove()
            except:
                pass
    '''

    df = df.sort_values('POS',axis=0, ascending=True)

    df['DIFF'] = df['POS'].diff()
    df['logDIFF'] = np.log10(df['DIFF'])

    hue_order  = ['C>T', 'C>A', 'C>G', 'T>C', 'T>A', 'T>G', 'others'] #,'G>A','G>C','G>T','T>A','T>C','T>G', 'not_snp']
    hue_dict = {'C>A': [72/256,69/256,156/256], 'C>G': [1/256,1/256,1/256], 'C>T': [228/256,41/256,38/256], 
                'T>A': [168/256,53/256,144/256], 'T>C': [232/256,229/256,56/256], 'T>G': [110/256,173/256,43/256], 'others': 'grey'}
    sns.scatterplot(x='POS', #range(len(vcf_ch)), 
                    y=y_column, 
                    hue=hue_column,
                    hue_order=hue_order,
                    palette=hue_dict, #'tab20_r', #gnuplot_r', #tab20',
                    ax=ax,
                    data=df)
