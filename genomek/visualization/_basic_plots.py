import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def genome_figure(chrom="all", figsize=(26.5, 6.5)):
    '''
    chrom: list of strings indicating which chromosome to return

    Usage:
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
    fig = plt.figure(figsize=figsize)
    ax = {}

    df_fai = pd.read_csv("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai", sep='\t', header=None, index_col=0, usecols=[0,1], nrows=24, names=['CHROM', 'SIZE'] )
    df_fai = df_fai / df_fai.sum(axis=0)
    df_fai['SIZE_SUM'] = pd.concat([pd.Series([0]), df_fai['SIZE'].cumsum()])[:-1].to_numpy()
    df_fai

    if chrom != "all":
        df_fai = df_fai.loc[chrom]

    for i, CHROM in enumerate(df_fai.index):
        start = df_fai['SIZE_SUM'][i]
        end = df_fai['SIZE'][i]
        if i == 0:
            ax[CHROM] = fig.add_axes([start, 0, end, 0.9])
        else:
            ax[CHROM] = fig.add_axes([start, 0, end, 0.9], sharey=ax['1'])
            ax[CHROM].axes.yaxis.set_visible(False)
        ax[CHROM].set_title(CHROM)

    return fig, ax 


'''
from matplotlib_venn import venn2

venn2([set(adata_GCSF.var.index[adata_GCSF.var['highly_variable_vst']]),
  set(adata_GCSF.var.index[adata_GCSF.var['highly_variable_scran']])],
  ['vst', 'scran'], 
  ax=ax)

from matplotlib_venn import venn3

venn3([set(df.index[df['highly_variable_vst']]),
      set(df.index[df['highly_variable_seurat']]),
      set(df.index[df['highly_variable_cellranger']])], ['vst', 'seurat', 'cellranger'])

from pyvenn.venn import venn5
from pyvenn.venn import venn4
from pyvenn.venn import get_labels
labels = get_labels([set(df.index[df['highly_variable_vst']]),
      set(df.index[df['highly_variable_seurat']]),
      set(df.index[df['highly_variable_cellranger']]),
      set(df.index[df['highly_variable_scran']]),
      set(df.index[df['highly_variable_poisson']])])
venn5(labels, names=['vst', 'seurat', 'cellranger', 'scran', 'poisson'])
plt.show()
labels = get_labels([set(df.index[df['highly_variable_vst']]),
      set(df.index[df['highly_variable_seurat']]),
      set(df.index[df['highly_variable_cellranger']]),
      set(df.index[df['highly_variable_scran']])])
venn4(labels, names=['vst', 'seurat', 'cellranger', 'scran'])
plt.show()
'''