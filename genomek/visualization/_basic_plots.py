import sys
from typing import Union
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.patches as patches
from matplotlib.axes._axes import Axes
from matplotlib.path import Path
from ..tools import print_err
import pyranges as pr

base_dict = {'A':'#35ab35', 'C':'#0000ff', 'G':'#d7862a', 'T':'#ff0000'}
snv_hue_order  = ['C>T', 'C>A', 'C>G', 'T>C', 'T>A', 'T>G', 'others'] #,'G>A','G>C','G>T','T>A','T>C','T>G', 'not_snp']
snv_hue_dict = {'C>A': [72/256,69/256,156/256], 'C>G': [1/256,1/256,1/256], 'C>T': [228/256,41/256,38/256], 
            'T>A': [168/256,53/256,144/256], 'T>C': [232/256,229/256,56/256], 'T>G': [110/256,173/256,43/256], 'others': 'grey'}
sv_color_dict = {'DUP': 'green', 'TRA':'purple', 'DEL':'red', 'INV':'blue'}

class Axes_genome(Axes):
    @classmethod
    def new(cls, ax, title, chrom, pos, leading=False):
        ax.__class__ = Axes_genome
        ax.genome_title = title
        ax.chrom = chrom
        ax.start = pos[0]
        ax.end = pos[1]
        ax.pos_range = range(pos[0]+1, pos[1]+1) # int being tested upon the range will be 1-based coordinate (not 0-based)
        ax.set_xlim(pos[0], pos[1])
        ax.leading = leading
        if leading:
            pass
        return ax

    def __repr__(self):
        return f"Genome Axes object of coordinate: {self.chrom}:{self.start}-{self.end}"


def process_chrom_pos_coord(coord: Union[str,list]=None, chrom: Union[str,list]=None, pos: Union[str,list]=None, reference_fai=''):
    title_result = []
    chrom_result = []
    pos_result = []
    if not len(reference_fai):
        sys.exit("please specify reference_fai")
    if coord is not None:
        if isinstance(coord, str):
            coord = [coord]
        assert isinstance(coord, list), "coordinate input should be either str or list"
        for c in coord:
            title_result.append(c)
            c_split = c.split(":")
            chrom_result.append(c_split[0])
            if len(c_split) > 1:
                pos_result.append(tuple([int(x) for x in c_split[1].split("-")]))
            else:
                end_pos = pd.read_csv(reference_fai, sep='\t', header=None, index_col=0, usecols=[0,1], nrows=24, names=['CHROM', 'size_int'] ).loc[c_split[0]].squeeze()
                pos_result.append((0, end_pos))    
    if chrom is not None and isinstance(chrom, str):
        chrom = [chrom]
    if pos is not None and isinstance(pos, str):
        pos = [tuple([int(x) for x in pos.split("-")])]
    if chrom is not None and pos is not None:
        assert len(chrom) == len(pos), "If you are specifying chrom and pos at the same time, make sure they have the same length"
        pos = [tuple([int(y) for y in x.split('-')]) if isinstance(x, str) else x for x in pos]
        title_result.extend([f"{x}:{y[0]}-{y[1]}" for x, y in zip(chrom, pos)])
        chrom_result.extend(chrom)
        pos_result.extend(pos)
    elif chrom is not None and pos is None:
        title_result.extend(chrom)
        chrom_result.extend(chrom)
        for c in chrom: 
            end_pos = pd.read_csv(reference_fai, sep='\t', header=None, index_col=0, usecols=[0,1], nrows=24, names=['CHROM', 'size_int'] ).loc[c].squeeze()
            pos_result.append((0, end_pos))
    return title_result, chrom_result, pos_result
                

def genome_figure(coord: Union[str,list]=None, chrom: Union[str,list]=None, pos: Union[str,list]=None, 
                    figsize=None, figheight=0.5, figx=25_000_000, padding: float=None, axes_ratio: list=None,
                    n_samples=1, heatmap=False, 
                    reference_fai="/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai"):
    '''
    chrom: list of strings where each string is a chromosome to draw
    pos: list of tuples. 
    padding: must be in inch! (same unit as figsize)
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
    title_result, chrom_result, pos_result = process_chrom_pos_coord(coord=coord, chrom=chrom, pos=pos, reference_fai=reference_fai)

    if len(title_result):
        plot_type = 'specific plot'
        pos_size_list = np.array(list(map(lambda x: abs(x[0]-x[1]), pos_result)))
    else : # coord is None and chrom is None and pos is None (default setting: draw all chromosomes with each of them covering whole ranges)
        plot_type = 'genome plot'
        df_fai = pd.read_csv("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai", sep='\t', header=None, index_col=0, usecols=[0,1], nrows=24, names=['CHROM', 'size_int'] )
        title_result = df_fai.index
        chrom_result = title_result
        pos_size_list = df_fai['size_int'].to_numpy()
        pos_result = [(0, x) for x in pos_size_list]
    if axes_ratio is not None:
        pos_size_list_normalized = np.array(axes_ratio)/np.array(axes_ratio).sum()
    else:
        pos_size_list_normalized = pos_size_list/pos_size_list.sum()
    pos_size_list_normalized_sum = np.insert(np.cumsum(pos_size_list_normalized), 0, 0)[:-1]

    if figsize:
        figsize = figsize
        if plot_type == 'specific plot' and padding is None:
            padding = max(0.1 / figsize[0] / len(title_result), 0.05 / figsize[0])
        elif padding is not None:
            padding = padding/figsize[0]
    else:
        padding = 0 if padding is None else padding
        figwidth = pos_size_list.sum()/figx + float(padding)*(len(title_result)-1)
        figsize = (figwidth, figheight)
        padding = padding/figsize[0]

    fig = plt.figure(figsize=figsize)
    axes_all = []
    # for i in range(n_samples):
    #     axes_all.append([])

    ax_height_total = 0.85
    ax_height_each = ax_height_total/n_samples
    for j in range(n_samples):
        axes = []
        for i, (title, chrom, pos) in enumerate(zip(title_result, chrom_result, pos_result)):
            start = pos_size_list_normalized_sum[i]
            end = pos_size_list_normalized[i]
            if i == 0:
                ax = Axes_genome.new(fig.add_axes([start, j*ax_height_each, end, ax_height_each]), title, chrom, pos, leading=True)
                ax.xaxis.set(ticks=[pos[0], pos[1]])
                if heatmap:
                    ax.yaxis.set(ticks=[])
                axes.append(ax)
            else:
                padding = 0 if plot_type == 'genome plot' else padding
                end = end - padding if end - padding > 0 else 0
                ax = Axes_genome.new(fig.add_axes([start + padding, j*ax_height_each, end, ax_height_each], sharey=axes[0]), title, chrom, pos)
                if plot_type == 'genome plot':
                    ax.xaxis.set(ticks=[pos_size_list[i]])
                ax.spines['left'].set_visible(False)
                ax.axes.yaxis.set_visible(False)
                axes.append(ax)
            if plot_type == 'specific plot':
                ax.xaxis.set_major_locator(tick.LinearLocator(numticks=max(round(pos_size_list_normalized[i]/0.1), 2)))
            if i % 2 == 1:
                ax.tick_params(pad=20, axis='x')
                if plot_type == 'genome plot' and not heatmap:
                    ax.set_facecolor('gainsboro')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if j > 0:
                ax.spines['bottom'].set_visible(False)
                ax.get_xaxis().set_visible(False)
            if j == n_samples - 1:
                ax.set_title(ax.genome_title)
            ax.xaxis.set_major_formatter(tick.FuncFormatter(reformat_large_tick_values))
        axes_all.append(axes)
    if n_samples == 1:
        axes_all = axes_all[0]

    return fig, axes_all


def genome_figure_(coord: Union[str,list]=None, chrom: Union[str,list]=None, pos: Union[str,list]=None, 
                    figsize=(20, 4), padding: float=None, axes_ratio: list=None,
                    reference_fai="/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai"):
    '''
    legacy function
    
    chrom: list of strings where each string is a chromosome to draw
    pos: list of tuples. 

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
    axes = []
    title_result, chrom_result, pos_result = process_chrom_pos_coord(coord=coord, chrom=chrom, pos=pos, reference_fai=reference_fai)

    if len(title_result):
        plot_type = 'specific plot'
        pos_size_list = np.array(list(map(lambda x: abs(x[0]-x[1]), pos_result)))
    else : # coord is None and chrom is None and pos is None (default setting: draw all chromosomes with each of them covering whole ranges)
        plot_type = 'genome plot'
        df_fai = pd.read_csv("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai", sep='\t', header=None, index_col=0, usecols=[0,1], nrows=24, names=['CHROM', 'size_int'] )
        title_result = df_fai.index
        chrom_result = title_result
        pos_size_list = df_fai['size_int'].to_numpy()
        pos_result = [(0, x) for x in pos_size_list]
    if axes_ratio is not None:
        pos_size_list = np.array(axes_ratio)
    pos_size_list_normalized = pos_size_list/pos_size_list.sum()
    pos_size_list_normalized_sum = np.insert(np.cumsum(pos_size_list_normalized), 0, 0)[:-1]
    if plot_type == 'specific plot':
        padding = min(0.05 / len(title_result), 0.005)
    elif padding is not None:
        padding = padding

    for i, (title, chrom, pos) in enumerate(zip(title_result, chrom_result, pos_result)):
        start = pos_size_list_normalized_sum[i]
        end = pos_size_list_normalized[i]
        if i == 0:
            ax = Axes_genome.new(fig.add_axes([start, 0, end, 0.85]), title, chrom, pos, leading=True)
            ax.xaxis.set(ticks=[pos[0], pos[1]])
            axes.append(ax)
        else:
            padding = 0 if plot_type == 'genome plot' else padding
            ax = Axes_genome.new(fig.add_axes([start + padding, 0, end - padding, 0.85], sharey=axes[0]), title, chrom, pos)
            if plot_type == 'genome plot':
                ax.xaxis.set(ticks=[pos_size_list[i]])
            ax.spines['left'].set_visible(False)
            ax.axes.yaxis.set_visible(False)
            axes.append(ax)
        if plot_type == 'specific plot':
            ax.xaxis.set_major_locator(tick.LinearLocator(numticks=max(round(pos_size_list_normalized[i]/0.1), 2)))
        if i % 2 == 1:
            ax.tick_params(pad=20, axis='x')
            if plot_type == 'genome plot':
                ax.set_facecolor('gainsboro')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(ax.genome_title)
        ax.xaxis.set_major_formatter(tick.FuncFormatter(reformat_large_tick_values))
    return fig, axes


def reformat_large_tick_values(tick_val, pos):
    """
    Turns large tick values (in the billions, millions and thousands) such as 4500 into 4.5K and also appropriately turns 4000 into 4K (no zero after the decimal).
    https://dfrieds.com/data-visualizations/how-format-large-tick-values.html
    """
    if tick_val >= 1000000000:
        val = round(tick_val/1000000000, 1)
        new_tick_format = '{:}B'.format(val)
    elif tick_val >= 1000000:
        val = round(tick_val/1000000, 1)
        new_tick_format = '{:}M'.format(val)
    elif tick_val >= 1000:
        val = round(tick_val/1000, 1)
        new_tick_format = '{:}K'.format(val)
    elif tick_val < 1000:
        new_tick_format = round(tick_val, 1)
    else:
        new_tick_format = tick_val

    # make new_tick_format into a string value
    new_tick_format = str(new_tick_format)
    
    # code below will keep 4.5M as is but change values such as 4.0M to 4M since that zero after the decimal isn't needed
    index_of_decimal = new_tick_format.find(".")
    
    if index_of_decimal != -1:
        value_after_decimal = new_tick_format[index_of_decimal+1]
        if value_after_decimal == "0":
            # remove the 0 after the decimal point since it's not needed
            new_tick_format = new_tick_format[0:index_of_decimal] + new_tick_format[index_of_decimal+2:]            
    return new_tick_format


def sanity_check(df, ax, *args):
    assert isinstance(df, pd.core.frame.DataFrame), "input df must be a form of pandas df"
    assert isinstance(ax, Axes_genome), "input ax must be a form of genome axes"
    for k in args:
        if isinstance(k, str):
            assert k in df.columns, f'the key {k} is not in the dataframe'
        elif isinstance(k, int):
            assert k < len(df.columns), f'specified key index {k} is not in the dataframe'
        else:
            raise ValueError("keys must be either a string or an integer")


def key_check(df, k):
    if isinstance(k, int):
        return df.columns[k]
    else:
        return k


def check_pos_in_ax(chrom, pos, *axes):
    result = False
    for ax in axes:
        result = (ax.chrom == chrom and pos in ax.pos_range) or result
    return result


def check_axes_overlap(ax1, ax2):
    return ax1.chrom == ax2.chrom and (ax1.start+1 in ax2.pos_range or ax1.end in ax2.pos_range or ax2.start+1 in ax1.pos_range or ax2.end in ax1.pos_range)


def draw_scatter(df, ax, chrom: Union[str, int]='CHROM', pos: Union[str, int]='POS', **kwargs):
    columns_to_check = [x for x in [chrom, pos] if x is not None]
    sanity_check(df, ax, *columns_to_check)
    chrom = key_check(df, chrom)
    pos = key_check(df, pos)
    filter_bool = (df[chrom] == ax.chrom) & (df[pos].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both'))
    sns.scatterplot(data=df[filter_bool], ax=ax, x=pos, **kwargs)
    if not ax.leading and ax.get_legend():
        ax.get_legend().remove()
    elif ax.leading and ax.get_legend():
        ax.legend(bbox_to_anchor=(1.0, 0.85), loc='upper left', bbox_transform=ax.get_figure().transFigure)
    elif not ax.get_legend():
        pass
    ax.set_xlabel("")
    return None


def draw_bed(df, ax, chrom: Union[str, int]='CHROM', start: Union[str, int]='start', end: Union[str, int]='end', y: Union[str, int]='CNt', **kwargs):
    columns_to_check = [x for x in [chrom, start, end, y] if x is not None]
    sanity_check(df, ax, *columns_to_check)
    chrom = key_check(df, chrom)
    start = key_check(df, start)
    end = key_check(df, end)
    y = key_check(df, y)
    filter_bool = (df[chrom] == ax.chrom) & (df[start].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both') | \
                                             df[end].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both') | \
                                             (df[start] <= min(ax.start,ax.end)) & (df[end] >= max(ax.start,ax.end)))
    X, Y = [], []
    for s, e, c in zip(df[filter_bool][start], df[filter_bool][end], df[filter_bool][y]):
        X.extend([s, e, None]), Y.extend([c, c, None])
    ax.plot(X, Y, **kwargs)
    return None


def draw_bnd(df, ax, chrom, pos, color_key=None, ymin=0, ymax=1, **kwargs):
    chrom = key_check(df, chrom)
    pos = key_check(df, pos)
    filter_bool = (df[chrom] == ax.chrom) & (df[pos].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both'))
    pos_array = df[filter_bool][pos]
    if color_key:
        color_array = df[filter_bool][color_key].map(sv_color_dict)
    else:
        color_array = np.repeat('black', len(pos_array))
    for x, c in zip(pos_array, color_array):
        ax.axvline(x=x, c=c, **kwargs)
    return


def draw_highlight(df, ax, chrom, start, end, color_key=None, cmap=None, ymin=0, ymax=1, **kwargs):
    chrom = key_check(df, chrom)
    start = key_check(df, start)
    end = key_check(df, end)
    filter_bool = (df[chrom] == ax.chrom) & (df[start].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both') | df[end].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both'))
    if cmap is None:
        cmap = sns.color_palette("coolwarm", as_cmap=True)
    if color_key is not None:
        color_key = key_check(df, color_key)
        for xmin, xmax, color in zip(df[filter_bool][start], df[filter_bool][end], cmap(df[filter_bool][color_key])):
            ax.axvspan(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, **kwargs)
        return
    else:
        for xmin, xmax in zip(df[filter_bool][start], df[filter_bool][end]):
            ax.axvspan(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, **kwargs)
        return




def coord_transform_data_to_fig(coord, ax, fig):
    ax_to_fig = ax.transData + fig.transFigure.inverted()
    coord = ax_to_fig.transform(coord)
    return coord

def coord_transform_data_to_dpi(coord, ax):
    return ax.transData.transform(coord)

def coord_transform_dpi_to_fig(coord, fig):
    return fig.transFigure.inverted().transform(coord)

def coord_to_beziers(coord1, coord2, curve_type='cubic', convex_up=True, normalize=False, scale=0.8):
    assert np.any(coord1 != coord2), 'you cannot draw an arc between the same points'
    coord1, coord2 = np.array(coord1), np.array(coord2)
    slope = coord2 - coord1
    if (slope[0] >= 0 and convex_up) or (slope[0] < 0 and not convex_up): # when to select counter clockwise norm
        norm = np.array((-slope[1], slope[0]))
    elif (slope[0] < 0 and convex_up) or (slope[0] >= 0 and not convex_up): # when to select counter clockwise norm
        norm = np.array((slope[1], -slope[0]))
    else:
        raise ValueError('something wrong with coordinates')
    if normalize:
        norm = norm/np.sqrt(np.sum(norm**2)) * 0.1
    if curve_type == 'cubic':
        coords = [coord1, coord1+0.2*slope+scale*norm, coord1+0.8*slope+scale*norm, coord2]
    elif curve_type == 'quadratic':
        coords = [coord1, coord1+0.5*slope+scale*norm, coord2]
    return coords


def beziers_to_arc_artist(fig=None, *coords, **kwargs):
    if len(coords) == 3:
        codes = [1, 3, 3]  
    elif len(coords) == 4:
        codes = [1, 4, 4, 4]
    else:
        raise ValueError("number of coords not supported")
    path_object = Path(vertices=coords, codes=codes)
    patch_return = patches.PathPatch(path_object, facecolor='none', transform=fig.transFigure , **kwargs)
    return patch_return


def coord_to_ticks(coord1, coord2, t, fig, scale=0.05, **kwargs):
    offset = np.array([scale, 0])
    if t == '3to5': #DEL
        coords = [coord1-offset, coord1, coord2, coord2+offset]
    elif t == '5to3': #DUP
        coords = [coord1, coord1+offset, coord2-offset, coord2]
    elif t == '3to3': #INV
        coords = [coord1-offset, coord1, coord2-offset, coord2]
    elif t == '5to5': #INV
        coords = [coord1, coord1+offset, coord2, coord2+offset]
    path_object = Path(vertices=coords, codes=[1,2,1,2])
    patch_return = patches.PathPatch(path_object, facecolor='none', transform=fig.transFigure, **kwargs)
    return patch_return


def draw_arcs_on_axes(df, ax1, ax2, chrom1='CHR1', pos1='POS1', y1=1, chrom2='CHR2', pos2='POS2', y2=1, hue=None, normalize=False, **kwargs):
    columns_to_check = [x for x in [chrom1, pos1, y1, chrom2, pos2, y2] if x is not None]
    sanity_check(df, ax1, *columns_to_check)
    sanity_check(df, ax2)
    assert ax1.get_figure() == ax2.get_figure(), "ax1 and ax2 are in different fig"
    if check_axes_overlap(ax1, ax2):
        print_err(f'{ax1.__repr__()} and {ax2.__repr__()} overlap')
        print_err('This overlap may cause an unexpected arc to show up')
    fig = ax1.get_figure()
    chrom1 = key_check(df, chrom1)
    pos1 = key_check(df, pos1)
    y1 = key_check(df, y1)
    chrom2 = key_check(df, chrom2)
    pos2 = key_check(df, pos2)
    y2 = key_check(df, y2)
    filter_bool = df.apply(lambda x: check_pos_in_ax(x[chrom1], x[pos1], ax1, ax2), axis=1) | df.apply(lambda x: check_pos_in_ax(x[chrom2], x[pos2], ax1, ax2), axis=1)
    df = df[filter_bool]
    for i in range(len(df)):
        coord1_x, coord1_y = int(df.iloc[i][pos1]), float(df.iloc[i][y1])
        coord2_x, coord2_y = int(df.iloc[i][pos2]), float(df.iloc[i][y2])
        for ax in [ax1, ax2]:
            if check_pos_in_ax(df.iloc[i][chrom1], coord1_x, ax):
                ax_of_coord1 = ax
            if check_pos_in_ax(df.iloc[i][chrom2], coord2_x, ax):
                ax_of_coord2 = ax
        coord1, coord2 = (coord1_x, coord1_y), (coord2_x, coord2_y)
        coord1, coord2 = coord_transform_data_to_fig(coord1, ax_of_coord1, fig), coord_transform_data_to_fig(coord2, ax_of_coord2, fig)
        coordinates = coord_to_beziers(coord1, coord2, normalize=normalize)
        patch = beziers_to_arc_artist(fig, *coordinates, **kwargs)
        fig.add_artist(patch)
    return None


def draw_sv_on_fig(df, fig, chrom1='CHR1', pos1='POS1', cn1=None,
                            chrom2='CHR2', pos2='POS2', cn2=None, 
                            hue=None, hue_dict=None, tick=None, tick_dict=None, 
                            curve=None, curve_dict=None, normalize=False, scale_arc=0.8, scale_tick=0.05, **kwargs):
    axes = fig.axes
    axes_pr = pr.PyRanges(chromosomes=[ax.chrom for ax in axes], starts=[ax.start for ax in axes], ends=[ax.end for ax in axes])
    assert len(axes) == len(axes_pr.merge()), "some axes overlap"
    filter_bool = df.apply(lambda x: check_pos_in_ax(x[chrom1], x[pos1], *axes), axis=1) & df.apply(lambda x: check_pos_in_ax(x[chrom2], x[pos2], *axes), axis=1)
    df = df[filter_bool]


    if hue is not None:
        if hue in df.columns and hue_dict is not None:
            df[hue] = df[hue].map(hue_dict)
        elif isinstance(hue, str):
            df[hue] = hue
        else:
            df['hue'] = hue
            hue = 'hue'
    else:
        df['hue'] = 'black'
        hue='hue'

    if tick_dict is not None:
        df[tick] = df[tick].map(dict_dict)

    if curve_dict is not None:
        df[curve] = df[curve].map(curve_dict)

    for c1, x1, y1, c2, x2, y2, h, t, u in zip(df[chrom1], df[pos1], df[cn1], df[chrom2], df[pos2], df[cn2], df[hue], df[tick], df[curve]):
        for ax in axes:
            if check_pos_in_ax(c1, x1, ax):
                ax1 = ax
            if check_pos_in_ax(c2, x2, ax):
                ax2 = ax
        coord1, coord2 = coord_transform_data_to_fig((x1, y1), ax1, fig), coord_transform_data_to_fig((x2, y2), ax2, fig)
        coordinates = coord_to_beziers(coord1, coord2, convex_up=u, scale=scale_arc, normalize=normalize)
        arc_patch = beziers_to_arc_artist(fig, *coordinates, edgecolor=h, **kwargs)
        fig.add_artist(arc_patch)
        tick_patch = coord_to_ticks(coord1, coord2, t, fig, edgecolor=h, scale=scale_tick, lw=3)
        fig.add_artist(tick_patch)
    return




from matplotlib_venn import venn2
from matplotlib_venn import venn3
from venn import venn
from venn import pseudovenn

def vennk(*args, label:list=[], mode:str='pyvenn'):
    '''
    args: input sets positionaly 
    label: label for each set
    mode: either 'pyvenn' or 'plt' or 'pseudo'
    '''
    assert len(args) == len(label) # The number of labels and the number of sets must be the same
    input_list = [set(x) for x in args]
    if mode == 'plt':
        if len(input_list) > 3:
            print("'plt' mode does not support venn diagram with more than 3 sets")
            return
        elif len(args) == 3:
            venn3(subsets=input_list, set_labels=label)
        elif len(args) == 2:
            venn2(subsets=input_list, set_labels=label)
        else:
            print("more than one set is required")
        return
    
    elif mode == 'pyvenn':
        venn({k:v for k,v in zip(label, input_list)})
        return
    elif mode == 'pseudo':
        pseudovenn({k:v for k,v in zip(label, input_list)})
        return
    else:
        print('unknown mode')
        return

        

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



# dictionary for signature colors
sig_color_dict = {'SBS1': '#e32636',
'SBS2': '#efdecd',
'SBS3': '#e52b50',
'SBS4': '#ffbf00', 
'SBS5': '#5d8aa8',
'SBS6': '#9966cc',
'SBS7a': '#a4c639',
'SBS7b': '#f2f3f4',
'SBS7c': '#cd9575', 
'SBS7d': '#915c83',
'SBS8': '#faebd7',
'SBS9': '#008000',
'SBS10a': '#00ffff', 
'SBS10b': '#4b5320',
'SBS10c': '#e9d66b',
'SBS10d': '#b2beb5',
'SBS11': '#87a96b',
'SBS12': '#ff9966',
'SBS13': '#a52a2a',
'SBS14': '#fdee00',
'SBS15': '#ff2052', 
'SBS16': '#007fff', 
'SBS17a': '#138808',
'SBS17b': '#f4c2c2', 
'SBS18': '#21abcd', 
'SBS19': '#ffe135',
'SBS20': '#98777b',
'SBS21': '#f5f5dc',
'SBS22': '#3d2b1f',
'SBS23': '#fe6f5e', 
'SBS24': '#0000ff',
'SBS25': '#0000ff', 
'SBS26': '#8a2be2', 
'SBS27': '#de5d83', 
'SBS28': '#e3dac9',
'SBS29': '#cc0000',
'SBS30': '#873260',
'SBS31': '#b5a642',
'SBS32': '#cb4154',
'SBS33': '#66ff00',
'SBS34': '#c32148',
'SBS35': '#ff007f',
'SBS36': '#f4bbff', 
'SBS37': '#fb607f', 
'SBS38': '#004225',
'SBS39': '#a52a2a',
'SBS40': '#0047ab',
'SBS41': '#f0dc82',
'SBS42': '#800020',
'SBS43': '#cc5500',
'SBS44': '#e97451',
'SBS45': '#5f9ea0',
'SBS46': '#91a3b0',
'SBS47': '#a67b5b',
'SBS48': '#4b3621', 
'SBS49': '#78866b',
'SBS50': '#ffff99',
'SBS51': '#ffef00',
'SBS52': '#00bfff',
'SBS53': '#c41e3a',
'SBS54': '#00cc99',
'SBS55': '#ff0038',
'SBS56': '#ed9121',
'SBS57': '#ace1af',
'SBS58': '#b2ffff',
'SBS59': '#4997d0',
'SBS60': '#7fff00',
'SBS84': '#cd5c5c',
'SBS85': '#fafad2',
'SBS86': '#d3d3d3',
'SBS87': '#ffb6c1',
'SBS88': '#ff9999',
'SBS89': '#20b2aa',
'SBS90': '#c8a2c8',
'SBS91': '#0abab5',
'SBS92': '#e08d3c',
'SBS93': '#ff6347',
'SBS94': '#746cc0',
'SBS_ox': '#ff6633',
'SBS_cisox': '#ff9933',
'ID1': '#cc3333',
'ID2': '#df00ff',
'ID3': '#000f89',
'ID4': '#01796f',
'ID5': '#93c572',
'ID6': '#1fcecb', 
'ID7': '#bc8f8f',
'ID8': '#7851a9',
'ID9': '#bb6528',
'ID10': '#a81c07',
'ID11': '#f4c430',
'ID12': '#ff91a4',
'ID13': '#507d2a',
'ID14': '#006d5b',
'ID15': '#321414',
'ID16': '#87ceeb',
'ID17': '#f28500',
'ID18': '#367588',
'ID-A': '#ff0033'
}