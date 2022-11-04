import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import FancyBboxPatch
from ._basic_plots import draw_highlight
from ._basic_plots import Axes_genome

chromosomes_37 = [str(i) for i in list(range(1, 23)) + ['X', 'Y']]
chromosomes_38 = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y']]
chrom_color = ['#6266ff', #1
               '#d6503b', #2
               '#88ab62', #3
               '#f3ea8b', #4
               '#587d98', #5
               '#c6763f', #6
               '#71bfe7', #7
               '#91317f', #8
               '#80dd6f', #9
               '#dda6b6', #10
               '#a35b23', #11
               '#968ea1', #12
               '#d16424', #13
               '#dca067', #14
               '#8c7aba', #15
               '#e9bd72', #16
               '#4c256b', #17
               '#d7e4c0', #18
               '#743b91', #19
               '#bc2e7a', #20
               '#ecd176', #21
               '#6d7870', #22
               '#d5a900', #X
               '#aad500', #Y
               ]
#https://software.broadinstitute.org/software/igv/interpreting_insert_size


def ideogram_ax_(ax, ymin=0, ymax=1, round_box=True, **kwargs):
    '''
    legacy function
    '''
    assert isinstance(ax, Axes_genome), "input ax must be a form of genome axes"
    df = pd.read_csv("/home/users/kimin/projects/00_Reference/cytoband/hg19.cytoband.processed.tsv.gz", sep='\t')
    chrom = 'chrom'
    start = 'start'
    end = 'end'
    color='color'
    filter_bool = (df[chrom] == ax.chrom) & (df[start].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both') | df[end].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both'))
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    if round_box:
        aspect_ratio = ax.transAxes.transform((1,1))
        aspect_ratio = aspect_ratio[0]/aspect_ratio[1]
        if aspect_ratio < 1:
            aspect_ratio = 1/aspect_ratio
        pad = 0.5
        y_pad = pad*aspect_ratio
        p_bbox = FancyBboxPatch((-pad, -y_pad), 1+2*pad, 1+2*y_pad,
                            boxstyle=f"round,pad=-{pad},rounding_size={pad/aspect_ratio}",
                            ec="black", fc="none", clip_on=False, lw=1,
                            mutation_aspect=aspect_ratio,
                            transform=ax.transAxes)
    else:
        p_bbox = FancyBboxPatch((0, 0), 1, 1,
                            boxstyle=f"square,pad=0",
                            ec="black", fc="none", clip_on=False, lw=1,
                            mutation_aspect=1,
                            transform=ax.transAxes)

    ax.add_patch(p_bbox)
    ax.patch = p_bbox
    for xmin, xmax, color, stain in zip(df[filter_bool][start], df[filter_bool][end], df[filter_bool][color], df[filter_bool]['stain']):
        if stain == 'acen1':
            ax.fill([xmin, xmax, xmin], [ymin, (ymin+ymax)/2, ymax], color, edgecolor=color, transform=trans)
        elif stain == 'acen2':
            ax.fill([xmin, xmax, xmax], [(ymin+ymax)/2, ymin, ymax], color, edgecolor=color, transform=trans)
        else:
            ax.axvspan(xmin=xmin, xmax=xmax, color=color, ymin=ymin, ymax=ymax, **kwargs)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axes.yaxis.set_visible(False)
    return


def ideogram_ax(ax, ymin=0, ymax=1, round_box=False, edge_width=2, use_chrom_color=True, ref_ver='37', **kwargs):
    '''

    '''
    # assert isinstance(ax, Axes_genome), "input ax must be a form of genome axes"
    df = pd.read_csv("/home/users/kimin/projects/00_Reference/cytoband/hg19.cytoband.processed.tsv.gz", sep='\t')
    chrom = 'chrom'
    start = 'start'
    end = 'end'
    color='color'
    filter_bool = (df[chrom] == ax.chrom) & (df[start].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both') | df[end].between(min(ax.start,ax.end), max(ax.start,ax.end), inclusive='both'))
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    if use_chrom_color:
        if ref_ver == '37':
            chrom_color_dict = {k:v for k,v in zip(chromosomes_37, chrom_color)}
        else:
            chrom_color_dict = {k:v for k,v in zip(chromosomes_38, chrom_color)}
        edge_color = chrom_color_dict[ax.chrom]
    else:
        edge_color = 'black'
    if round_box:
        aspect_ratio = ax.transAxes.transform((1,1))
        aspect_ratio = aspect_ratio[0]/aspect_ratio[1]
        if aspect_ratio < 1:
            aspect_ratio = 1/aspect_ratio
        pad = 0.5
        y_pad = pad*aspect_ratio
        p_bbox = FancyBboxPatch((-pad, -y_pad), 1+2*pad, 1+2*y_pad,
                            boxstyle=f"round,pad=-{pad},rounding_size={pad/aspect_ratio}",
                            ec=edge_color, fc="none", clip_on=False, lw=edge_width,
                            mutation_aspect=aspect_ratio,
                            transform=ax.transAxes)
        ax.patch = p_bbox
    else:
        p_bbox = FancyBboxPatch((0, ymin), 1, ymax-ymin,
                            boxstyle=f"square,pad=0",
                            ec=edge_color, fc="none", clip_on=False, lw=edge_width,
                            mutation_aspect=1,
                            transform=ax.transAxes)
    for xmin, xmax, color, stain in zip(df[filter_bool][start], df[filter_bool][end], df[filter_bool][color], df[filter_bool]['stain']):
        if stain == 'acen1':
            ax.fill([xmin, xmax, xmin], [ymin, (ymin+ymax)/2, ymax], color, edgecolor='none', transform=trans)
        elif stain == 'acen2':
            ax.fill([xmin, xmax, xmax], [(ymin+ymax)/2, ymin, ymax], color, edgecolor='none', transform=trans)
        else:
            ax.axvspan(xmin=xmin, xmax=xmax, facecolor=color, ymin=ymin, ymax=ymax, **kwargs)
    ax.add_patch(p_bbox)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_facecolor('none')
    return


def ideogram(fig, **kwargs):
    for ax in fig.get_axes():
        ideogram_ax(ax, **kwargs)
    return