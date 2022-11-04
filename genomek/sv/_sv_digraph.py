import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import pyranges as pr
import tempfile
from collections import defaultdict
from ..tools._gadgets import chrom_sort_dict_37 as chrom_sort_dict
from ._sv import get_svtype

sv_hue_dict = {'DUP': 'green', 'TRA':'purple', 'DEL':'red', 'INV':'blue'}

def sv_to_bps(df, fai_path="/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai"):
    result = list(zip(df['CHR1'], df['POS1'])) + list(zip(df['CHR2'], df['POS2']))
    df_fai = pd.read_csv(fai_path, sep='\t', header=None, index_col=None, usecols=[0,1], nrows=24, names=['CHR', 'POS'] )
    chrom_ends = list(zip(df_fai['CHR'], df_fai['POS']))
    chrom_ends = chrom_ends + [(x,0) for x,y in chrom_ends]
    result += chrom_ends
    result = list(set(result))
    result = sorted(result, key=lambda x: (chrom_sort_dict[str(x[0])], int(x[1])))
    return result

def bps_to_pr(node_list):
    range_list = list(zip(*[(node_list[i][0], node_list[i][1], node_list[i+1][1]) for i in range(len(node_list)-1) if node_list[i][0] == node_list[i+1][0]]))
    P = pr.PyRanges(chromosomes=range_list[0], starts=range_list[1], ends=range_list[2])
    return P

def bps_to_G(node_list):
    G = nx.DiGraph()
    i = 0
    for chrom in chrom_sort_dict:
        i_start = i
        for CHR, POS in filter(lambda x: x[0] == chrom, node_list):            
            node1 = (2*i, {'CHR': str(CHR), 'POS': int(POS), 'coord': f'{CHR}:{POS}', 'dir': '5'})
            node2 = (2*i+1,  {'CHR': str(CHR), 'POS': int(POS), 'coord': f'{CHR}:{POS}', 'dir': '3'})
            G.add_nodes_from([node1, node2])
            i += 1
        G.add_edges_from([(2*i, 2*(i+1)) for i in range(i_start, i-1)]+[(2*(i+1)+1, 2*i+1) for i in range(i_start, i-1)], segment=True, read=False) 
    return G

def find_node(G, CHR, POS):
    if len(result := [n for n,d in G.nodes(data=True) if d['CHR'] == CHR and d['POS'] == POS]) == 2:
        return result
    else:
        raise ValueError("Can't find the node")

def add_sv_from_df(G, df, *args):
    for row in df.itertuples():
        node1 = find_node(G, row.CHR1, row.POS1)
        n_1_5 = [x for x in node1 if x % 2 == 0][0]
        n_1_3 = [x for x in node1 if x % 2 == 1][0]
        node2 = find_node(G, row.CHR2, row.POS2)
        n_2_5 = [x for x in node2 if x % 2 == 0][0]
        n_2_3 = [x for x in node2 if x % 2 == 1][0]
        CT = row.CT
        common_dict = {k:row[df.columns.tolist().index(k) + 1] for k in args}
        if CT == '5to3':
            G.add_edge(n_2_5, n_1_5, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
            G.add_edge(n_1_3, n_2_3, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
        elif CT == '3to5':
            G.add_edge(n_1_5, n_2_5, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
            G.add_edge(n_2_3, n_1_3, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
        elif CT == '3to3':
            G.add_edge(n_1_5, n_2_3, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
            G.add_edge(n_2_5, n_1_3, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
        elif CT == '5to5':
            G.add_edge(n_2_3, n_1_5, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
            G.add_edge(n_1_3, n_2_5, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
        else:
            raise ValueError("Wrong CT")        
    return G

def bed_to_cov_df(bed_path, bam_path_t, bam_path_n, df_key):
    df_mos_t, avg_t, mt_t = genomek.cov.bam_to_cov(bam_path_t, bed_path)
    df_mos_n, avg_n, mt_n = genomek.cov.bam_to_cov(bam_path_n, bed_path)
    df_mos_t['cov_n'] = df_mos_n['cov']
    df_mos_t[df_key] = df_mos_t.apply(lambda x: genomek.cov.absolute_copy(x['cov'], x['cov_n'], avg_t, avg_n, 2, 1), axis=1)
    return df_mos_t

def add_weight_from_df(G, df, df_key, weight_key):
    '''
    df should have CHROM start end
    '''
    for row in df.itertuples():
        node1 = find_node(G, row.CHROM, row.start)
        n_1_5 = [x for x in node1 if x % 2 == 0][0]
        n_1_3 = [x for x in node1 if x % 2 == 1][0]
        node2 = find_node(G, row.CHROM, row.end)
        n_2_5 = [x for x in node2 if x % 2 == 0][0]
        n_2_3 = [x for x in node2 if x % 2 == 1][0]
        weight = row[df.columns.tolist().index(df_key) + 1]
        G[n_1_5][n_2_5][weight_key] = weight
        G[n_2_3][n_1_3][weight_key] = weight     
    return G

def df_to_diG(df, *args):
    node_list = sv_to_bps(df)
    G = bps_to_G(node_list)
    G = add_sv_from_df(G, df, *args)
    return G

def df_to_cov_df(df, bam_path_t, bam_path_n, df_key):
    node_list = sv_to_bps(df)
    with tempfile.TemporaryDirectory() as tempdir:
        bed_path = f"{tempdir}/temp.bed"
        bps_to_pr(node_list).to_bed(bed_path, keep=False)
        cov_df = bed_to_cov_df(bed_path, bam_path_t=bam_path_t, bam_path_n=bam_path_n, df_key=df_key)
    return cov_df

def draw_diG(G, figsize=(6,6), pos_linear=False):
    fig, ax = plt.subplots(figsize=figsize)
    if pos_linear:
        pos_dict = defaultdict(set)
        for n,d in G.nodes(data=True):
            pos_dict[d['CHR']].add(d['POS'])
        pos_dict = pd.Series({k:max(v) for k,v in pos_dict.items()}).shift(1).fillna(0)
        adjust_pad = pos_dict.mean()/len(pos_dict)
        pos_dict = (pos_dict.cumsum() + pd.Series([i*adjust_pad for i in range(len(pos_dict))], index=pos_dict.index)).to_dict()
        node_pos = {**{x:(G.nodes[x]['POS']+pos_dict[G.nodes[x]['CHR']], 2) for x in filter(lambda x: x % 2 == 0, G.nodes)},
                    **{x:(G.nodes[x]['POS']+pos_dict[G.nodes[x]['CHR']], 1) for x in filter(lambda x: x % 2 == 1, G.nodes)}}
        node_labels = None
    else:
        node_pos = {**nx.circular_layout(filter(lambda x: x % 2 == 0, G.nodes), scale=1.1, center=(0,0)),
                    **nx.circular_layout(filter(lambda x: x % 2 == 1, G.nodes), scale=0.9, center=(0,0))}
        node_labels = {n:d['coord'] for n,d in G.nodes(data=True)}
#     node_colors = [chrom_color_dict[d['CHR']] for n, d in G.nodes(data=True)]
    node_colors = 'gray'
    seg_list = [(u,v) for u,v,e in G.edges(data=True) if e['segment']]
    read_list = [(u,v) for u,v,e in G.edges(data=True) if e['read']]
    read_colors = [sv_hue_dict[e['SVTYPE']] for u,v,e in G.edges(data=True) if e['read']]

    nx.draw_networkx_nodes(G, node_pos, ax=ax, node_color=node_colors, node_size=50, alpha=1)
    nx.draw_networkx_labels(G, node_pos, ax=ax, labels=node_labels, font_size=8)
    nx.draw_networkx_edges(G, node_pos, ax=ax, edgelist=seg_list, node_size=50, edge_color='black', width=0.5)
    for read, read_color in zip(read_list, read_colors):
        arrowprop=dict(linestyle='-', color=read_color, connectionstyle="arc3,rad=-0.3", width=0.5, headwidth=4, headlength=4)
        ax.annotate("", xy=node_pos[read[1]], xytext=node_pos[read[0]], arrowprops=arrowprop)
    plt.show()
    return

# G = df_to_diG(df_sv_test)
# draw_diG(G.subgraph((n for n,d in G.nodes(data=True) if d['CHR'] == '1' or d['CHR'] == '2' or d['CHR'] == '3')), pos_linear=False, figsize=(6,6))
# draw_diG(G.subgraph((n for n,d in G.nodes(data=True) if d['CHR'] == '1' or d['CHR'] == '2' or d['CHR'] == '3')), pos_linear=True, figsize=(12,1))
# G = df_to_diG(df_sv, 'score:::O', 'TIER1:::O', 'vaf_1:::O', 'vaf_2:::O')
# draw_diG(G.subgraph((n for n,d in G.nodes(data=True) if d['CHR'] == '1' or d['CHR'] == '2' or d['CHR'] == '3')), pos_linear=False, figsize=(6,6))


def G_to_seg_df(G):
    result = {}
    G_5prime = G.subgraph([n for n,d in G.nodes(data=True) if d['dir'] == '5'])
    for u,v in ((u,v) for u,v,d in G_5prime.edges(data=True) if d['segment']):
        assert G_5prime.nodes[u]['CHR'] == G_5prime.nodes[v]['CHR']
        result[u] = {'chromosome':G_5prime.nodes[u]['CHR'], 'start':G_5prime.nodes[u]['POS'], 'end':G_5prime.nodes[v]['POS']}
    return pd.DataFrame.from_dict(result, orient='index')