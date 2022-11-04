import re
import pandas as pd
import networkx as nx
from itertools import combinations
from ..tools._gadgets import chromosomes_37 as chromosomes
from ..tools._gadgets import chrom_cat_type_37 as chrom_cat_type
from ..tools._gadgets import chrom_sort_dict_37 as chrom_sort_dict


def check_bps_make_semgent(n1, n2):
    if 'knot' in n1 and 'knot' in n2:
        return False
    if n1['CHR'] == n2['CHR']:
        if n1['3_open'] and n2['5_open'] and n1['POS'] < n2['POS']:
            return True
        elif n1['5_open'] and n2['3_open'] and n1['POS'] > n2['POS']:
            return True
    return False


def get_svtype(chr1, chr2, ct):
    if chr1 != chr2:
        svtype = 'TRA'
    else:
        if ct == '5to3':
            svtype = 'DUP'
        elif ct == '3to5':
            svtype = 'DEL'
        else:
            svtype = 'INV'
    return svtype


def add_1or2_to_keys(k, n):
    if ":::" in k:
        return re.sub(r":::", f"{n}:::", k)
    else:
        return f"{k}{n}"


def test_first_coord_is_proximal(chr1, pos1, chr2, pos2):
    if chrom_sort_dict[chr1] < chrom_sort_dict[chr2]:
        return True
    elif chrom_sort_dict[chr1] > chrom_sort_dict[chr2]:
        return False
    else:
        if int(pos1) <= int(pos2):
            return True
        else: 
            return False


def G_add_chrom_ends(G, fai_path="/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta.fai"):
    n_nodes = G.number_of_nodes()
    df_fai = pd.read_csv(fai_path, sep='\t', header=None, index_col=None, usecols=[0,1], nrows=24, names=['CHR', 'POS'] )
    n_index_shift = len(df_fai)
    for i, n in df_fai.to_dict(orient='index').items():
        index = n_nodes+i
        G.add_node(index, **{'CHR':n['CHR'], 'POS':int(n['POS']), 'coord': f"{n['CHR']}:{n['POS']}", '5_open':True, '3_open':False, 'knot':True})
        index = n_nodes+i+n_index_shift
        G.add_node(index, **{'CHR':n['CHR'], 'POS':0, 'coord': f"{n['CHR']}:0", '5_open':False, '3_open':True, 'knot':True})
    return G


def df_sv_to_G(df):
    G = nx.Graph()
    # G = nx.MultiGraph()
    for i, row in df.iterrows():
        index1 = i * 2
        index2 = index1+1
        bp1 = f"{row.CHR1}:{row.POS1}"
        bp2 = f"{row.CHR2}:{row.POS2}"
        ct1, ct2 = row.CT.split('to')
        ct1_5_open = False if ct1 == '5' else True
        ct1_3_open = False if ct1 == '3' else True
        ct2_5_open = False if ct2 == '5' else True
        ct2_3_open = False if ct2 == '3' else True
        node1 = (index1, {**{re.sub(r"_*1", "", k):v for k,v in row.filter(regex="^CHR|^POS|^vaf").filter(like='1').to_dict().items()}, **{'coord': bp1, '5_open': ct1_5_open, '3_open': ct1_3_open}})
        node2 = (index2, {**{re.sub(r"_*2", "", k):v for k,v in row.filter(regex="^CHR|^POS|^vaf").filter(like='2').to_dict().items()}, **{'coord': bp2, '5_open': ct2_5_open, '3_open': ct2_3_open}})
        common_dict = row.filter(regex="^TIER|^score|^CT").to_dict()
        G.add_nodes_from([node1, node2])
        G.add_edge(index1, index2, read=True, segment=False, SVTYPE=get_svtype(row.CHR1, row.CHR2, row.CT), **common_dict)
    G = G_add_chrom_ends(G)
    for n1, n2 in combinations(G.nodes, 2):
        if not G.has_edge(n1, n2) and check_bps_make_semgent(G.nodes[n1], G.nodes[n2]):
        # if check_bps_make_semgent(G.nodes[n1], G.nodes[n2]):
            G.add_edge(n1, n2, read=False, segment=True, CT=None)
    return G


def G_to_df_sv(G):
    result = {}
    for u,v,d in G.edges(data=True):
        if not d['read']: continue
        node1, node2 = G.nodes[u], G.nodes[v]
        if test_first_coord_is_proximal(node1['CHR'], node1['POS'], node2['CHR'], node2['POS']):
            prox_bp = {add_1or2_to_keys(k,"1"):v for k,v in node1.items()}
            dist_bp = {add_1or2_to_keys(k,"2"):v for k,v in node2.items()}
            index = u
        else:
            prox_bp = {add_1or2_to_keys(k,"1"):v for k,v in node2.items()}
            dist_bp = {add_1or2_to_keys(k,"2"):v for k,v in node1.items()}
            index = v
        result[index] = {**prox_bp, **dist_bp, **d}
    df = pd.DataFrame.from_dict(result, orient='index')
    df = df.astype({'CHR1': chrom_cat_type, 'CHR2': chrom_cat_type})
    essential_colnames = ['CHR1', 'POS1', 'CHR2', 'POS2', 'CT', 'SVTYPE']
    df = df[essential_colnames + [x for x in df.columns if x not in essential_colnames]].sort_values(by=['CHR1', 'CHR2', 'POS1', 'POS2'])
    return df


def find_all_cycles(G, source=None):
    """
    https://gist.github.com/joe-jordan/6548029
    edit for sv G
    """
    if source is None:
        # produce edges for all components
        # nodes=[i for i in G.nodes][0] # if n connected components of G is already 1
        nodes=[list(i)[0] for i in nx.connected_components(G)]
    else:
        # produce edges for components with source
        nodes=[source]
    # extra variables for cycle detection:
    cycle_stack = []
    output_cycles = set()
    
    def get_hashable_cycle(cycle):
        """cycle as a tuple in a deterministic order."""
        m = min(cycle)
        mi = cycle.index(m)
        mi_plus_1 = mi + 1 if mi < len(cycle) - 1 else 0
        if cycle[mi-1] > cycle[mi_plus_1]:
            result = cycle[mi:] + cycle[:mi]
        else:
            result = list(reversed(cycle[:mi_plus_1])) + list(reversed(cycle[mi_plus_1:]))
        return tuple(result)

    def change_filter_node_key(k):
        if k == 'read': return 'segment'
        elif k == 'segment': return 'read'
        else: raise ValueError()
    
    for start in nodes:
        if start in cycle_stack:
            continue
        cycle_stack.append(start)
        
        stack = [(start, ((u,v,e) for u,v,e in G.edges(start, data=True)))]
        stack_is_seg = []
        edge_list = []
        while stack:
            parent,children = stack[-1]
            try:
                child = next(children)
                if child[1] not in cycle_stack:
                    cycle_stack.append(child[1])
                    stack_is_seg.append(child[2]['segment'])
                    edge_list.append(child)
                    stack.append((child[1],((u,v,e) for u,v,e in G.edges(child[1], data=True) if e['segment'] is not child[2]['segment'])))
                else:
                    i = cycle_stack.index(child[1])
                    if i < len(cycle_stack) - 1: 
                        output_cycles.add(get_hashable_cycle(cycle_stack[i:]))
                
            except StopIteration:
                cycle_stack.pop()
                stack_is_seg.pop()
                edge_list.pop()
                stack.pop()
    
    return [list(i) for i in output_cycles]


def find_all_cycles_original(G, source=None, cycle_length_limit=None):
    """
    https://gist.github.com/joe-jordan/6548029
    forked from networkx dfs_edges function. Assumes nodes are integers, or at least
    types which work with min() and > .
    """
    if source is None:
        # produce edges for all components
        nodes=[list(i)[0] for i in nx.connected_components(G)]
    else:
        # produce edges for components with source
        nodes=[source]
    # extra variables for cycle detection:
    cycle_stack = []
    output_cycles = set()
    
    def get_hashable_cycle(cycle):
        """cycle as a tuple in a deterministic order."""
        m = min(cycle)
        mi = cycle.index(m)
        mi_plus_1 = mi + 1 if mi < len(cycle) - 1 else 0
        if cycle[mi-1] > cycle[mi_plus_1]:
            result = cycle[mi:] + cycle[:mi]
        else:
            result = list(reversed(cycle[:mi_plus_1])) + list(reversed(cycle[mi_plus_1:]))
        return tuple(result)
    
    for start in nodes:
        if start in cycle_stack:
            continue
        cycle_stack.append(start)
        
        stack = [(start,iter(G[start]))]
        while stack:
            parent,children = stack[-1]
            try:
                child = next(children)
                
                if child not in cycle_stack:
                    cycle_stack.append(child)
                    stack.append((child,iter(G[child])))
                else:
                    i = cycle_stack.index(child)
                    if i < len(cycle_stack) - 2: 
                      output_cycles.add(get_hashable_cycle(cycle_stack[i:]))
                
            except StopIteration:
                stack.pop()
                cycle_stack.pop()
    
    return [list(i) for i in output_cycles]



def draw_G(G, figsize=(6,6), sort_coord=False, pos_linear=False):
    fig, ax = plt.subplots(figsize=figsize)
    if sort_coord:
        sorted_nodes = sorted([k for k in G.nodes], key=lambda x: (chrom_sort_dict[G.nodes(data=True)[x]['CHR']], G.nodes(data=True)[x]['POS']))
    else:
        sorted_nodes = G.nodes
    if pos_linear:
        node_pos = {x:(G.nodes[x]['POS'], 1) for x in sorted_nodes}
        node_labels = None
    else:
        node_pos = nx.circular_layout(sorted_nodes)
        node_labels = {k:G.nodes(data=True)[k]['coord'] for k in sorted_nodes}
    node_colors = ['blue' if d['5_open'] else 'green' for n, d in G.nodes(data=True)]
    seg_list = [(u,v) for u,v,e in G.edges(data=True) if e['segment']]
    read_list = [(u,v) for u,v,e in G.edges(data=True) if e['read']]
    read_colors = [sv_hue_dict[e['SVTYPE']] for u,v,e in G.edges(data=True) if e['read']]
    
    nx.draw_networkx_nodes(G, node_pos, ax=ax, node_color=node_colors, node_size=100, alpha=1)
    nx.draw_networkx_labels(G, node_pos, ax=ax, labels=node_labels)
    nx.draw_networkx_edges(G, node_pos, ax=ax, edgelist=seg_list, edge_color='black')
    for read, read_color in zip(read_list, read_colors):
        arrowprop=dict(arrowstyle="-", linestyle='-', color=read_color, connectionstyle="arc3,rad=-0.3")
        ax.annotate("", xy=node_pos[read[0]], xytext=node_pos[read[1]], arrowprops=arrowprop)

