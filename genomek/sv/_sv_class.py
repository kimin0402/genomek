from ._sv_annotation import bam_to_readdict
from ..vcf._filter import variant_in_bed
from ..tools._gadgets import chrom_sort_dict_37 as chrom_sort_dict
from ..tools._gadgets import mask_path_37 as mask_path
from ..tools._gadgets import satellite_path_37 as satellite_path
from ..annotation._readclass import cigarstring_to_cigartuples
from collections import Counter
import pyranges as pr
p_satellite = pr.read_bed(satellite_path).set_union(pr.read_bed(mask_path))

def read_fails_qc(read):
    if read.is_duplicate or read.is_qcfail or read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < 30: 
        return True
    return False

def adjust_bp_from_cigar(pos, cigarstring):
    '''get read.reference_end for sa read if clip is at 3prime end (if the clip is at the right side)''' 
    cigartuples = cigarstring_to_cigartuples(cigarstring)
    adjusted_pos = pos + sum(y for x,y in cigartuples if x in [0,2,3,7,8]) -1 #from sam format cigarstring add op that consumes reference. -1 for excluding the first match
    return int(adjusted_pos)

def get_candidate_bp_list_sc(chr1, pos1, c1, bam, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
    '''
    for chr1 and pos1, identify candidate bp1 and bp2 from soft clip and supplementary alignemtns
    return two counter objects
    '''
    c_dict = {'3':-1, '5':0}
    if c1 == '3':
        fetch_lt = fetch_rearrangement_side # if c1 is 3' clip is at 3' and 5' direction (left) is where the rearranged chromosome exists
        fetch_rt = fetch_non_rearrangement_side
    else:
        fetch_lt = fetch_non_rearrangement_side
        fetch_rt = fetch_rearrangement_side
    bp_list_from_sc = []
    start = max(0, pos1-1-fetch_lt)
    for read in bam.fetch(contig=chr1, start=start, end=pos1+fetch_rt):
        if read_fails_qc(read):
            continue
        if read.cigartuples[c_dict[c1]][0] == 4 or read.cigartuples[c_dict[c1]][0] == 5: # search clip depending on c
            if c_dict[c1]: # if the clip is right
                bp_list_from_sc.append((read.reference_name, read.reference_end, c1))
            else: # if the clip is left
                bp_list_from_sc.append((read.reference_name, read.reference_start+1, c1))
    return Counter(bp_list_from_sc)

def get_candidate_bp_list_sa(chr1, pos1, c1, chr2, pos2, c2, bam, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
    '''
    for chr1 and pos1, identify candidate bp1 and bp2 from soft clip and supplementary alignemtns
    return two counter objects
    '''
    c_dict = {'3':-1, '5':0}
    if c1 == '3':
        fetch_lt = fetch_rearrangement_side # if c1 is 3' clip is at 3' and 5' direction (left) is where the rearranged chromosome exists
        fetch_rt = fetch_non_rearrangement_side
    else:
        fetch_lt = fetch_non_rearrangement_side
        fetch_rt = fetch_rearrangement_side
    bp_list_from_sa = []
    start = max(0, pos1-1-fetch_lt)
    for read in bam.fetch(contig=chr1, start=start, end=pos1+fetch_rt):
        if read_fails_qc(read):
            continue
        if read.has_tag('SA'):
            chrom_sa, pos_sa, strand_sa, cigar_sa, mq_sa, nm_sa, *_ = read.get_tag('SA').split(',')
            pos_sa, mq_sa = int(pos_sa), int(mq_sa)
            if mq_sa >= 30 and chrom_sa == chr2: # simple filter to check sa
                if c_dict[c2]:
                    bp_list_from_sa.append((chrom_sa, adjust_bp_from_cigar(pos_sa, cigar_sa), c2))
                else:
                    bp_list_from_sa.append((chrom_sa, pos_sa, c2))
    return Counter(bp_list_from_sa)

def get_candidate_bp_list_dc(chr1, pos1, c1, chr2, pos2, c2, bams: list, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
    def get_readdict(chr1, pos1, c1, bam, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
        if c1 == '3':
            fetch_lt = fetch_rearrangement_side # if c1 is 3' clip is at 3' and 5' direction (left) is where the rearranged chromosome exists
            fetch_rt = fetch_non_rearrangement_side
        else:
            fetch_lt = fetch_non_rearrangement_side
            fetch_rt = fetch_rearrangement_side
        readdict = bam_to_readdict(bam, chr1, pos1, left_margin=fetch_lt, right_margin=fetch_rt)
        return readdict

    def find_bp_from_dc_reads(readdict1, c1, readdict2, c2):
        '''
        find discordant fragments
        '''
        sv_candidates = set()
        ct_list = [{'5': True, '3': False}[x] for x in (c1, c2)]
        pos1_dc_list = list()
        pos2_dc_list = list()
        for read_name in set(readdict1.keys()).intersection(set(readdict2.keys())):
            if len(frag_list := set(readdict1[read_name]).union(set(readdict2[read_name]))) == 2:
                frag_list = sorted(list(frag_list), key=lambda x: (chrom_sort_dict.get(x.reference_name, 25), x.reference_start))
                if c1 == '5':
                    pos1_dc_list.append((frag_list[0].reference_start))
                else:
                    pos1_dc_list.append((frag_list[0].reference_end))
                if c2 == '5':
                    pos2_dc_list.append((frag_list[1].reference_start))
                else:
                    pos2_dc_list.append((frag_list[1].reference_end))
        return pos1_dc_list, pos2_dc_list

    for i, bam in enumerate(bams):
        if i == 0:
            bp1_readdict = get_readdict(chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            bp2_readdict = get_readdict(chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
        else:
            bp1_readdict = {**bp1_readdict, **get_readdict(chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)}
            bp2_readdict = {**bp2_readdict, **get_readdict(chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)}
    bp1, bp2, = find_bp_from_dc_reads(bp1_readdict, c1, bp2_readdict, c2)
    if len(bp1):
        if c1 == '5':
            bp1_return = min(bp1)
        else:
            bp1_return = max(bp1)
    else:
        bp1_return = None
    if len(bp2):
        if c2 == '5':
            bp2_return = min(bp2)
        else:
            bp2_return = max(bp2)
    else:
        bp2_return = None
    return (chr1, bp1_return, c1), (chr2, bp2_return, c2)


def adjust_bp_rowop(chr1, pos1, c1, chr2, pos2, c2, bams, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
    for i, bam in enumerate(bams):
        if i == 0:
            sc1_counter = get_candidate_bp_list_sc(chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sc2_counter = get_candidate_bp_list_sc(chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sa2_counter = get_candidate_bp_list_sa(chr1, pos1, c1, chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sa1_counter = get_candidate_bp_list_sa(chr2, pos2, c2, chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
        else:
            sc1_counter += get_candidate_bp_list_sc(chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sc2_counter += get_candidate_bp_list_sc(chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sa2_counter += get_candidate_bp_list_sa(chr1, pos1, c1, chr2, pos2, c2, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
            sa1_counter += get_candidate_bp_list_sa(chr2, pos2, c2, chr1, pos1, c1, bam, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)
    (chr1_dc, pos1_dc, c1_dc), (chr2_dc, pos2_dc, c2_dc) = get_candidate_bp_list_dc(chr1, pos1, c1, chr2, pos2, c2, bams, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side)

    scsa1_counter = sc1_counter + Counter({k:d*3 for k,d in sa1_counter.items()})
    scsa2_counter = sc2_counter + Counter({k:d*3 for k,d in sa2_counter.items()})


    if len(scsa1_counter):
        (chr1_sc, pos1_sc, c1_sc), n1_sc = [(k,v) for k,v in sorted(scsa1_counter.items(), key=lambda x: x[1], reverse=True)][0]
    else:
        (chr1_sc, pos1_sc, c1_sc), n1_sc = (None, None, c1), 0
    if len(scsa2_counter):
        (chr2_sc, pos2_sc, c2_sc), n2_sc = [(k,v) for k,v in sorted(scsa2_counter.items(), key=lambda x: x[1], reverse=True)][0]
    else:
        (chr2_sc, pos2_sc, c2_sc), n2_sc = (None, None, c2), 0
    
    if n1_sc > 4:
        chr1, pos1 = chr1_sc, pos1_sc
    else:
        chr1, pos1 = chr1_dc, pos1_dc
    if n2_sc > 4:
        chr2, pos2 = chr2_sc, pos2_sc
    else:
        chr2, pos2 = chr2_dc, pos2_dc
    # Last sanity check
    if chr1 is not None and chr2 is not None and pos1 is not None and pos2 is not None and chr1 == chr2 and pos1 >= pos2:
        chr1, chr2, pos1, pos2 = None, None, None, None
    return chr1, pos1, c1, n1_sc, chr2, pos2, c2, n2_sc

    # if len(sa1_counter):
    #     (chr1_sa, pos1_sa, c1_sa), n1_sa = [(k,v) for k,v in sorted(sa1_counter.items(), key=lambda x: x[1], reverse=True)][0]
    # else:
    #     (chr1_sa, pos1_sa, c1_sa), n1_sa = (None, None, None), 0
    # if len(sa2_counter):
    #     (chr2_sa, pos2_sa, c2_sa), n2_sa = [(k,v) for k,v in sorted(sa2_counter.items(), key=lambda x: x[1], reverse=True)][0]
    # else:
    #     (chr2_sa, pos2_sa, c2_sa), n2_sa = (None, None, None), 0
    # return chr1_sc, pos1_sc, n1_sc, chr2_sc, pos2_sc, n2_sc, chr1_sa, pos1_sa, n1_sa, chr2_sa, pos2_sa, n2_sa, chr1_dc, pos1_dc, chr2_dc, pos2_dc


def df_adjust_bp(df, bams: list, ref_key, fetch_rearrangement_side=0, fetch_non_rearrangement_side=0):
    filter_bool1 = variant_in_bed(df, p_satellite, chrom_key='CHR1', pos_key='POS1', ref_key=ref_key)
    filter_bool2 = variant_in_bed(df, p_satellite, chrom_key='CHR2', pos_key='POS2', ref_key=ref_key)
    df = df[(~filter_bool1 & ~filter_bool2) & ((df['CHR1'] != df['CHR2']) | (abs(df['POS1'] - df['POS2']) > 1000))]
    df_result = df.apply(lambda row: adjust_bp_rowop(row['CHR1'], row['POS1'], row['CT'].split('to')[0], \
                                                     row['CHR2'], row['POS2'], row['CT'].split('to')[1], \
                                                     bams, fetch_rearrangement_side=fetch_rearrangement_side, fetch_non_rearrangement_side=fetch_non_rearrangement_side), 
                         axis=1, result_type='expand')
    return df_result
