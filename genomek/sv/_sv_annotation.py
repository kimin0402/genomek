import pyranges as pr
from collections import defaultdict
from datetime import datetime
from ..tools import print_err
from ..tools import get_vaf
from ..tools._gadgets import chromosomes_37 as chromosomes
from ..tools._gadgets import chrom_sort_dict_37 as chrom_sort_dict
from ..tools._gadgets import satellite_path_37 as satellite_path

satellite_pr = pr.read_bed(satellite_path)
bad_pr = pr.PyRanges(chromosomes='2', starts=[33141145], ends=[33141699])
filter_pr = satellite_pr.set_union(bad_pr)

def read_filter_basic(read):
    '''
    basic filter that returns True for unqualified reads (reads that break the code when searching for mates)
    '''
    if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_unmapped or read.mate_is_unmapped:
        return True
    return False

def find_mate(read, bam):
    '''
    for a given read, returns the mate without disrupting the pysam bam object
    '''
    for mate_read in bam.fetch(read.next_reference_name, max(read.next_reference_start-2, 0), read.next_reference_start+2):
        if read_filter_basic(mate_read):
            continue
        if mate_read.query_name == read.query_name:
            if mate_read.is_read1 != read.is_read1 or mate_read.reference_start != read.reference_start or mate_read.reference_end != read.reference_end or mate_read.reference_name != read.reference_name:
                return mate_read
    return None

def add_mate_to_readdict(readdict, bam):
    '''
    readdict is a dictionary with read names as a key and list of AlignmentObject as values
    this function first checks if there is only one AlignmentObject in the value, and appends the mate to it
    AlignmentObjects are sorted according to the genome coordinate
    '''
    for read_name, read_list in readdict.items():
        if len(read_list) == 1:
            original_read = read_list[0]
            if read_filter_basic(original_read):
                continue
            mate_read = find_mate(original_read, bam)
            if not mate_read:
                continue
            read_list.append(mate_read)
            readdict[read_name] = sorted(read_list, key=lambda x: (chrom_sort_dict.get(x.reference_name, 25), x.reference_start)) # if mate is mapped at a special contig such as GL000191.1, chrom_sort_dict will return integer 25
    return readdict

def sort_readdict(readdict):
    for read_name, read_list in readdict.items():
        readdict[read_name] = sorted(read_list, key=lambda x: (chrom_sort_dict.get(x.reference_name, 25), x.reference_start)) # if mate is mapped at a special contig such as GL000191.1, chrom_sort_dict will return integer 25
    return readdict

def bam_to_readdict(bam, chr1, pos1, left_margin, right_margin, readdict_limit=20000):
    '''
    given a pysam alignmentfile and genome coordinate, returns a readdict
    readdict is a dictionary with read names as a key and list of AlignmentObject as values
    '''
    readdict = defaultdict(list)
    start = max(0, pos1-left_margin)
    end = pos1+right_margin
    reads = bam.fetch(contig=chr1, start=start, end=end)
    for read in (x for x in reads if not read_filter_basic(x) and x.reference_start in range(start, end+1) and x.reference_end in range(start, end+1)):
        readdict[read.query_name].append(read)
        if len(readdict) > readdict_limit:
            break        
    else:
        readdict = sort_readdict(readdict)
        return readdict
    print_err(f"# readdict exceeded {readdict_limit=} keys") # here we can skip add_mate_to_readdict using a for-else block
    readdict = sort_readdict(readdict)
    return readdict

def find_ref_frag(readdict, chr1, pos1, isize=None, ref_padding=3):
    '''
    among all sets of fragments in readdict, find fragments that spans through a certain breakpoint
    fragments with given isize is filtered out, but if isize is not specified, fragments with any size is included in the result
    '''
    ref_frag_set = set()
    for read_name, read_list in readdict.items():
        if len(read_list) != 2:
            continue
        if read_list[0].is_proper_pair and read_list[0].reference_name == read_list[1].reference_name and read_list[0].reference_name == chr1:
            if read_list[0].reference_start <= pos1-ref_padding  and pos1+ref_padding <= read_list[1].reference_end:
                if isize and abs(read_list[0].template_length) <= isize:
                    ref_frag_set.add(read_name)
                    continue
                else:
                    ref_frag_set.add(read_name)
    return ref_frag_set

def find_sa_readname(readdict, pos_of_readdict, chr_other, pos_other, isize):
    '''
    find fragments with SA tag
    '''
    sv_candidates = set()
    for read_name, read_list in readdict.items():
        for read in read_list:
            if read.reference_start <= pos_of_readdict and pos_of_readdict <= read.reference_end and read.has_tag("SA"):
                chr_sa, pos_sa = read.get_tag("SA").split(",")[0:2]
                if chr_other == chr_sa and int(pos_sa) in range(pos_other-isize, pos_other+isize):
                    sv_candidates.add(read_name)
    return sv_candidates

def find_discordant_readname(readdict1, pos_of_bp1, readdict2, pos_of_bp2, ct, discordant_padding=100):
    '''
    find discordant fragments
    '''
    sv_candidates = set()
    ct_list = [{'5': True, '3': False}[x] for x in ct.split('to')]
    for read_name in set(readdict1.keys()).intersection(set(readdict2.keys())):
        if len(frag_list := set(readdict1[read_name]).union(set(readdict2[read_name]))) == 2:
            frag_list = sorted(list(frag_list), key=lambda x: (chrom_sort_dict.get(x.reference_name, 25), x.reference_start))
            if frag_list[0].is_reverse == ct_list[0] and frag_list[1].is_reverse == ct_list[1]:
                if ct_list[0]: # if first bp is 5' test whether discordant read is located right side of bp
                    bp_1_is_correct = pos_of_bp1 - discordant_padding <= frag_list[0].reference_start
                else: 
                    bp_1_is_correct = frag_list[0].reference_end <=  pos_of_bp1 + discordant_padding 
                if ct_list[1]: # if second bp is 5' 
                    bp_2_is_correct = pos_of_bp2 - discordant_padding <= frag_list[1].reference_start
                else:
                    bp_2_is_correct = frag_list[1].reference_end <=  pos_of_bp2 + discordant_padding 
                if bp_1_is_correct and bp_2_is_correct:
                    sv_candidates.add(read_name)
    return sv_candidates

def calculate_sv_vaf(chr1, pos1, chr2, pos2, ct, isize, bam):
    '''
    bam: pysam.AlignmentFile
    '''
    if chr1 == chr2 and pos2 in range(pos1-isize, pos1+isize):
        print_err(chr1, pos1, chr2, pos2, 'WARNING:pos1_and_pos2_overlaps', sep='\t')
        is_overlap = True
    else:
        is_overlap = False
    if len(pr.PyRanges(chromosomes=[chr1, chr2], starts=[pos1, pos2], ends=[pos1, pos2]).intersect(filter_pr, how='containment')) > 0:
        print_err(chr1, pos1, chr2, pos2, 'satellite', sep='\t')
        return None, None, is_overlap
    
    left_margin, right_margin = isize, isize
    start_time = datetime.now()
    bp_1 = bam_to_readdict(bam=bam, chr1=chr1, pos1=pos1, left_margin=left_margin, right_margin=right_margin)
    bp_1_time = datetime.now() - start_time

    start_time = datetime.now()
    bp_2 = bam_to_readdict(bam=bam, chr1=chr2, pos1=pos2, left_margin=right_margin, right_margin=left_margin)
    bp_2_time = datetime.now() - start_time

    start_time = datetime.now()
    bp_1_ref = find_ref_frag(bp_1, chr1, pos1, isize=isize)
    bp_2_ref = find_ref_frag(bp_2, chr2, pos2, isize=isize)
    sa_reads = find_sa_readname(bp_1, pos1, chr2, pos2, isize=isize).union(find_sa_readname(bp_2, pos2, chr1, pos1, isize=isize))
    dc_reads = find_discordant_readname(bp_1, pos1, bp_2, pos2, ct)
    sv_reads = sa_reads.union(dc_reads)
    num_sv_reads = len(sv_reads)
    vaf1 = get_vaf(num_sv_reads, len(bp_1_ref))
    vaf2 = get_vaf(num_sv_reads, len(bp_2_ref))
    other_time = datetime.now() - start_time

    print_err(chr1, pos1, chr2, pos2, bp_1_time, bp_2_time, other_time, sep='\t')
    return vaf1, vaf2, is_overlap



## For gremlin processing
def coord_to_svcols(string):
    '''
    10_100887832_13_19825866_3to5 -> 10, 100887832, 13, 19825866, 3to5
    bp sorting implemented
    
    '''
    chr1, pos1, chr2, pos2, ct = string.split("_")
    pos1, pos2 = int(pos1), int(pos2)
    if (chrom_sort_dict[chr1] > chrom_sort_dict[chr2]) or (chrom_sort_dict[chr1] == chrom_sort_dict[chr2] and pos1>pos2):
        chr1, pos1, chr2, pos2 = chr2, pos2, chr1, pos1
        ct = 'to'.join(ct.split('to')[::-1])
    return chr1, pos1, chr2, pos2, ct