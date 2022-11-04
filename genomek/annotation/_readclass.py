import re
import itertools
import pysam
import numpy as np
from pysam.libcalignedsegment import AlignedSegment 
from ..tools import print_err
from ..tools._gadgets import satellite_path_37 as satellite_path

CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
CIGAR2CODE = dict([y, x] for x, y in enumerate("MIDNSHP=XB"))

def cigarstring_to_cigartuples(cigarstring):
    parts = CIGAR_REGEX.findall(cigarstring)
    return [(CIGAR2CODE[y], int(x)) for x,y in parts]


def cigartuples_to_number(cigar_tuples, start=True, read_length=151):
    if cigar_tuples is None:
        return 0
    cigar_list = list(zip(*cigar_tuples))
    try:
        index = cigar_list[0].index(0)
        if start:
            return sum(cigar_list[1][:index])
        else:
            return sum(cigar_list[1][index:])
    except:
        if start:
            return 0
        else:
            return read_length


class variant():
    def __init__(self, CHROM: str, POS: int, REF: str, ALT: str):
        '''the variants must be left aligned and normalized'''
        self.CHROM = CHROM
        self.POS = POS
        self.REF = REF
        self.ALT = ALT
        self.mttype = self.get_mttype()
        self.var_seq = self.get_variant_sequence()
        self.start, self.end = self.get_zero_coordinate()
        
    def __repr__(self):
        return f"CHROM:{self.CHROM} POS:{self.POS} REF:{self.REF} ALT:{self.ALT}"
        
    def get_mttype(self) -> str:
        ref = self.REF
        alt = self.ALT
        if len(ref) == len(alt):
            return 'snv' if len(ref) == 1 else 'mnv'
        else:
            if len(ref) < len(alt):
                return 'ins'
            elif len(ref) > len(alt):
                return 'del'
            else:
                raise ValueError('cannot specify mttype')
                # return 'cindel' #complex indel

    def get_variant_sequence(self) -> str:
        if self.mttype == 'snv' or self.mttype == 'mnv':
            return self.ALT
        elif self.mttype == 'ins':
            return self.ALT[len(self.REF):]
        elif self.mttype == 'del':
            return self.REF[len(self.ALT):]

    def get_zero_coordinate(self):
        # It is important to remember zero coordinates are used for fetching pysam
        # pysam uses 0-baed half-open coordinate system
        # for snv and mnv, I'm going to fetch reads precisely overlapping the mismatching regions
        # for indel, I'm going to fetch reads extending +1 base of that indel (both 5' and 3')
        if self.mttype == 'snv':
            start, end = self.POS - 1, self.POS
        elif self.mttype == 'mnv':
            start, end = self.POS - 1, self.POS - 1 + len(self.REF) 
        elif self.mttype == 'ins':
            start, end = self.POS - 1, self.POS + 1 #reads extending 0-based coordinates POS-1 and POS are fetched (0-based POS+1 is half-open and not included)
        elif self.mttype == 'del':
            start, end = self.POS - 1, self.POS + len(self.REF) 
        else:
            raise ValueError('mttype is not one of snv mnv ins del')
        return max(start, 0), end  # 0 based coordinates are used for fetching pysam alignment segments
            
    def fetch_reads_pysam(self, bam: pysam.AlignmentFile):
        reads = bam.fetch(contig=self.CHROM, start=self.start, end=self.end)
        return reads

    def fetch_reads(self, bam: pysam.AlignmentFile):
        for read in bam.fetch(contig=self.CHROM, start=self.start, end=self.end):
            yield readclass(read, self)


class readclass():
    def __init__(self, read, variant):
        self.read = read
        self.query_length = read.query_length
        self.read_aligned_pairs = read.get_aligned_pairs()
        self.cigarstats = read.get_cigar_stats()[0]
        for attribute in ['mttype', 'start', 'end', 'var_seq', 'REF', 'ALT']:
            setattr(self, attribute, getattr(variant, attribute))
        self.read_start = self.genome_pos_to_read_pos(self.start, start=True)
        self.read_end = self.genome_pos_to_read_pos(self.end, start=False)

    def genome_pos_to_read_pos(self, genome_pos, start):
        ''' 
        genome pos and read pos are both 0-based 
        start: boolean indicating whether the search is for the read_start or the read_end
        read_pos can be None or -1 (None indicates 'other supporting'  while -1 indicates non-informative)
        '''
        if start:
            for x in self.read_aligned_pairs:
                if x[1] is not None and genome_pos == x[1]:
                    read_pos = x[0]
                    return read_pos
        else:
            for x in reversed(self.read_aligned_pairs):
                # Because the coordinate is half-open we need to query genome pos -1 and then add +1 for read_pos
                if x[1] is not None and genome_pos - 1 == x[1]:
                    read_pos = x[0] + 1 if x[0] is not None else x[0]
                    return read_pos # max read_pos is 151
        # If looping over read aligned pairs fail (genome_pos is out of range in respect to read_aligned_pairs)
        return -1 # This read is uninformative (read is unmapped, read does not fully cover fetched region)
            
    def get_read_type(self):
        '''
        returns read_type and mean_BQ
        read_type:
          0: reference
          1: alt
          2: other
          3: non-informative
          3: deletion in variant site (only for snv)
        '''
        if self.read_start == -1 or self.read_end == -1:
            self.read_type = 3
            return 3
        if self.read_start is None or self.read_end is None:
            self.read_type = 2
            return 2
        if self.mttype in ['snv', 'mnv']:
            # read_sliced_sequence = "".join([self.read.query_sequence[i[0]] for i in self.read_aligned_pairs if i[1] in range(self.start, self.end) and i[0] is not None])
            read_sliced_sequence = self.read.query_sequence[self.read_start:self.read_end]
            if read_sliced_sequence == self.REF:
                read_type = 0
            elif read_sliced_sequence == self.var_seq:
                read_type = 1
            # elif len(read_sliced_sequence) != len(self.REF):
            #     read_type = 3 # indel in variant site
            else:
                read_type = 2
        elif self.mttype == 'ins':
            read_sliced_sequence = self.read.query_sequence[self.read_start+len(self.REF):self.read_end-1]
            if read_sliced_sequence == self.var_seq:
                read_type = 1
            elif read_sliced_sequence == '':
                read_type = 0
            else:
                read_type = 2
        elif self.mttype == 'del':
            read_sliced_sequence = self.read.query_sequence[self.read_start+len(self.ALT):self.read_end-1]
            if len(read_sliced_sequence) == 0:
                read_type = 1
            elif len(read_sliced_sequence) == len(self.var_seq):
                read_type = 0
            else:
                read_type = 2
        else:
            raise ValueError("function: get_read_type error")
        self.read_type = read_type
        return read_type

    def get_read_type_and_BQ(self):
        '''
        returns read_type and mean_BQ
        read_type:
          0: reference
          1: alt
          2: other
          3: deletion in variant site (only for snv)
        '''
        if self.read_start == -1 or self.read_end == -1:
            self.read_type = 3
            return 3, 0
        if self.read_start is None or self.read_end is None:
            self.read_type = 2
            return 2, 0
        if self.mttype in ['snv', 'mnv']:
            # read_sliced_sequence = "".join([self.read.query_sequence[i[0]] for i in self.read_aligned_pairs if i[1] in range(self.start, self.end) and i[0] is not None])
            # read_sliced_mean_BQ = np.nan_to_num(np.mean([self.read.query_qualities[i[0]] for i in self.read_aligned_pairs if i[1] in range(self.start, self.end) and i[0] is not None], dtype=np.float16))
            read_sliced_sequence = self.read.query_sequence[self.read_start:self.read_end]
            read_sliced_mean_BQ = np.nan_to_num(np.mean(self.read.query_qualities[self.read_start:self.read_end], dtype=np.float16))
            if read_sliced_sequence == self.REF:
                read_type = 0
            elif read_sliced_sequence == self.var_seq:
                read_type = 1
            else:
                read_type = 2
        elif self.mttype == 'ins':
            read_sliced_sequence = self.read.query_sequence[self.read_start+len(self.REF):self.read_end-1]
            read_sliced_mean_BQ = np.nan_to_num(np.mean(self.read.query_qualities[self.read_start+len(self.REF):self.read_end-1], dtype=np.float16))
            if read_sliced_sequence == self.var_seq:
                read_type = 1
            elif read_sliced_sequence == '':
                read_type = 0
            else:
                read_type = 2
        elif self.mttype == 'del':
            read_sliced_sequence = self.read.query_sequence[self.read_start+len(self.ALT):self.read_end-1]
            read_sliced_mean_BQ = np.nan_to_num(np.mean(self.read.query_qualities[self.read_start+len(self.ALT):self.read_end-1], dtype=np.float16))
            if len(read_sliced_sequence) == 0:
                read_type = 1
            elif len(read_sliced_sequence) == len(self.var_seq):
                read_type = 0
            else:
                read_type = 2
        else:
            raise ValueError("function: get_read_type error")
        self.read_type = read_type
        return read_type, read_sliced_mean_BQ

    def get_read_basic_filter_flag(self):
        if self.read.is_duplicate or self.read.is_qcfail or self.read.is_unmapped: #or self.read.mate_is_unmapped : 
            return True
        return False

    def get_read_name_dup_flag(self):
        if self.read.is_secondary or self.read.is_supplementary:
            return True
        return False

    def get_read_sv_flag(self):
        if not self.read.is_proper_pair or (self.read.reference_name != self.read.next_reference_name) or self.read.has_tag("SA"):
            return True
        return False

    def get_read_proximal_flag(self):
        if self.read.is_unmapped or self.read.mate_is_unmapped:
            return False
        if not self.read.has_tag("MC"): # if mate does not have MC tag
            return False
        self_cigarstring, mate_cigarstring = self.read.cigarstring, self.read.get_tag("MC")
        self_start_clip_length, mate_start_clip_length = cigartuples_to_number(cigarstring_to_cigartuples(self_cigarstring)), cigartuples_to_number(cigarstring_to_cigartuples(mate_cigarstring))
        proximal_flag = self.read.reference_start - self_start_clip_length < self.read.next_reference_start - mate_start_clip_length
        return proximal_flag

    def get_read_MM(self):
        MM = self.read.get_tag("NM") - self.cigarstats[1] - self.cigarstats[2] if self.read.has_tag("NM") else 0
        if self.read_type == 1 and (self.mttype == 'snv' or self.mttype == 'mnv'):
            MM -= len(self.REF)
        return MM

    def get_five_prime_distance(self):
        if self.read_start is None or self.read_end is None or self.read_start == 0 or self.read_end == 0:
            return 0
        if self.read.is_reverse:
            distance = self.query_length - self.read_end
        else:
            distance = self.read_start
        return distance


def readInfoBasic(readclass):
    read_type = readclass.get_read_type()
    read_filter_flag = readclass.get_read_basic_filter_flag()
    read_name_dup_flag = readclass.get_read_name_dup_flag()
    return read_type, read_filter_flag, read_name_dup_flag


def readInfoBasic2(readclass):
    read_type = readclass.get_read_type()
    read_filter_flag = readclass.get_read_basic_filter_flag()
    read_name_dup_flag = readclass.get_read_name_dup_flag()
    read_sv_flag = readclass.get_read_sv_flag()
    read_MQ = readclass.read.mapping_quality
    return read_type, read_filter_flag, read_name_dup_flag, read_sv_flag, read_MQ


def readInfoFull(readclass: readclass):
    read_type, read_mean_BQ = readclass.get_read_type_and_BQ()
    read_filter_flag = readclass.get_read_basic_filter_flag()
    read_name_dup_flag = readclass.get_read_name_dup_flag()
    read_sv_flag = readclass.get_read_sv_flag()
    read_proximal_flag = readclass.get_read_proximal_flag()
    read_reverse_flag = readclass.read.is_reverse
    read_read1_flag = readclass.read.is_read1
    read_MQ = readclass.read.mapping_quality
    read_len_ins = readclass.cigarstats[1]
    read_len_del = readclass.cigarstats[2]
    read_len_clip = readclass.cigarstats[4]
    read_template_length = abs(readclass.read.template_length)
    read_MM = readclass.get_read_MM()
    read_five_prime_distance = readclass.get_five_prime_distance()
    return read_type, read_mean_BQ, read_filter_flag, read_name_dup_flag, read_sv_flag, read_proximal_flag, read_reverse_flag, read_read1_flag, read_MQ, read_len_ins, read_len_del, read_len_clip, read_template_length, read_MM, read_five_prime_distance


def get_read_names(reads, reads_limit=10000):
    # _, read_names = np.unique(np.array([read.query_name for read in reads]), return_inverse=True)
    _, read_names = np.unique(np.array([read_enumerated[1].query_name for read_enumerated in itertools.takewhile(lambda x: x[0] < reads_limit, enumerate(reads))]), return_inverse=True)
    return read_names


def readsDataBasic(variant, bam, reads_limit=10000):
    reads = variant.fetch_reads_pysam(bam)
    index = get_read_names(reads, reads_limit=reads_limit)
    n_reads = len(index)
    reads = variant.fetch_reads(bam)
    
    array_type = np.full(n_reads, 0, dtype=np.uint8)
    array_basic_filter = np.full(n_reads, False)
    array_name_dup_filter = np.full(n_reads, False)
    for i, read in enumerate(reads):
        if i >= reads_limit: break
        array_type[i], array_basic_filter[i], array_name_dup_filter[i]= readInfoBasic(read)

    data_dict = {
    'read_index':index,
    'type':array_type,
    'basic_filter':array_basic_filter ,
    'name_dup_filter':array_name_dup_filter,
    }
    return data_dict


def readsDataBasic2(variant, bam, reads_limit=10000):
    reads = variant.fetch_reads_pysam(bam)
    index = get_read_names(reads, reads_limit=reads_limit)
    n_reads = len(index)
    reads = variant.fetch_reads(bam)
    
    array_type = np.full(n_reads, 0, dtype=np.uint8)
    array_basic_filter = np.full(n_reads, False)
    array_name_dup_filter = np.full(n_reads, False)
    array_sv_flag = np.full(n_reads, False)
    array_MQ = np.full(n_reads, 0, dtype=np.uint8)
    for i, read in enumerate(reads):
        if i >= reads_limit: break
        array_type[i], array_basic_filter[i], array_name_dup_filter[i], array_sv_flag[i], array_MQ[i]= readInfoBasic2(read)

    data_dict = {
    'read_index':index,
    'type':array_type,
    'basic_filter':array_basic_filter ,
    'name_dup_filter':array_name_dup_filter,
    'sv_flag':array_sv_flag ,
    'MQ':array_MQ 
    }
    return data_dict


def readsDataFull(variant, bam, reads_limit=10000):
    reads = variant.fetch_reads_pysam(bam)
    index = get_read_names(reads, reads_limit=reads_limit)
    n_reads = len(index)
    reads = variant.fetch_reads(bam)

    array_type = np.full(n_reads, 0, dtype=np.uint8)
    array_BQ = np.full(n_reads, 0, dtype=np.float16)
    array_basic_filter = np.full(n_reads, False)
    array_name_dup_filter = np.full(n_reads, False)
    array_sv_flag = np.full(n_reads, False)
    array_proximal_flag = np.full(n_reads, False)
    array_reverse_flag = np.full(n_reads, False)
    array_read1_flag = np.full(n_reads, False)
    array_MQ = np.full(n_reads, 0, dtype=np.uint8)
    array_len_ins = np.full(n_reads, 0, dtype=np.uint8)
    array_len_del = np.full(n_reads, 0, dtype=np.uint8)
    array_len_clip = np.full(n_reads, 0, dtype=np.uint8)
    array_template_length = np.full(n_reads, 0, dtype=np.int64)
    array_MM = np.full(n_reads, 0, dtype=np.uint8)
    array_five_prime_distance = np.full(n_reads, 0, dtype=np.uint8)
    for i, read in enumerate(reads):
        if i >= reads_limit: break
        array_type[i], array_BQ[i], array_basic_filter[i], array_name_dup_filter[i], array_sv_flag[i], array_proximal_flag[i], array_reverse_flag[i], array_read1_flag[i], array_MQ[i], array_len_ins[i], array_len_del[i], array_len_clip[i], array_template_length[i], array_MM[i], array_five_prime_distance[i] = readInfoFull(read)

    data_dict = {
    'read_index':index,
    'type':array_type,
    'BQ':array_BQ,
    'basic_filter':array_basic_filter ,
    'name_dup_filter':array_name_dup_filter,
    'sv_flag':array_sv_flag ,
    'proximal_flag':array_proximal_flag ,
    'reverse_flag':array_reverse_flag ,
    'read1_flag':array_read1_flag ,
    'MQ':array_MQ ,
    'len_ins':array_len_ins ,
    'len_del':array_len_del ,
    'len_clip':array_len_clip ,
    'template_length':array_template_length ,
    'MM':array_MM ,
    'five_prime_distance':array_five_prime_distance
    }
    return data_dict


