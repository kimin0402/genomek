import numpy as np
import pysam
import pyranges as pr
from collections import Counter
from collections import defaultdict
from ..tools import print_err

satellite_pr_dict = {}
satellite_pr_dict['37'] = pr.read_bed('/home/users/pjh/References/ucsc_RepeatMasker_files/hg19/custom_files/satellite_hg19_rename_slop_subt_sort.bed.gz')
satellite_pr_dict['38'] = pr.read_bed('/home/users/pjh/References/ucsc_RepeatMasker_files/hg38/satellite_hg38_slop_subt_sort.bed.gz')
stackDict_path = '/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt'
stackDict = {k:v for k, v in np.genfromtxt('/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt', dtype=int)}


def calculate_mean_med(x):
    if x == None or x == []:
        return None, None
    else:
        return np.mean(x), int(np.median(x))


def get_mttype(ref: str, alt:str) -> str:
    if len(ref) == len(alt):
        return 'snv' if len(ref) == 1 else 'mnv'
    else:
        if len(ref) != 1 and len(alt) != 1:
            return 'cindel'
        elif len(ref) == 1:
            return 'ins'
        else:
            return 'del'


def read_filter_basic(read):
    if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_unmapped or read.mate_is_unmapped:
        return True
    return False

def find_mate(read, bam):
    for mate_read in bam.fetch(read.next_reference_name, read.next_reference_start, read.next_reference_start + 1):
        if mate_read.query_name == read.query_name:
            if mate_read.is_read1 != read.is_read1 or mate_read.reference_start != read.reference_start or mate_read.reference_end != read.reference_end or mate_read.reference_name != read.reference_name:
                return mate_read
    return None

def add_mate_to_readdict(readdict, bam):
    for read_name, read_list in readdict.items():
        if len(read_list) == 1:
            original_read = read_list[0]
            if read_filter_basic(original_read):
                continue
            mate_read = find_mate(original_read, bam)
            if read_filter_basic(mate_read):
                continue
            read_list.append(mate_read)
    return None   
    

class readplus():
    def __init__(self, read):
        self.read = read
        try:
            self.pairs = read.get_aligned_pairs(with_seq=True)
        except ValueError: # when MD tag is absent, get_aligned_pairs method fails with ValueError
            self.pairs = None
            return
            
        self.readclass = None
        self.is_irrelevant = False
        self.irrelevant_cause = list()
        self.is_lowqual = False
        self.lowqual_cause = list()
        self.is_SVread = False
        self.SVread_cause = list()
        self.cigarstats = read.get_cigar_stats()[0]
        self.lenIns = self.cigarstats[1]
        self.lenDel = self.cigarstats[2]
        self.lenClip = self.cigarstats[4]
        self.MQ = read.mapping_quality
        self.BQ = None
        self.margins = ( read.reference_start, read.reference_end )
        self.orientation = None


    def set_irrelevant_reads(self, start: int, end: int, flankLen: int=1):
        if not (start - self.read.reference_start >= flankLen and self.read.reference_end - end >= flankLen ):
            self.is_irrelevant = True
            self.irrelevant_cause.append('insufficient_flanking_sequences')
        elif not self.read.has_tag('MD'):
            self.is_irrelevant = True
            self.irrelevant_cause.append('absent_MD_tag')
        elif self.read.is_unmapped:
            self.is_irrelevant = True
            self.irrelevant_cause.append('flag_unmapped')
        elif self.read.is_duplicate:
            self.is_irrelevant = True
            self.irrelevant_cause.append('flag_duplicate')
        elif self.read.is_supplementary:
            self.is_irrelevant = True
            self.irrelevant_cause.append('flag_supplementary')
        elif self.read.is_secondary:
            self.is_irrelevant = True
            self.irrelevant_cause.append('flag_secondary')
        elif self.read.is_qcfail:
            self.is_irrelevant = True
            self.irrelevant_cause.append('flag_qcfail')
        else:
            pass


    def set_irrelevant_reads_by_repeat(self, relevant_repeats):
        for repeat in relevant_repeats:
            if not ( self.read.reference_start < repeat[0][0] and self.read.reference_end > repeat[0][1] ):
                self.is_irrelevant = True
                self.irrelevant_cause.append('insufficient_repeat_span')
                break


    def set_SVread(self, clip_threshold: int=0, template_length_threshold: int=2000):
        if self.lenClip > clip_threshold:
            self.is_SVread= True
            self.SVread_cause.append('clipping')
        if self.read.is_proper_pair == False:
            self.is_SVread= True
            self.SVread_cause.append('discordant pair')
        if self.read.reference_name != self.read.next_reference_name:
            self.is_SVread= True
            self.SVread_cause.append('Mate on difference chromosome')
        if abs(self.read.template_length) > template_length_threshold:
            self.is_SVread= True
            self.SVread_cause.append(f'template length > {template_length_threshold}')
        if self.read.has_tag('SA'):
            self.is_SVread= True
            self.SVread_cause.append('SA tag')


    def set_pairs_list(self, start: int, end: int, flankLen: int=0):
        if self.is_irrelevant:
            self.pairs_before_REF = None
            self.pairs_after_REF  = None
            self.pairs_within_REF = None
            self.query_pos_within_REF = None
            self.seq_within_REF = None
        else:
            idx_REFstart = list(zip(*self.pairs))[1].index(start)
            idx_nexttoREF = idx_REFstart + (end - start)
                
            try:
                self.pairs_after_REF  = self.pairs[idx_nexttoREF : idx_nexttoREF + flankLen]
            except:
                self.pairs_after_REF = []
            self.pairs_before_REF = self.pairs[idx_REFstart - flankLen : idx_REFstart]
            self.pairs_within_REF = self.pairs[idx_REFstart : idx_nexttoREF]
            self.query_pos_within_REF = [ x[0] for x in self.pairs_within_REF if x[0] != None ] 
                # 0-based
                # this can be an empty list
            self.seq_within_REF = ''.join([ self.read.query_sequence[idx] for idx in self.query_pos_within_REF ])
                # this can be an empty string


    def set_alleleClass(self, REF: str, ALT: str):
        # determine a read supports which allele ; executed only when readclass is neither -2 nor -1
        # 0: ref supporting ; 1: alt supporting ; 2: other supporting
        def check_flanking(pairs_before_REF, pairs_after_REF):
            # check if all flanking sequences are cigar M with identical bases
            for x in pairs_before_REF + pairs_after_REF:
                if x[0] == None or x[1] == None or x[2].islower():
                    return False
            return True

        if not check_flanking(self.pairs_before_REF, self.pairs_after_REF):
            #self.is_lowqual = True
            self.readclass = 2
        else:
            if self.seq_within_REF == REF:
                self.readclass = 0
            elif self.seq_within_REF == ALT:
                self.readclass = 1
            else:
                self.readclass = 2


    def set_MM_BQ(self, var_type: str):
        # NM tag value includes indels as well as base mismatch ; indel length is subtracted
        self.MM = self.read.get_tag('NM') - self.cigarstats[1] - self.cigarstats[2] if self.read.has_tag('NM') else None
        # values attributed to the mutation itself are subtracted
        if self.readclass == 1 and var_type == 'snv':
            self.MM -= 1
        if len(self.query_pos_within_REF) != 0:
            self.BQ = np.mean( [ self.read.query_qualities[idx] for idx in self.query_pos_within_REF ] )
        

    def set_Loca(self, REF: str, ALT: str, var_type: str):
        if self.readclass == 1:
            if var_type == 'ins':
                self.lenIns -= ( len(ALT) - len(REF) )
            elif var_type == 'del':
                self.lenDel -= ( len(REF) - len(ALT) )
            # calculating LocaLt/LocaRt (1-based, including softclip)
            # LocaLt: index of the target base counted from the left (index of the leftmost base is 1)
                # target base: snv -> the substituted base ; indel -> the base corresponding to REF
            # LocaRt: counting from the rightmost base (index of which is 1)
                # target base: snv -> the substituted base ; indel -> the base on right to the one corresponding to REF
            if len(self.query_pos_within_REF) == 0:
                self.LocaLt = None ; self.LocaRt = None
            else:
                self.LocaLt = self.query_pos_within_REF[0] + 1
                self.LocaRt = self.read.query_length - self.query_pos_within_REF[-1] 


    def set_five_prime_distance(self):
        if len(self.query_pos_within_REF):
            if self.read.is_reverse:
                self.Distance = self.read.query_length - self.query_pos_within_REF[-1] 
            else:
                self.Distance = self.query_pos_within_REF[0] + 1
        else:
            self.Distance = None


    def set_orientation(self):
        if self.read.mate_is_unmapped or not self.read.is_paired or self.read.is_secondary or self.read.is_supplementary:
            return
        if self.read.reference_name != self.read.next_reference_name:
            return

        if self.read.is_reverse: 
            self_orientation = 'R'
        else: 
            self_orientation = 'F'
        if self.read.mate_is_reverse:
            mate_orientation = 'R'
        else:
            mate_orientation = 'F'

        if self.read.is_read1:
            self_orientation += '1'
            mate_orientation += '2'
        elif self.read.is_read2:
            self_orientation += '2'
            mate_orientation += '1'
        else: # if read is not paired orientation is set to 0
            self_orientation += '0'
            mate_orientation += '0'

        if self.read.reference_start <= self.read.next_reference_start:
            self.orientation = self_orientation + mate_orientation
        else:
            self.orientation = mate_orientation + self_orientation
        return


    def set_lowqual_by_SV(self):
        if self.is_SVread:
            self.is_lowqual = True
            self.lowqual_cause.append('SVread')

                
    def set_lowqual_by_MM(self, MMfrac: float=0.3):
        if self.MM >= self.read.query_length * MMfrac:
            self.is_lowqual = True
            self.lowqual_cause.append('mismatch')
            #self.readclass = -1


    def set_lowqual_by_MQ(self, limMQ: int=20):
        if self.MQ != None:
            if self.read.mapping_quality <= limMQ:
                self.is_lowqual = True
                self.lowqual_cause.append('MQ')
            

    def set_lowqual_by_BQ(self, limBQ: int=20):
        if self.BQ != None:
            if self.BQ <= limBQ:
                self.is_lowqual = True
                self.lowqual_cause.append('BQ')
                #self.readclass = -1


class variant_edit():
    def __init__(self, CHROM: str, POS: int, REF: str, ALT: str, ref_ver: str='37', **args):
        self.CHROM = CHROM
        self.POS = POS
        self.REF = REF
        self.ALT = ALT
        self.ref_ver = ref_ver
        self.start = POS - 1
        self.end =  POS - 1 + max(len(REF), len(ALT))
        self.var_type = get_mttype(REF, ALT)
        if self.var_type == 'ins':
            self.indelSeq = self.ALT[1:]
        elif self.var_type == 'del':
            self.indelSeq = self.REF[1:]
        else:
            self.indelSeq = None
        for key, val in args.items():
            setattr(self, key, val)


    def check_satellite(self) -> bool:
        satellite_pr = satellite_pr_dict[self.ref_ver]
        variant_pr = pr.PyRanges(chromosomes = self.CHROM, starts = [self.start], ends = [self.end])
        self.is_satellite = not satellite_pr.intersect(variant_pr).empty
        return self.is_satellite


    def find_repeats(self, fasta: pysam.FastaFile, pre: int=10, post: int=40) -> list:
        seq_all = fasta.fetch(self.CHROM, max(0, self.POS-pre-1), self.POS+post)
        switch = max(0, self.POS-pre-1)
        L = len(seq_all)
        repeat_list = []
        repeat_coord_list = []
        for i in range(L):
            if True in [ (i in range(x[0], x[1])) for x in repeat_coord_list ]:
                continue
            
            for j in range(i, L):
                if True in [ (i in range(x[0], x[1])) for x in repeat_coord_list ]:
                    pass #continue
                subseq = seq_all[i:j+1]
                if subseq * 2 == seq_all[i:i+(j+1-i)*2]:
                    quot = 2 # quotient
                    for k in range(3, (L-i)//len(subseq)+1):
                        if subseq * k == seq_all[i:i+len(subseq)*k]:
                            quot = k
                        else:
                            break
                            
                    repeat_coord = (i, i + len(subseq)*quot)
                    repeat_coord_list.append(repeat_coord)
                    
                    repeat_coord_posedit = (repeat_coord[0]+switch, repeat_coord[1]+switch)
                    repeat_list.append((repeat_coord_posedit, subseq, quot))
        self.repeat_list = repeat_list
        return repeat_list # 0-based half-open


    def find_relevant_repeats(self, distance_threshold:int=1) -> list:
        '''
        1) Divide objects within repeat_list into upstream, including, downstream repeats.
        2) For upstream and downstrem, only the closest repeats are selected.
        3) Among all the repeats (maximum three) if distance between the variant and the repeat are smaller than the threshold, repeats are returned as a list. 
        '''
        repeat_including = []
        repeat_upstream = []
        repeat_downstream = []
        try:
            repeat_list = self.repeat_list
        except:
            raise ValueError("run find_repeats first!")
        
        for x in repeat_list:
            if x[0][1] > self.start and x[0][0] < self.end:
                repeat_including.append(x)
            elif x[0][1] <= self.start:
                repeat_upstream.append(x)
            elif x[0][0] >= self.end:
                repeat_downstream.append(x)
        
        repeat_upstream_chosen = repeat_upstream[-1] if len(repeat_upstream) > 0 else None
        repeat_downstream_chosen = repeat_downstream[0] if len(repeat_downstream) > 0 else None
        
        result = []
        
        if self.indelSeq != None: # with indel
            def indelSeq_in_repeat(repeat):
                return self.indelSeq in repeat[1]*repeat[2]
            
            for repeat in repeat_including:
                if indelSeq_in_repeat(repeat):
                    result.append(repeat)
            if repeat_upstream_chosen != None:
                if \
                self.start - repeat_upstream_chosen[0][1] + 1 <= distance_threshold and \
                indelSeq_in_repeat(repeat_upstream_chosen): # distance between upstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_upstream_chosen) 
            if repeat_downstream_chosen != None:
                if \
                repeat_downstream_chosen[0][0] - self.end + 1 <= distance_threshold and \
                indelSeq_in_repeat(repeat_downstream_chosen): # distance between downstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_downstream_chosen)
                    
        else: # not indel
            for repeat in repeat_including:
                result.append(repeat)
            if repeat_upstream_chosen != None:
                if self.start - repeat_upstream_chosen[0][1] + 1 <= distance_threshold: # distance between upstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_upstream_chosen)
            if repeat_downstream_chosen != None:
                if repeat_downstream_chosen[0][0] - self.end + 1 <= distance_threshold: # distance between downstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_downstream_chosen)
            
        self.relevant_repeats = result    
        return result


    def fetch_reads(self, bam: pysam.AlignmentFile, nreads: int=0):
        if nreads == 0:
            self.readlist = list(bam.fetch(contig=self.CHROM, start=self.start, end=self.end))
        else:
            self.readlist = list(bam.fetch(contig=self.CHROM, start=self.start, end=self.end))[0:nreads]
        self.nread = len(self.readlist)
        print_err(f'{self.CHROM}:{self.start}-{self.end} number of reads: {self.nread}')


    def reads_to_readplus(self):
        if self.nread == 0:
            self.rplist = None
            return
        self.rplist = []
        for read in self.readlist:
            rp = readplus(read)
            if rp.pairs != None:
                self.rplist.append(rp)
        return self.rplist


    def update_readplus_list(self, flankLen: int=1, 
        clip_threshold: int=0, 
        template_length_threshold: int=1000,
        MMfrac: float=0.1, 
        limMQ: int=20,
        limBQ: int=20):
        for rp in self.rplist:
            rp.set_irrelevant_reads(start=self.start, end=self.end, flankLen=flankLen)
            rp.set_irrelevant_reads_by_repeat(relevant_repeats=self.relevant_repeats)
            if not rp.is_irrelevant:
                rp.set_SVread(clip_threshold=clip_threshold, template_length_threshold=template_length_threshold)
                rp.set_lowqual_reads(limMQ)

            rp.set_pairs_list(start=self.start, end=self.end, flankLen=flankLen)

            if not rp.is_irrelevant:
                rp.set_alleleClass(REF=self.REF, ALT=self.ALT)
                rp.set_MM_BQ(var_type=self.var_type)
                rp.set_Loca(REF=self.REF, ALT=self.ALT, var_type=self.var_type)
                rp.set_lowqual_by_MM(MMfrac=MMfrac)
                rp.set_lowqual_by_BQ(limBQ=limBQ)


    def modify_readclass_by_mateoverlap(self):
        qname_list = [ rp.read.query_name for rp in self.rplist if not rp.read.is_supplementary ]
        mateovlp_qnames = [ qname for qname, cnt in Counter(qname_list).items() if cnt > 1 ]
        for qname in mateovlp_qnames:
            rp_list = [ rp for rp in self.rplist if rp.read.query_name == qname ]
            readclasses = [ rp.readclass for rp in rp_list ]
            if len(set(readclasses)) == 1:
                for rp in rp_list[1:]:
                    rp.is_irrelevant = True
                    rp.irrelevant_cause.append('mateoverlap_same_readclass')
                    #rp.readclass = -2
            else:
                for rp in rp_list:
                    rp.is_irrelevant = True
                    rp.irrelevant_cause.append('mateoverlap_different_readclass')
                    #rp.readclass = -2


    def modify_readclass_by_marginoverlap(self):
        if not self.nread in stackDict:
            return
        cutoff = stackDict[self.nread]
        if cutoff == '0':
            return
        
        margins_counter = Counter( [rp.margins for rp in self.rplist] )
        for rp in self.rplist:
            if margins_counter[rp.margins] >= cutoff:
                rp.is_irrelevant = True
                rp.irrelevant_cause.append('marginoverlap')
                    
                    
    def make_readclass_indices(self):
        self.readclass_indices_highqual = dict()
        self.readclass_indices_lowqual = dict()
        for i in [0, 1, 2]:
            self.readclass_indices_highqual[i] = list()
            self.readclass_indices_lowqual[i] = list()
        for rp in self.rplist:
            if not rp.is_irrelevant:
                if rp.is_lowqual:
                    self.readclass_indices_lowqual[rp.readclass].append( self.rplist.index(rp) )
                else:
                    self.readclass_indices_highqual[rp.readclass].append( self.rplist.index(rp) )


    def get_read_count(self):
        self.read_count = {}
        self.read_count['ref_read_all'] = 0
        self.read_count['var_read_all'] = 0
        self.read_count['other_read_all'] = 0
        self.read_count['ref_read_highqual'] = 0
        self.read_count['var_read_highqual'] = 0
        self.read_count['other_read_highqual'] = 0
        self.read_count['ref_read_lowqual'] = 0
        self.read_count['var_read_lowqual'] = 0
        self.read_count['other_read_lowqual'] = 0

        for rp in self.rplist:
            if rp.is_irrelevant:
                continue

            class_string = 'ref' if rp.readclass == 0 else ('var' if rp.readclass == 1 else 'other')
            qual_string = 'lowqual' if rp.is_lowqual else 'highqual'

            self.read_count[f'{class_string}_read_all'] += 1
            self.read_count[f'{class_string}_read_{qual_string}'] += 1
        
        Nread_all = \
                self.read_count['ref_read_highqual'] + self.read_count['var_read_highqual'] + self.read_count['other_read_highqual'] + \
                self.read_count['ref_read_lowqual'] + self.read_count['var_read_lowqual'] + self.read_count['other_read_lowqual']
        Nread_highqual = \
                self.read_count['ref_read_highqual'] + self.read_count['var_read_highqual'] + self.read_count['other_read_highqual']

        if Nread_all != 0:
            self.read_count['vaf_all'] = ( self.read_count['var_read_highqual'] + self.read_count['var_read_lowqual'] ) / Nread_all
        else: 
            self.read_count['vaf_all'] = np.nan
            
        if Nread_highqual != 0:
            self.read_count['vaf_highqual'] = self.read_count['var_read_highqual'] / Nread_highqual
        else: 
            self.read_count['vaf_highqual'] = np.nan

        return self.read_count


    def get_read_count_reverse(self):
        self.read_count_reverse = {}
        self.read_count_reverse['ref_read_all'] = 0
        self.read_count_reverse['ref_read_all_reverse'] = 0
        self.read_count_reverse['var_read_all'] = 0
        self.read_count_reverse['var_read_all_reverse'] = 0
        self.read_count_reverse['other_read_all'] = 0
        self.read_count_reverse['other_read_all_reverse'] = 0
        self.read_count_reverse['ref_read_highqual'] = 0
        self.read_count_reverse['ref_read_highqual_reverse'] = 0
        self.read_count_reverse['var_read_highqual'] = 0
        self.read_count_reverse['var_read_highqual_reverse'] = 0
        self.read_count_reverse['other_read_highqual'] = 0
        self.read_count_reverse['other_read_highqual_reverse'] = 0
        self.read_count_reverse['ref_read_lowqual'] = 0
        self.read_count_reverse['ref_read_lowqual_reverse'] = 0
        self.read_count_reverse['var_read_lowqual'] = 0
        self.read_count_reverse['var_read_lowqual_reverse'] = 0
        self.read_count_reverse['other_read_lowqual'] = 0
        self.read_count_reverse['other_read_lowqual_reverse'] = 0

        for rp in self.readsinfo.rplist:
            if rp.irrelevant:
                continue

            class_string = 'ref' if rp.readclass == 0 else ('var' if rp.readclass == 1 else 'other')
            qual_string = 'lowqual' if rp.lowqual else 'highqual'

            self.read_count_reverse[f'{class_string}_read_all'] += 1
            self.read_count_reverse[f'{class_string}_read_{qual_string}'] += 1

            if rp.read.is_reverse:
                self.read_count_reverse[f'{class_string}_read_all_reverse'] += 1
                self.read_count_reverse[f'{class_string}_read_{qual_string}_reverse'] += 1
        
        Nread_all = \
                self.read_count_reverse['ref_read_highqual'] + self.read_count_reverse['var_read_highqual'] + self.read_count_reverse['other_read_highqual'] + \
                self.read_count_reverse['ref_read_lowqual'] + self.read_count_reverse['var_read_lowqual'] + self.read_count_reverse['other_read_lowqual']
        Nread_highqual = \
                self.read_count_reverse['ref_read_highqual'] + self.read_count_reverse['var_read_highqual'] + self.read_count_reverse['other_read_highqual']

        if Nread_all != 0:
            self.read_count_reverse['vaf_all'] = ( self.read_count_reverse['var_read_highqual'] + self.read_count_reverse['var_read_lowqual'] ) / Nread_all
        else: 
            self.read_count_reverse['vaf_all'] = format_MISSING_float
            
        if Nread_highqual != 0:
            self.read_count_reverse['vaf_highqual'] = self.read_count_reverse['var_read_highqual'] / Nread_highqual
        else: 
            self.read_count_reverse['vaf_highqual'] = format_MISSING_float

        return self.read_count_reverse


    def get_MQ_BQ_MM(self):
        self.MQ_BQ_MM = {}
        self.MQ_BQ_MM['ref_all_meanMQ'] = 0
        self.MQ_BQ_MM['var_all_meanMQ'] = 0
        self.MQ_BQ_MM['ref_lowqual_meanMQ'] = 0
        self.MQ_BQ_MM['var_lowqual_meanMQ'] = 0
        self.MQ_BQ_MM['ref_all_medMQ'] = 0
        self.MQ_BQ_MM['var_all_medMQ'] = 0
        self.MQ_BQ_MM['ref_lowqual_medMQ'] = 0
        self.MQ_BQ_MM['var_lowqual_medMQ'] = 0
        self.MQ_BQ_MM['ref_all_meanBQ'] = 0
        self.MQ_BQ_MM['var_all_meanBQ'] = 0
        self.MQ_BQ_MM['ref_lowqual_meanBQ'] = 0
        self.MQ_BQ_MM['var_lowqual_meanBQ'] = 0
        self.MQ_BQ_MM['ref_all_medBQ'] = 0
        self.MQ_BQ_MM['var_all_medBQ'] = 0
        self.MQ_BQ_MM['ref_lowqual_medBQ'] = 0
        self.MQ_BQ_MM['var_lowqual_medBQ'] = 0
        self.MQ_BQ_MM['ref_all_meanMM'] = 0
        self.MQ_BQ_MM['var_all_meanMM'] = 0
        self.MQ_BQ_MM['ref_lowqual_meanMM'] = 0
        self.MQ_BQ_MM['var_lowqual_meanMM'] = 0
        self.MQ_BQ_MM['ref_all_medMM'] = 0
        self.MQ_BQ_MM['var_all_medMM'] = 0
        self.MQ_BQ_MM['ref_lowqual_medMM'] = 0
        self.MQ_BQ_MM['var_lowqual_medMM'] = 0

        def select_readclass_highqual(field, readclass):
            tmp = list()
            for idx in self.readclass_indices_highqual[readclass]:
                val = getattr(self.rplist[idx], field)
                if val != None:
                    tmp.append(val)
            return tmp

        def select_readclass_lowqual(field, readclass):
            tmp = list()
            for idx in self.readclass_indices_lowqual[readclass]:
                val = getattr(self.rplist[idx], field)
                if val != None:
                    tmp.append(val)
            return tmp

        self.MQ_BQ_MM['ref_all_meanMQ'], self.MQ_BQ_MM['ref_all_medMQ'] = calculate_mean_med(select_readclass_highqual('MQ', 0) + select_readclass_lowqual('MQ', 0))
        self.MQ_BQ_MM['ref_lowqual_meanMQ'], self.MQ_BQ_MM['ref_lowqual_medMQ'] = calculate_mean_med(select_readclass_lowqual('MQ', 0))
        self.MQ_BQ_MM['var_all_meanMQ'], self.MQ_BQ_MM['var_all_medMQ'] = calculate_mean_med(select_readclass_highqual('MQ', 1) + select_readclass_lowqual('MQ', 1))
        self.MQ_BQ_MM['var_lowqual_meanMQ'], self.MQ_BQ_MM['var_lowqual_medMQ'] = calculate_mean_med(select_readclass_lowqual('MQ', 1))
        self.MQ_BQ_MM['ref_all_meanBQ'], self.MQ_BQ_MM['ref_all_medBQ'] = calculate_mean_med(select_readclass_highqual('BQ', 0) + select_readclass_lowqual('BQ', 0))
        self.MQ_BQ_MM['ref_lowqual_meanBQ'], self.MQ_BQ_MM['ref_lowqual_medBQ'] = calculate_mean_med(select_readclass_lowqual('BQ', 0))
        self.MQ_BQ_MM['var_all_meanBQ'], self.MQ_BQ_MM['var_all_medBQ'] = calculate_mean_med(select_readclass_highqual('BQ', 1) + select_readclass_lowqual('BQ', 1))
        self.MQ_BQ_MM['var_lowqual_meanBQ'], self.MQ_BQ_MM['var_lowqual_medBQ'] = calculate_mean_med(select_readclass_lowqual('BQ', 1))
        self.MQ_BQ_MM['ref_all_meanMM'], self.MQ_BQ_MM['ref_all_medMM'] = calculate_mean_med(select_readclass_highqual('MM', 0) + select_readclass_lowqual('MM', 0))
        self.MQ_BQ_MM['ref_lowqual_meanMM'], self.MQ_BQ_MM['ref_lowqual_medMM'] = calculate_mean_med(select_readclass_lowqual('MM', 0))
        self.MQ_BQ_MM['var_all_meanMM'], self.MQ_BQ_MM['var_all_medMM'] = calculate_mean_med(select_readclass_highqual('MM', 1) + select_readclass_lowqual('MM', 1))
        self.MQ_BQ_MM['var_lowqual_meanMM'], self.MQ_BQ_MM['var_lowqual_medMM'] = calculate_mean_med(select_readclass_lowqual('MM', 1))

        return self.MQ_BQ_MM


    def get_Loca_Distance_Other(self):
        self.get_Loca_Distance_Other = {}
        self.get_Loca_Distance_Other['meanLocaLt_all'] = 0
        self.get_Loca_Distance_Other['medLocaLt_all'] = 0
        self.get_Loca_Distance_Other['meanLocaRt_all'] = 0
        self.get_Loca_Distance_Other['medLocaRt_all'] = 0
        self.get_Loca_Distance_Other['meanDistance_all'] = 0
        self.get_Loca_Distance_Other['medDistance_all'] = 0
        self.get_Loca_Distance_Other['ref_ClipCnt_all'] = 0
        self.get_Loca_Distance_Other['var_ClipCnt_all'] = 0
        self.get_Loca_Distance_Other['ref_InsCnt_all'] = 0
        self.get_Loca_Distance_Other['var_InsCnt_all'] = 0
        self.get_Loca_Distance_Other['ref_DelCnt_all'] = 0
        self.get_Loca_Distance_Other['var_DelCnt_all'] = 0

        def select_readclass_highqual(field, readclass):
            tmp = list()
            for idx in self.readclass_indices_highqual[readclass]:
                val = getattr(self.rplist[idx], field)
                if val != None:
                    tmp.append(val)
            return tmp

        def select_readclass_lowqual(field, readclass):
            tmp = list()
            for idx in self.readclass_indices_lowqual[readclass]:
                val = getattr(self.rplist[idx], field)
                if val != None:
                    tmp.append(val)
            return tmp

        self.get_Loca_Distance_Other['meanLocaLt_all'], self.get_Loca_Distance_Other['medLocaLt_all'] = calculate_mean_med(select_readclass_highqual('LocaLt', 1) + select_readclass_lowqual('LocaLt', 1))
        self.get_Loca_Distance_Other['meanLocaRt_all'], self.get_Loca_Distance_Other['medLocaRt_all'] = calculate_mean_med(select_readclass_highqual('LocaRt', 1) + select_readclass_lowqual('LocaRt', 1))
        self.get_Loca_Distance_Other['meanDistance_all'], self.get_Loca_Distance_Other['medDistance_all'] = calculate_mean_med(select_readclass_highqual('Distance', 1) + select_readclass_lowqual('Distance', 1))

        self.get_Loca_Distance_Other['ref_ClipCnt_all'] = len( [ y for y in select_readclass_highqual('lenClip', 0) + select_readclass_lowqual('lenClip', 0) if y > 0 ] )
        self.get_Loca_Distance_Other['var_ClipCnt_all'] = len( [ y for y in select_readclass_highqual('lenClip', 1) + select_readclass_lowqual('lenClip', 1) if y > 0 ] )
        self.get_Loca_Distance_Other['ref_InsCnt_all'] = len( [ y for y in select_readclass_highqual('lenIns', 0) + select_readclass_lowqual('lenIns', 0) if y > 0 ] )
        self.get_Loca_Distance_Other['var_InsCnt_all'] = len( [ y for y in select_readclass_highqual('lenIns', 1) + select_readclass_lowqual('lenIns', 1) if y > 0 ] )
        self.get_Loca_Distance_Other['ref_DelCnt_all'] = len( [ y for y in select_readclass_highqual('lenDel', 0) + select_readclass_lowqual('lenDel', 0) if y > 0 ] )
        self.get_Loca_Distance_Other['var_DelCnt_all'] = len( [ y for y in select_readclass_highqual('lenDel', 1) + select_readclass_lowqual('lenDel', 1) if y > 0 ] )

        return self.get_Loca_Distance_Other


    def get_info_values(self):
        pass


def position_to_vafall(CHROM, POS, REF, ALT, bam: pysam.AlignmentFile, ref_ver:str='37'):
    '''
    Chrom, Pos, Ref, Alt 만 알 때 빠르게 vaf_all 만 bam file 에서 알고 싶을 때 쓰면 좋을 것 같다.
    
    returns: dictionary of {'vaf_all': some value}

    '''
    variant = variant_edit(CHROM, POS, REF, ALT, ref_ver=ref_ver)
    variant.fetch_reads(bam)
    variant.reads_to_readplus()
    if not variant.rplist:
        result = {'vaf_all': 0}
        return 0
    variant.modify_readclass_by_mateoverlap()
    ref_var_other = {'ref':0, 'var':0, 'other':0}
    for rp in variant.rplist:
        rp.set_irrelevant_reads(start=variant.start, end=variant.end, flankLen=0)
        rp.set_pairs_list(start=variant.start, end=variant.end, flankLen=0)
        if rp.is_irrelevant:
            continue
        rp.set_alleleClass(REF=variant.REF, ALT=variant.ALT)
        class_string = 'ref' if rp.readclass == 0 else ('var' if rp.readclass == 1 else 'other')
        ref_var_other[class_string] += 1
    Nread_all = ref_var_other['ref'] + ref_var_other['var'] + ref_var_other['other']
    if Nread_all != 0:
        vaf_all = ref_var_other['var'] / Nread_all
    else:
        vaf_all = 0
    result = {'vaf_all': vaf_all}
    return vaf_all