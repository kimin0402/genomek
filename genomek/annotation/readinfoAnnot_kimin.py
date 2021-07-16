#!/home/users/pjh/conda_bin/python
##!/home/users/pjh/python/bin/python3

#for BWA
#updates
#180503: cigar loop update: error correction for deletion spanning the alt
#180518: cigar loop multiple error correction;add var_loca_reverse;change output format; get rlength from read
#180518: add informations: presence of clipping(soft of hard), INS, DEL in ref of var reads 
#180531: add informations: min, max location, refNM, varNM
#180717: indel is not available now ( write NA at all columns)
#180830: remove indel nt count from nm count
#180903: support gzip input
#181010: print NA when variant read count was zero.
#181022: indel can be annotated
#181101: unmapped read, duplicate read -> 'blank'
#181114: correction for multiple nucleotides in mutect2 such as CTT > CT, CT>AG
#181128: deletion length calculation error correction
#181129: deletion and insertion variant read detection error correction
#190221: count by readname
#190619: ins, del, clip count by readname.
#190729: detecting mnv and complex indel, median;min basequality
#190806: mnv error correction, variant condition for unknown inserted sequence
#190920: make mapq_info and mismatch_info to show min, median, and max value.
#191005: print sys.argv[0] at start, print 'done' at end
#191111: correct distance error when ins is located at the first base
#200409: rlength=read.infer_query_length() -> infer_read_length()
#200414: (pjh)
#    Revises alreadly present annotation columns. 
#    Finds columns named 'ref_medMQ var_medMQ' and changes it into 'ref_minMQ;med;max var_minMQ;med;max ref_med;meanBQ var_med;meanBQ' 
#    Finds columns named 'ref_medMismatch var_medMismatch' and changes it into 'ref_minMismatch;med;max var_minMismatch;med;max'
#    Saves output in a new file whose name includes 'readinfo2' instead of 'readinfo'
#200513: 1) Not patch, for new annotation ; 2) includes read-overlap and non-mapq0 annotations
#200518: count of each of ACGT in each variant position

#200521: modified version 1
    #annotate normal and tumor reads information at once

#200714: v2 initial
#200717: 1) context format changed into 'trinucleotide'>'alt' ; 2) other_base1, other_base2

#200818: v3 initial
    # added vaf only mode

#200903: init unified version
    # unified snv and indel 
    # revised program structure (using class, cyvcf2, etc.)
    
#200922: readsinfo.dict['MM'] value does not include indel length ; only includes base substitution count (excluding snv itself)

#200925 v2:
    # adds sample-specific information to FORMAT fields
    
#v2.2
    #201021
        # Uses pysam.AlignmentFile.fetch instead of pileup (fetch method by itself does not filter reads)
        # Made a new class which wraps pysam.AlignedSegment object and additional attributes ; this is used instead of readsinfo.dict
    #201023 : executes bcftools index on the resultant file
    #201024 : added nan_to_dot step at the end
    #201104 : when snv is within a repeat region, set all attributes as NA
    #201108
        # revised readclass system: non_reads are not a separate class but just another attribute. All reads except exclude-reads are given one of ref/var/other readclass.
        # reads not spanning over repeat region is considered irrelevant
    #201118
        # in get_pairs_list function, changed details as follows:
            #if x[1] == self.REF_start:
                #idx_REFstart = rp.pairs.index(x)
            #if x[1] == self.REF_end:
                #idx_nexttoREF = rp.pairs.index(x)
            #It does not matter if the base at REF_start or REF_end is deletion
        # get_repeat_coord function: self.repeatBases = set( list( self.fasta.fetch(self.chrom, self.REF_start, self.REF_end) ) )
    #210225
        # added 'N' base to compl_base dictionary
        # edited line 683 (in context annotation, fetched sequences are converted to uppercase)
        
    #210302
        # (v2.2.1) if check_flanking is False, allele class is set to 2 (previously treated as lowqual read)
        # v2.2.3 init
            # revised repeat search function
    #210303
        # vaf is now calculated from highqual & lowqual reads
        # vaf_highqual is calculated from only highqual reads
    #210304
        # default -O option value is change from 'b' to 'z'
    #210308
        # if ALT contains more than one alleles, split the line and annotate for each
    #210309
        # modifiy_readclass_by_mateoverlap function edited (supplementary reads are not considered)

    #210312
        # v2.2.4 init
        # added garbage collection
        # if fetch finds no read, ref_read, var_read, ... are set as 0 instead of MISSING value

    #210317
        # -f option accepts only '19' or '38'

    #210322
        # kimin: fasta.fetch error 해결
        # stackdict loading 한줄로 변경



import os
import sys
import re
import copy
import gc
import numpy as np
from collections import Counter
import pysam
import argparse
import cyvcf2
import pyranges as pr

info_MISSING = '.'
format_MISSING_int = -1 * 2**31
format_MISSING_float = np.nan
bcftools_path = '/home/users/pjh/tools/bcftools-1.10.2/bin/bcftools'
nantodot_path = '/home/users/pjh/scripts/vcf_editing/nan_to_dot.py'
#satellite_path_hg19 = '/home/users/pjh/References/ucsc_RepeatMasker_files/hg19/satellite_hg19_gencodeSubtract_renameChr.bed'
satellite_path_hg19 = '/home/users/pjh/References/ucsc_RepeatMasker_files/hg19/custom_files/satellite_hg19_rename_slop_subt_sort.bed.gz'
satellite_path_hg38 = '/home/users/pjh/References/ucsc_RepeatMasker_files/hg38/satellite_hg38_slop_subt_sort.bed.gz'
## Important ## Fasta file changed!!
# fasta_path_hg19 = '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
fasta_path_hg19 = '/home/users/kimin/projects/12_MT/MT_C5.fasta'
fasta_path_hg38 = '/home/users/pjh/References/reference_genome/GRCh38/GCA_for_alignment_pipelines/no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta'

# import satellite region file
pr_satellite = dict()
pr_satellite['hg19'] = pr.read_bed(satellite_path_hg19)
#pr_satellite['hg19'] = pr_satellite['hg19'].slack(500)
pr_satellite['hg38'] = pr.read_bed(satellite_path_hg38)
#pr_satellite['hg38'] = pr_satellite['hg38'].slack(500)


def printErr(*args, **kwargs):
    print(*args, file = sys.stderr, flush = True, **kwargs)
    
def summary(x):
    if x == None or x == []:
        return [format_MISSING_int, format_MISSING_int, format_MISSING_float, format_MISSING_int]
    else:
        return (min(x), int(np.median(x)), np.mean(x), max(x))
    
    
class readplus:
    def __init__(self, read):
        self.read = read
        
        try:
            self.pairs = self.read.get_aligned_pairs(with_seq = True)
        except ValueError: # when MD tag is absent, get_aligned_pairs method fails with ValueError
            self.pairs = None
            return
            
        self.readclass = None
        self.irrelevant = False
        self.irrelevant_cause = list()
        self.lowqual = False
        self.lowqual_cause = list()
        self.SVread = False
        self.SVread_cause = list()
        self.cigarstats = self.read.get_cigar_stats()[0]
        self.lenIns = self.cigarstats[1]
        self.lenDel = self.cigarstats[2]
        self.lenClip = self.cigarstats[4]
        self.MQ = self.read.mapping_quality
        self.BQ = None
        self.margins = ( self.read.reference_start, self.read.reference_end )
    
    
class makeReadsinfo:
    
    ############
    ### INIT ###
    ############
    
    def __init__(self, **args):
    # all lim_* cutoffs are inclusive (lowqual if identical to cutoff)
    # flankLen : length (one side) of bases flanking inserted/deleted site which must be cigar M to be selected as a valid read
    # MMfrac : if a read has 'MM' value higher than readlen * MMfrac, its readclass will be -1
    
        for key, val in args.items():
            setattr(self, key, val)
            
            
        # check satellite area
        pr_variant = pr.PyRanges(chromosomes = self.chrom, starts = [self.REF_start], ends = [self.REF_end])
        if self.fopt == '19':
            self.satellite = not pr_satellite['hg19'].intersect(pr_variant).empty
        elif self.fopt == '38':
            self.satellite = not pr_satellite['hg38'].intersect(pr_variant).empty
        else:
            self.satellite = False
            
            
        # set basic attributes
        if self.mttype == 'ins':
            self.indelSeq = self.alt[1:]
        elif self.mttype == 'del':
            self.indelSeq = self.ref[1:]
        else:
            self.indelSeq = None
            
        self.repeat_list = self.find_repeats(self.fasta, self.chrom, self.pos)
        self.relevant_repeats = self.get_relevant_repeats(self.repeat_list)
        self.snv_withinRepeat = False
        
        
        # stop if satellite
        if self.satellite == True:
            self.rplist = None
            printErr(f'{self.chrom} {self.pos} satellite')
            return
            
            
        # fetch reads
        fetch = self.bam.fetch(region = f'{self.chrom}:{self.pos}-{self.pos}')
        self.readlist = list(fetch)
        self.nread = len(self.readlist)
        
        printErr(f'{self.chrom} {self.pos} {self.nread}')
    
        # create rplist
        if self.nread == 0:
            self.rplist = None
        else:
            self.rplist = list()
            for read in self.readlist:
                rp = readplus(read)
                if rp.pairs != None:
                    self.add_information(rp)
                    self.rplist.append(rp)
                
                    
            self.modify_readclass_by_mateoverlap()
            # self.modify_readclass_by_marginoverlap() ## kimin: botseq 위해서 주석처리 함
                
            self.make_readclass_indices()
            
                
                
    def find_repeats(self, fasta, chrom, pos, pre=10, post=40): # pos: 1-based
        seq_all = fasta.fetch(chrom, max(0, pos-pre-1), pos+post) # fasta.fetch raises error when fetch position is negative integer
        switch = max(0, pos-pre-1)
        L = len(seq_all)
        repeat_list = list()
        repeat_coord_list = list()
        
        for i in range(L):
            if True in [ (i in range(x[0], x[1])) for x in repeat_coord_list ]:
                continue
            
            for j in range(i, L):
                if True in [ (i in range(x[0], x[1])) for x in repeat_coord_list ]:
                    break
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
                    
        return(repeat_list) # 0-based half-open
                
        
    def get_relevant_repeats(self, repeat_list):
        repeat_including = list()
        repeat_upstream = list()
        repeat_downstream = list()
        
        for x in repeat_list:
            if x[0][1] > self.REF_start and x[0][0] < self.REF_end:
                repeat_including.append(x)
            elif x[0][1] <= self.REF_start:
                repeat_upstream.append(x)
            elif x[0][0] >= self.REF_end:
                repeat_downstream.append(x)
        
        repeat_upstream_chosen = repeat_upstream[-1] if len(repeat_upstream) > 0 else None
        repeat_downstream_chosen = repeat_downstream[0] if len(repeat_downstream) > 0 else None
        
        
        result = list()
        
        if self.indelSeq != None: # with indel
            def indelSeq_in_repeat(repeat):
                if self.indelSeq in repeat[1]*repeat[2]:
                    return True
                else:
                    return False
            
            for repeat in repeat_including:
                if indelSeq_in_repeat(repeat):
                    result.append(repeat)
            if repeat_upstream_chosen != None:
                if \
                self.REF_start - repeat_upstream_chosen[0][1] + 1 <= 1 and \
                indelSeq_in_repeat(repeat_upstream_chosen): # distance between upstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_upstream_chosen) 
            if repeat_downstream_chosen != None:
                if \
                repeat_downstream_chosen[0][0] - self.REF_end + 1 <= 1 and \
                indelSeq_in_repeat(repeat_downstream_chosen): # distance between downstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_downstream_chosen)
                    
            return(result)
        
        else: # not indel
            for repeat in repeat_including:
                result.append(repeat)
            if repeat_upstream_chosen != None:
                if self.REF_start - repeat_upstream_chosen[0][1] + 1 <= 1: # distance between upstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_upstream_chosen)
            if repeat_downstream_chosen != None:
                if repeat_downstream_chosen[0][0] - self.REF_end + 1 <= 1: # distance between downstream repeat and snv site ; equals to 1 with adjacency
                    result.append(repeat_downstream_chosen)
                
            return(result)
        
        
    def modify_readclass_by_mateoverlap(self):
        qname_list = [ rp.read.query_name for rp in self.rplist if not rp.read.is_supplementary ]
        mateovlp_qnames = [ qname for qname, cnt in Counter(qname_list).items() if cnt > 1 ]
        for qname in mateovlp_qnames:
            rp_list = [ rp for rp in self.rplist if rp.read.query_name == qname ]
            readclasses = [ rp.readclass for rp in rp_list ]
            if len(set(readclasses)) == 1:
                for rp in rp_list[1:]:
                    rp.irrelevant = True
                    rp.irrelevant_cause.append('mateoverlap_same_readclass')
                    #rp.readclass = -2
            else:
                for rp in rp_list:
                    rp.irrelevant = True
                    rp.irrelevant_cause.append('mateoverlap_different_readclass')
                    #rp.readclass = -2
                    
                    
    def modify_readclass_by_marginoverlap(self):
        if not self.nread in self.stackDict:
            return
        cutoff = self.stackDict[self.nread]
        if cutoff == '0':
            return
        
        margins_counter = Counter( [rp.margins for rp in self.rplist] )
        for rp in self.rplist:
            if margins_counter[rp.margins] >= cutoff:
                rp.irrelevant = True
                rp.irrelevant_cause.append('marginoverlap')
                    
                    
    def make_readclass_indices(self):
        self.readclass_indices = dict()
        self.readclass_indices_lowqual = dict()
        for i in [0, 1, 2]:
            self.readclass_indices[i] = list()
            self.readclass_indices_lowqual[i] = list()
        for rp in self.rplist:
            if rp.irrelevant == False:
                if rp.lowqual:
                    self.readclass_indices_lowqual[rp.readclass].append( self.rplist.index(rp) )
                else:
                    self.readclass_indices[rp.readclass].append( self.rplist.index(rp) )
            
            
    def add_information(self, rp):
        
        def set_irrelevant_reads():
            if not (self.REF_start - rp.read.reference_start >= self.flankLen and rp.read.reference_end - self.REF_end >= self.flankLen ):
                rp.irrelevant = True
                rp.irrelevant_cause.append('insufficient_flanking_sequences')
            elif rp.read.has_tag('MD') == False:
                rp.irrelevant = True
                rp.irrelevant_cause.append('absent_MD_tag')
            elif rp.read.is_unmapped == True:
                rp.irrelevant = True
                rp.irrelevant_cause.append('flag_unmapped')
            elif rp.read.is_duplicate == True:
                rp.irrelevant = True
                rp.irrelevant_cause.append('flag_duplicate')
            elif rp.read.is_supplementary == True:
                rp.irrelevant = True
                rp.irrelevant_cause.append('flag_supplementary')
            elif rp.read.is_secondary == True:
                rp.irrelevant = True
                rp.irrelevant_cause.append('flag_secondary')
            elif rp.read.is_qcfail == True:
                rp.irrelevant = True
                rp.irrelevant_cause.append('flag_qcfail')
            else:
                for repeat in self.relevant_repeats:
                    if not ( rp.read.reference_start < repeat[0][0] and rp.read.reference_end > repeat[0][1] ):
                        rp.irrelevant = True
                        rp.irrelevant_cause.append('insufficient_repeat_span')
                        break
                
                
        def set_SVread():
            if rp.lenClip > 0:
                rp.SVread = True
                rp.SVread_cause.append('clipping')
            if rp.read.is_proper_pair == False:
                rp.SVread = True
                rp.SVread_cause.append('discordant pair')
            if rp.read.reference_name != rp.read.next_reference_name:
                rp.SVread = True
                rp.SVread_cause.append('Mate on difference chromosome')
            if abs(rp.read.template_length) > 1000:
                rp.SVread = True
                rp.SVread_cause.append('template length > 1000')
            if rp.read.has_tag('SA'):
                rp.SVread = True
                rp.SVread_cause.append('SA tag')
                
        def set_lowqual_reads():
            # if rp.SVread == True:  03.23 kimin: template 사이즈가 큰 read 의 snp 도 잡기 위해서
            #     rp.lowqual = True
            #     rp.lowqual_cause.append('SVread')
            if rp.read.mapping_quality <= self.limMQ:
                rp.lowqual = True
                rp.lowqual_cause.append('MQ')
                        
        def get_pairs_list():
            if rp.irrelevant == True:
                rp.pairs_before_REF = None
                rp.pairs_after_REF  = None
                rp.pairs_within_REF = None
                rp.query_pos_within_REF = None
                rp.seq_within_REF = None
                
            else:
                for x in rp.pairs:
                    if x[1] is None:
                        continue
                    if x[1] == self.REF_start:
                        idx_REFstart = rp.pairs.index(x)
                    if x[1]+1 == self.REF_end:
                        idx_nexttoREF = rp.pairs.index(x)+1
                        break
                
                rp.pairs_before_REF = rp.pairs[idx_REFstart - self.flankLen : idx_REFstart]
                rp.pairs_after_REF  = rp.pairs[idx_nexttoREF : idx_nexttoREF + self.flankLen]
                rp.pairs_within_REF = rp.pairs[idx_REFstart : idx_nexttoREF]
                
                rp.query_pos_within_REF = [ x[0] for x in rp.pairs_within_REF if x[0] != None ] 
                    # 0-based
                    # this can be an empty list
                    
                rp.seq_within_REF = ''.join([ rp.read.query_sequence[idx] for idx in rp.query_pos_within_REF ])
                    # this can be an empty string
            
                        
        def set_alleleClass():
            # determine a read supports which allele ; executed only when readclass is neither -2 nor -1
            # 0: ref supporting ; 1: alt supporting ; 2: other supporting
            def check_flanking(pairs_before_REF, pairs_after_REF):
                # check if all flanking sequences are cigar M with identical bases
                for x in pairs_before_REF + pairs_after_REF:
                    if x[0] == None or x[1] == None or x[2].islower():
                        return False
                return True
            
            if check_flanking(rp.pairs_before_REF, rp.pairs_after_REF) == False:
                #rp.lowqual = True
                rp.readclass = 2
            else:
                if rp.seq_within_REF == self.ref:
                    rp.readclass = 0
                elif rp.seq_within_REF == self.alt:
                    rp.readclass = 1
                else:
                    rp.readclass = 2
                
        def addInfo_basic():
            # NM tag value includes indels as well as base mismatch ; indel length is subtracted
            rp.MM = rp.read.get_tag('NM') - rp.cigarstats[1] - rp.cigarstats[2] if rp.read.has_tag('NM') else None
            # values attributed to the mutation itself are subtracted
            if rp.readclass == 1 and self.mttype == 'snv':
                rp.MM -= 1
            
            if len(rp.query_pos_within_REF) != 0:
                rp.BQ = np.mean( [ rp.read.query_qualities[idx] for idx in rp.query_pos_within_REF ] )
            
        def addInfo_detail1():
            if rp.readclass == 1:
                if self.mttype == 'ins':
                    rp.lenIns -= ( len(self.alt) - len(self.ref) )
                elif self.mttype == 'del':
                    rp.lenDel -= ( len(self.ref) - len(self.alt) )
                    
                # calculating LocaLt/LocaRt (1-based, including softclip)
                # LocaLt: index of the target base counted from the left (index of the leftmost base is 1)
                    # target base: snv -> the substituted base ; indel -> the base corresponding to REF
                # LocaRt: counting from the rightmost base (index of which is 1)
                    # target base: snv -> the substituted base ; indel -> the base on right to the one corresponding to REF
                if len(rp.query_pos_within_REF) == 0:
                    rp.LocaLt = None ; rp.LocaRt = None
                else:
                    rp.LocaLt = rp.query_pos_within_REF[0] + 1
                    rp.LocaRt = rp.read.query_length - rp.query_pos_within_REF[-1] 
                    
        def set_lowqual_by_MM():
            if rp.MM >= rp.read.query_length * self.MMfrac:
                rp.lowqual = True
                rp.lowqual_cause.append('mismatch')
                #rp.readclass = -1
                
        def set_lowqual_by_BQ():
            if rp.BQ != None:
                if rp.BQ <= self.limBQ:
                    rp.lowqual = True
                    rp.lowqual_cause.append('BQ')
                    #rp.readclass = -1
                    
        
        def main():
            set_irrelevant_reads() 
            if rp.irrelevant == False:
                set_SVread()
                set_lowqual_reads()
                
            get_pairs_list()
            #if self.mttype == 'snv':
            #    check_indelWithinRepeat()
                
            if rp.irrelevant == False:
                set_alleleClass()
                addInfo_basic()
                addInfo_detail1()
                set_lowqual_by_MM()
                set_lowqual_by_BQ()
                
        main()
                
                
                
            
            
            
class makeINFO:
    def __init__(self, **args):
        # default parameters
        self.makeMeta = False
        self.makeNull = False
        # parse arguments
        for key, val in args.items():
            setattr(self, key, val)
        
        if self.makeMeta == True:
            self.info_meta = list()
            self.format_meta = list()
        elif self.makeMeta == False:
            self.info_values = dict()
            self.format_values = dict()
        
        # add information to dict
        if self.level == 2:
            self.addInfo_readcount()
            self.addInfo_detail1()
            self.addInfo_context()
            self.addInfo_repeatRegion()
            self.addInfo_satellite()
        elif self.level == 1:
            self.addInfo_readcount()
            self.addInfo_satellite()
            
    def set_missing_all(self, info_meta_list, format_meta_list):
        for field in [ x['ID'] for x in info_meta_list ]:
            self.info_values[field] = info_MISSING
        for dic in format_meta_list:
            if dic['Type'] == 'Integer':
                self.format_values[dic['ID']] = format_MISSING_int
            elif dic['Type'] == 'Float':
                self.format_values[dic['ID']] = format_MISSING_float
                
                
    def addInfo_readcount(self):
        info_meta_list = []
        format_meta_list = [
            {'ID':'ref_read', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting ref allele'},
            {'ID':'var_read', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting alt allele'},
            {'ID':'other_read', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting another allele different from ref and alt'},
            #{'ID':'non_read', 'Type':'Integer', 'Number':1, 'Description':'Number of reads not supporting any haplotype'},
            {'ID':'ref_read_lowqual', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting ref allele'},
            {'ID':'var_read_lowqual', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting alt allele among lowqual reads'},
            {'ID':'other_read_lowqual', 'Type':'Integer', 'Number':1, 'Description':'Number of reads supporting another allele different from ref and alt among lowqual reads'},
            {'ID':'vaf', 'Type':'Float', 'Number':1, 'Description':'(var_read + var_read_lowqual) / (ref_read+var_read+other_read+ref_read_lowqual+var_read_lowqual+other_read_lowqual) '},
            {'ID':'vaf_highqual', 'Type':'Float', 'Number':1, 'Description':'var_read / (ref_read + var_read + other_read) '},
        ]
        
        if self.makeMeta:
            self.info_meta.extend(info_meta_list)
            self.format_meta.extend(format_meta_list)
            
        else:
            if self.makeNull:
                self.set_missing_all(info_meta_list, format_meta_list)

            elif self.readsinfo.rplist == None:
                self.format_values['ref_read'] = 0
                self.format_values['var_read'] = 0
                self.format_values['other_read'] = 0
                self.format_values['ref_read_lowqual'] = 0
                self.format_values['var_read_lowqual'] = 0
                self.format_values['other_read_lowqual'] = 0
                self.format_values['vaf'] = format_MISSING_float
                self.format_values['vaf_highqual'] = format_MISSING_float

            else:
                self.format_values['ref_read'] = len(self.readsinfo.readclass_indices[0])
                self.format_values['var_read'] = len(self.readsinfo.readclass_indices[1])
                self.format_values['other_read'] = len(self.readsinfo.readclass_indices[2])
                self.format_values['ref_read_lowqual'] = len(self.readsinfo.readclass_indices_lowqual[0])
                self.format_values['var_read_lowqual'] = len(self.readsinfo.readclass_indices_lowqual[1])
                self.format_values['other_read_lowqual'] = len(self.readsinfo.readclass_indices_lowqual[2])
                
                Nread = self.format_values['ref_read'] + self.format_values['var_read'] + self.format_values['other_read'] + self.format_values['ref_read_lowqual'] + self.format_values['var_read_lowqual'] + self.format_values['other_read_lowqual']
                Nread_highqual = self.format_values['ref_read'] + self.format_values['var_read'] + self.format_values['other_read']
                if Nread != 0:
                    self.format_values['vaf'] = ( self.format_values['var_read'] + self.format_values['var_read_lowqual'] ) / Nread
                else: 
                    self.format_values['vaf'] = format_MISSING_float
                    
                if Nread_highqual != 0:
                    self.format_values['vaf_highqual'] = self.format_values['var_read'] / Nread_highqual
                else: 
                    self.format_values['vaf_highqual'] = format_MISSING_float
                    
                    
    def addInfo_detail1(self):
        info_meta_list = []
        format_meta_list = [
            {'ID':'ref_minMQ', 'Type':'Integer', 'Number':1, 'Description':'minimum MQ among all ref-supporting reads'},
            {'ID':'ref_medMQ', 'Type':'Integer', 'Number':1, 'Description':'median MQ among all ref-supporting reads'},
            {'ID':'ref_meanMQ', 'Type':'Float', 'Number':1, 'Description':'mean MQ among all ref-supporting reads'},
            {'ID':'ref_maxMQ', 'Type':'Integer', 'Number':1, 'Description':'maximum MQ among all ref-supporting reads'},
            {'ID':'var_minMQ', 'Type':'Integer', 'Number':1, 'Description':'minimum MQ among all alt-supporting reads'},
            {'ID':'var_medMQ', 'Type':'Integer', 'Number':1, 'Description':'median MQ among all alt-supporting reads'},
            {'ID':'var_meanMQ', 'Type':'Float', 'Number':1, 'Description':'mean MQ among all alt-supporting reads'},
            {'ID':'var_maxMQ', 'Type':'Integer', 'Number':1, 'Description':'maximum MQ among all alt-supporting reads'},
            #{'ID':'non_minMQ', 'Type':'Integer', 'Number':1, 'Description':'minimum MQ among all non-supporting reads'},
            #{'ID':'non_medMQ', 'Type':'Integer', 'Number':1, 'Description':'median MQ among all non-supporting reads'},
            #{'ID':'non_meanMQ', 'Type':'Float', 'Number':1, 'Description':'mean MQ among all non-supporting reads'},
            #{'ID':'non_maxMQ', 'Type':'Integer', 'Number':1, 'Description':'maximum MQ among all non-supporting reads'},
            
            {'ID':'ref_minBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) minimum BQ among all ref-supporting reads'},
            {'ID':'ref_medBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) median BQ among all ref-supporting reads'},
            {'ID':'ref_meanBQ', 'Type':'Float', 'Number':1, 'Description':'(snv only) mean BQ among all ref-supporting reads'},
            {'ID':'ref_maxBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) maximum BQ among all ref-supporting reads'},
            {'ID':'var_minBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) minimum BQ among all alt-supporting reads'},
            {'ID':'var_medBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) median BQ among all alt-supporting reads'},
            {'ID':'var_meanBQ', 'Type':'Float', 'Number':1, 'Description':'(snv only) mean BQ among all alt-supporting reads'},
            {'ID':'var_maxBQ', 'Type':'Integer', 'Number':1, 'Description':'(snv only) maximum BQ among all alt-supporting reads'},
        
            {'ID':'minLocaLt', 'Type':'Integer', 'Number':1, 'Description':'minimum of position of the target base counted from left end of the read among alt-supporting reads'},
            {'ID':'medLocaLt', 'Type':'Integer', 'Number':1, 'Description':'median of position of the target base counted from left end of the read among alt-supporting reads'},
            {'ID':'meanLocaLt', 'Type':'Float', 'Number':1, 'Description':'mean of position of the target base counted from left end of the read among alt-supporting reads'},
            {'ID':'maxLocaLt', 'Type':'Integer', 'Number':1, 'Description':'maximum of position of the target base counted from left end of the read among alt-supporting reads'},
            {'ID':'minLocaRt', 'Type':'Integer', 'Number':1, 'Description':'minimum of position of the target base counted from right end of the read among alt-supporting reads'},
            {'ID':'medLocaRt', 'Type':'Integer', 'Number':1, 'Description':'median of position of the target base counted from right end of the read among alt-supporting reads'},
            {'ID':'meanLocaRt', 'Type':'Float', 'Number':1, 'Description':'mean of position of the target base counted from right end of the read among alt-supporting reads'},
            {'ID':'maxLocaRt', 'Type':'Integer', 'Number':1, 'Description':'maximum of position of the target base counted from right end of the read among alt-supporting reads'},
            
            {'ID':'ref_ClipCnt', 'Type':'Integer', 'Number':1, 'Description':'number of ref-supporting reads harboring soft clipping'},
            {'ID':'ref_ClipFrac', 'Type':'Float', 'Number':1, 'Description':'ref_ClipCnt / ref_read'},
            {'ID':'var_ClipCnt', 'Type':'Integer', 'Number':1, 'Description':'number of alt-supporting reads harboring soft clipping'},
            {'ID':'var_ClipFrac', 'Type':'Float', 'Number':1, 'Description':'var_ClipCnt / var_read'},
            #{'ID':'non_ClipCnt', 'Type':'Integer', 'Number':1, 'Description':'number of non-supporting reads harboring soft clipping'},
            #{'ID':'non_ClipFrac', 'Type':'Float', 'Number':1, 'Description':'non_ClipCnt / non_read'},
            {'ID':'ref_InsCnt', 'Type':'Integer', 'Number':1, 'Description':'number of ref-supporting reads harboring insertion'},
            {'ID':'ref_InsFrac', 'Type':'Float', 'Number':1, 'Description':'ref_InsCnt / ref_read'},
            {'ID':'var_InsCnt', 'Type':'Integer', 'Number':1, 'Description':'number of alt-supporting reads harboring insertion'},
            {'ID':'var_InsFrac', 'Type':'Float', 'Number':1, 'Description':'var_InsCnt / var_read'},
            #{'ID':'non_InsCnt', 'Type':'Integer', 'Number':1, 'Description':'number of non-supporting reads harboring insertion'},
            #{'ID':'non_InsFrac', 'Type':'Float', 'Number':1, 'Description':'non_InsCnt / non_read'},
            {'ID':'ref_DelCnt', 'Type':'Integer', 'Number':1, 'Description':'number of ref-supporting reads harboring deletion'},
            {'ID':'ref_DelFrac', 'Type':'Float', 'Number':1, 'Description':'ref_DelCnt / ref_read'},
            {'ID':'var_DelCnt', 'Type':'Integer', 'Number':1, 'Description':'number of alt-supporting reads harboring deletion'},
            {'ID':'var_DelFrac', 'Type':'Float', 'Number':1, 'Description':'var_DelCnt / var_read'},
            #{'ID':'non_DelCnt', 'Type':'Integer', 'Number':1, 'Description':'number of non-supporting reads harboring deletion'},
            #{'ID':'non_DelFrac', 'Type':'Float', 'Number':1, 'Description':'non_DelCnt / non_read'},
        
            {'ID':'ref_minMM', 'Type':'Integer', 'Number':1, 'Description':'minimum mismatch (modified NM tag value) among all ref-supporting reads'},
            {'ID':'ref_medMM', 'Type':'Integer', 'Number':1, 'Description':'median mismatch (modified NM tag value) among all ref-supporting reads'},
            {'ID':'ref_meanMM', 'Type':'Float', 'Number':1, 'Description':'mean mismatch (modified NM tag value) among all ref-supporting reads'},
            {'ID':'ref_maxMM', 'Type':'Integer', 'Number':1, 'Description':'maximum mismatch (modified NM tag value) among all ref-supporting reads'},
            {'ID':'var_minMM', 'Type':'Integer', 'Number':1, 'Description':'minimum mismatch (modified NM tag value) among all alt-supporting reads'},
            {'ID':'var_medMM', 'Type':'Integer', 'Number':1, 'Description':'median mismatch (modified NM tag value) among all alt-supporting reads'},
            {'ID':'var_meanMM', 'Type':'Float', 'Number':1, 'Description':'mean mismatch (modified NM tag value) among all alt-supporting reads'},
            {'ID':'var_maxMM', 'Type':'Integer', 'Number':1, 'Description':'maximum mismatch (modified NM tag value) among all alt-supporting reads'},
            #{'ID':'non_minMM', 'Type':'Integer', 'Number':1, 'Description':'minimum mismatch (modified NM tag value) among all non-supporting reads'},
            #{'ID':'non_medMM', 'Type':'Integer', 'Number':1, 'Description':'median mismatch (modified NM tag value) among all non-supporting reads'},
            #{'ID':'non_meanMM', 'Type':'Float', 'Number':1, 'Description':'mean mismatch (modified NM tag value) among all non-supporting reads'},
            #{'ID':'non_maxMM', 'Type':'Integer', 'Number':1, 'Description':'maximum mismatch (modified NM tag value) among all non-supporting reads'},
        
            {'ID':'N_Loca', 'Type':'Integer', 'Number':1, 'Description':'number of different target positions of var-supporting reads'},
        ]
        
        if self.makeMeta:
            self.info_meta.extend(info_meta_list)
            self.format_meta.extend(format_meta_list)
        
        else:
            if self.makeNull:
                self.set_missing_all(info_meta_list, format_meta_list)
            elif self.readsinfo.rplist == None:
                self.set_missing_all(info_meta_list, format_meta_list)
            else:
                def select_readclass(field, readclass):
                    tmp = list()
                    for idx in self.readsinfo.readclass_indices[readclass]:
                        val = getattr(self.readsinfo.rplist[idx], field)
                        if val != None:
                            tmp.append(val)
                    return tmp
                    
                self.format_values['ref_minMQ'], self.format_values['ref_medMQ'], self.format_values['ref_meanMQ'], self.format_values['ref_maxMQ'] = summary( select_readclass('MQ', 0) )
                self.format_values['var_minMQ'], self.format_values['var_medMQ'], self.format_values['var_meanMQ'], self.format_values['var_maxMQ'] = summary( select_readclass('MQ', 1) )
                #self.format_values['non_minMQ'], self.format_values['non_medMQ'], self.format_values['non_meanMQ'], self.format_values['non_maxMQ'] = summary( select_readclass('MQ', -1) )
                
                if self.mttype == 'snv':
                    self.format_values['ref_minBQ'], self.format_values['ref_medBQ'], self.format_values['ref_meanBQ'], self.format_values['ref_maxBQ'] = summary( select_readclass('BQ', 0) )
                    self.format_values['var_minBQ'], self.format_values['var_medBQ'], self.format_values['var_meanBQ'], self.format_values['var_maxBQ'] = summary( select_readclass('BQ', 1) )
                
                self.format_values['minLocaLt'], self.format_values['medLocaLt'], self.format_values['meanLocaLt'], self.format_values['maxLocaLt'] = summary( select_readclass('LocaLt', 1) )
                self.format_values['minLocaRt'], self.format_values['medLocaRt'], self.format_values['meanLocaRt'], self.format_values['maxLocaRt'] = summary( select_readclass('LocaRt', 1) )
                
                self.format_values['ref_minMM'], self.format_values['ref_medMM'], self.format_values['ref_meanMM'], self.format_values['ref_maxMM'] = summary( select_readclass('MM', 0) )
                self.format_values['var_minMM'], self.format_values['var_medMM'], self.format_values['var_meanMM'], self.format_values['var_maxMM'] = summary( select_readclass('MM', 1) )
                #self.format_values['non_minMM'], self.format_values['non_medMM'], self.format_values['non_meanMM'], self.format_values['non_maxMM'] = summary( select_readclass('MM', -1) )
                
                self.format_values['ref_ClipCnt'] = len( [ y for y in select_readclass('lenClip', 0) if y > 0 ] )
                self.format_values['ref_ClipFrac'] = self.format_values['ref_ClipCnt'] / self.format_values['ref_read'] if self.format_values['ref_read'] != 0 else format_MISSING_float
                self.format_values['var_ClipCnt'] = len( [ y for y in select_readclass('lenClip', 1) if y > 0 ] )
                self.format_values['var_ClipFrac'] = self.format_values['var_ClipCnt'] / self.format_values['var_read'] if self.format_values['var_read'] != 0 else format_MISSING_float
                #self.format_values['non_ClipCnt'] = len( [ y for y in select_readclass('lenClip', -1) if y > 0 ] )
                #self.format_values['non_ClipFrac'] = self.format_values['non_ClipCnt'] / self.format_values['non_read'] if self.format_values['non_read'] != 0 else format_MISSING_float
                
                self.format_values['ref_InsCnt'] = len( [ y for y in select_readclass('lenIns', 0) if y > 0 ] )
                self.format_values['ref_InsFrac'] = self.format_values['ref_InsCnt'] / self.format_values['ref_read'] if self.format_values['ref_read'] != 0 else format_MISSING_float
                self.format_values['var_InsCnt'] = len( [ y for y in select_readclass('lenIns', 1) if y > 0 ] )
                self.format_values['var_InsFrac'] = self.format_values['var_InsCnt'] / self.format_values['var_read'] if self.format_values['var_read'] != 0 else format_MISSING_float
                #self.format_values['non_InsCnt'] = len( [ y for y in select_readclass('lenIns', -1) if y > 0 ] )
                #self.format_values['non_InsFrac'] = self.format_values['non_InsCnt'] / self.format_values['non_read'] if self.format_values['non_read'] != 0 else format_MISSING_float
                
                self.format_values['ref_DelCnt'] = len( [ y for y in select_readclass('lenDel', 0) if y > 0 ] )
                self.format_values['ref_DelFrac'] = self.format_values['ref_DelCnt'] / self.format_values['ref_read'] if self.format_values['ref_read'] != 0 else format_MISSING_float
                self.format_values['var_DelCnt'] = len( [ y for y in select_readclass('lenDel', 1) if y > 0 ] )
                self.format_values['var_DelFrac'] = self.format_values['var_DelCnt'] / self.format_values['var_read'] if self.format_values['var_read'] != 0 else format_MISSING_float
                #self.format_values['non_DelCnt'] = len( [ y for y in select_readclass('lenDel', -1) if y > 0 ] )
                #self.format_values['non_DelFrac'] = self.format_values['non_DelCnt'] / self.format_values['non_read'] if self.format_values['non_read'] != 0 else format_MISSING_float
                
                self.format_values['N_Loca'] = len( set( select_readclass('LocaLt', 1) ) )
        
            
    def addInfo_context(self):
        info_meta_list = [
            {'ID':'context', 'Type':'String', 'Number':1, 'Description':'(snv only) trinucleotide context around the variant position. Format: NNN>N'}
        ]
        format_meta_list = []
        
        if self.makeMeta:
            self.info_meta.extend(info_meta_list)
            self.format_meta.extend(format_meta_list)
        else:
            if self.mttype == 'snv':
                tri = self.fasta.fetch(reference=self.chrom, start=max(0,self.pos-2), end=self.pos+1).upper()  # variant.POS is 1-based coordinate
                if self.ref in 'CT':
                    self.info_values['context'] = f'{tri}>{self.alt}'
                else:
                    self.info_values['context'] = f'{ "".join([ self.compl_base[x] for x in tri[::-1] ]) }>{ self.compl_base[self.alt] }'
                
                
    def addInfo_repeatRegion(self):
        info_meta_list = [
            {'ID':'repeatRegion', 'Type':'Integer', 'Number':1, 'Description':'If 1, current position is within a possible repeat region. 0 if not.'}
        ]
        format_meta_list = [
            {'ID':'snv_withinRepeat', 'Type':'Integer', 'Number':1, 'Description':'If 1, current snv variant is likely to be a mapping artifact. 0 if not.'}
        ]
        
        if self.makeMeta:
            self.info_meta.extend(info_meta_list)
            self.format_meta.extend(format_meta_list)
        else:
            if (not self.makeNull) and self.readsinfo.rplist != None:
                #if len(self.readsinfo.repeat_coord) >= 5:
                if len(self.readsinfo.relevant_repeats) > 0:
                    self.info_values['repeatRegion'] = 1
                else:
                    self.info_values['repeatRegion'] = 0
                    
                if self.readsinfo.snv_withinRepeat == True:
                    self.format_values['snv_withinRepeat'] = 1
                else:
                    self.format_values['snv_withinRepeat'] = 0
                
                
    def addInfo_satellite(self):
        info_meta_list = [
            {'ID':'satellite', 'Type':'Integer', 'Number':1, 'Description':'If 1, current position is within a satellite area. Areas overlapping with Gencode features are not included here.'}
        ]
        format_meta_list = []
        
        if self.makeMeta == True:
            self.info_meta.extend(info_meta_list)
            self.format_meta.extend(format_meta_list)
        else:
            self.info_values['satellite'] = 1 if self.readsinfo.satellite else 0
                
                
def get_mttype(ref, alt):
    if len(ref) == len(alt):
        return 'snv' if len(ref) == 1 else 'mnv'
    else:
        if len(ref) != 1 and len(alt) != 1:
            return 'cindel'
        elif len(ref) == 1:
            return 'ins'
        else:
            return 'del'
                
                
def write_results(**args):
    
    def getMetaLines(level):
        INFO = makeINFO(makeMeta = True, level = level)
        return INFO.info_meta, INFO.format_meta
    
    def write_modified_variants(vcf, writer, args, new_info_meta, new_format_meta):
        # loop over variant lines and add new INFO fields
        #NR = 0
        for variant in vcf:
            #NR += 1
            #if NR % 100 == 0:
            #    print(f'{NR}th line processed', end = '\r', file = sys.stderr, flush = True)

            alt_list = copy.deepcopy(variant.ALT)
            for alt in alt_list:
                new_args = dict()
                new_args['chrom'], new_args['pos'], new_args['ref'], new_args['alt'], new_args['REF_start'], new_args['REF_end'] = variant.CHROM, variant.POS, variant.REF, alt, variant.start, variant.end
                
                # filter non-ACGT ref or alt alleles (does not write the current line to the output)
                if re.match('^[ACGT]+$', new_args['alt']) == None or re.match('^[ACGT]+$', new_args['ref']) == None:
                    continue
                
                ### for debugging ###
                if args['debug'] == True:
                    global readsinfo
                    global INFO
                #####################
    
                # set mttype
                new_args['mttype'] = get_mttype(new_args['ref'], new_args['alt'])
                
                # get readsinfo and INFO
                readsinfo = makeReadsinfo(**new_args, **args)
                INFO = makeINFO(makeMeta = False, readsinfo = readsinfo, **new_args, **args)
                
                # set 'RdInfoAnnot' flag
                INFO.format_values['RdInfoAnnot'] = args['level']
                    
                # add new INFO and FORMAT values to the variant line and write it to output
                for field, val in INFO.info_values.items():
                    variant.INFO.__setitem__(field, val)
                for field, val in INFO.format_values.items():
                    old_value = variant.format(field)
                    
                    if isinstance(old_value, np.ndarray) == False:
                        for dic in new_format_meta:
                            if dic['ID'] == field:
                                field_type = dic['Type'] ; break
                        if field_type == 'Integer':
                            new_value = np.array( [[format_MISSING_int]] * len(vcf.samples) )
                        elif field_type == 'Float':
                            new_value = np.array( [[format_MISSING_float]] * len(vcf.samples) )
                    else:
                        new_value = old_value
                    
                    new_value[ args['sampleID_idx'], 0 ] = val
                        
                    variant.set_format(field, new_value)

                # edit ALT of variant object
                variant.ALT = [alt]
                    
                writer.write_record(variant)

                gc.collect()
            
        
        
    def main():
        import subprocess, tempfile

        vcf = cyvcf2.VCF(args['vcf_path'])
        
        # check if specified sample name is present in the vcf file ; create a tmp file which has an additional empty sample column
        flag_sampleID_absent = False
        
        if args['sampleID'] not in vcf.samples:
            flag_sampleID_absent = True
            
            # make a temporary file
            tf = tempfile.NamedTemporaryFile(dir = os.path.dirname(os.path.abspath(args['vcf_path'])), prefix = os.path.basename(args['vcf_path']) + '.', delete = False)
            tmpfile_path = tf.name
            tf.close()
            
            # write temporary vcf content to the temporary file
            sample_adder = '/home/users/pjh/scripts/vcf_editing/add_empty_sample.awk'
                # arg1 : input vcf file
                # arg2 : sample ID to add
                # arg3 : output file path
                # output file file is of compressed bcf format.
            subprocess.run([ sample_adder, args['vcf_path'], args['sampleID'], tmpfile_path ])
            
            # close old vcf_handler and make a new one
            vcf.close()
            vcf = cyvcf2.VCF(tmpfile_path)
            
        # get the index of sample ID
        args['sampleID_idx'] = vcf.samples.index(args['sampleID'])
        
        # add new INFO and FORMAT meta lines
        new_info_meta, new_format_meta = getMetaLines(args['level'])
        for dic in new_info_meta:
            vcf.add_info_to_header(dic)
        for dic in new_format_meta:
            vcf.add_format_to_header(dic)
        
        # add a meta line for the flag showing that readinfo annotation is done
        vcf.add_format_to_header( {'ID':'RdInfoAnnot', 'Type':'Integer', 'Number':1, 'Description':'Flag to show whether readinfo annotation was performed. 1: readCountOnly ; 2: everything'} )
            
        # make another temporary file
        tmpfile2 = tempfile.NamedTemporaryFile(dir = os.path.dirname(os.path.abspath(args['vcf_path'])))
        tmpfile2_path = tmpfile2.name
        tmpfile2.close()
            
        # write resultant header lines to output file
        writer = cyvcf2.Writer(tmpfile2_path, vcf, 'wb') 
            # mode 'wbu': uncompressed BCF, 'wb': compressed BCF, 'wz': compressed VCF, 'w': uncompressed VCF
        writer.write_header()
        
        # add new information to each variant line and write to output file
        write_modified_variants(vcf, writer, args, new_info_meta, new_format_meta)
            
        writer.close()
        vcf.close()
            
        # cleanup of temp file with samplename addition
        if flag_sampleID_absent == True:
            os.remove(tmpfile_path)
            try:
                os.remove(tmpfile_path + '.csi')
            except:
                pass
                
        # remove 'nan' letters
        subprocess.run([nantodot_path, tmpfile2_path, '-o', args['ofile_path'], '-O', args['ofile_format_rawarg']])
        os.remove(tmpfile2_path)
        try:
            os.remove(tmpfile2_path + '.csi')
        except:
            pass
        
    main()
        
        
def argument_parsing(cmdargs = None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', required = True)
    parser.add_argument('-s', required = True)
    #parser.add_argument('--readlen', type = int, required = True)
    parser.add_argument('--level')
    parser.add_argument('-o')
    parser.add_argument('-O', default = 'z')
    parser.add_argument('-p', type = int, default = 1)
    parser.add_argument('--log', action = 'store_true')
    parser.add_argument('--limBQ', type = int, default = 20)
    parser.add_argument('--limMQ', type = int, default = 20)
    parser.add_argument('--MMfrac', type = float, default = 0.1)
    parser.add_argument('--flankLen', type = int, default = 1)
    parser.add_argument('--debug', action = 'store_true')
    parser.add_argument('vcf')
    parser.add_argument('bam')
    raw_args = parser.parse_args(cmdargs)
    
    # sanity checks
    #if re.match('^[A-Za-z_][0-9A-Za-z_.]*$', raw_args.s) == None:
    #    printErr('sample ID must match pattern ^[A-Za-z_][0-9A-Za-z_.]*$') ; sys.exit()
    if raw_args.level not in [None, 'readCountOnly', 'everything']:
        printErr('possible values of -l are: "readCountOnly", "everything"') ; sys.exit(1)
    if raw_args.O not in 'vzub':
        printErr('possible values of -O are v,z,u,b') ; sys.exit(1)
        
    args = dict()
    
    args['fopt'] = raw_args.f
    
    if raw_args.f == '19':
        args['fasta'] = pysam.FastaFile(fasta_path_hg19)
    elif raw_args.f == '38':
        args['fasta'] = pysam.FastaFile(fasta_path_hg38)
    else:
        printErr('-f option must be "19" or "38".') ; sys.exit(1)
        #args['fasta'] = pysam.FastaFile(raw_args.f)
        
    args['sampleID']   = raw_args.s # sample id
    #raw_args['readlen']    = args.readlen
    
    level_dic = {'readCountOnly':1, 'everything':2, None:2}
    args['level']      = level_dic[raw_args.level]
    
    args['vcf_path']   = raw_args.vcf
    args['bam']        = pysam.AlignmentFile(raw_args.bam)
    
    if raw_args.o == None:
        args['ofile_path'] = '-'
    else:
        args['ofile_path'] = raw_args.o
        
    args['ofile_format_rawarg'] = raw_args.O
    ofile_format_table = { 'v':'w', 'z':'wz', 'u':'wbu', 'b':'wb' } # keys: used in bcftools ; values: used in cyvcf2.Writer
    args['ofile_format_cyvcf2'] = ofile_format_table[raw_args.O]
    
    args['p'] = raw_args.p
        
    if raw_args.log == True:
        args['log_path'] = f'{args["vcf_path"]}.readannot.log'
    else:
        args['log_path'] = '/dev/stderr'
        
    args['limBQ'] = raw_args.limBQ
    args['limMQ'] = raw_args.limMQ
    args['MMfrac'] = raw_args.MMfrac
    args['flankLen'] = raw_args.flankLen
    args['debug'] = raw_args.debug
    
    return args


def helpMessage():
    if sys.argv[1] == '-h':
        printErr('Usage: <program> -f <reference fasta> -s <sample name> [--level readCountOnly|everything] [-o <outfile path>] [-O v|z|u|b] [--limBQ <limBQ>] [--limMQ <limMQ>] [--MMfrac <MMfrac>] [--flankLen <flanking length>] <vcf file> <bam file>')
        printErr('')
        printErr('-f : Reference genome version. Must be "19" or "38". Fasta file is used for getting trinucleotide context.')
        printErr('-s : Sample name. This will be added as a prefix to all INFO field names. Only alphanumeric characters and underbar allowed.')
        #printErr('--readlen : read length. Must be an integer.')
        printErr('--level : level of annotation. Default is full annotation. If this value is "readCountOnly", only annotate numbers of ref/alt/other/total reads and vaf.')
        printErr('-o : Path of result output file. Default: stdout')
        printErr('-O : Format of output file. One of v,z,u,b (identical to bcftools). Default: "z"')
            # mode 'wbu': uncompressed BCF, 'wb': compressed BCF, 'wz': compressed VCF, 'w': uncompressed VCF
        #printErr('-p : Number of subprocesses for parallelization. Default: 1')
        #printErr('--log : If set, writes log for nonfatal errors to <vcf file>.readannot.log')
        printErr('--limBQ : reads with BQ of target base below this value will be treated as non-supporting reads. Default: 10')
        printErr('--limMQ : reads with MQ below this value will be treated as non-supporting reads. Default: 10')
        printErr('--MMfrac : reads with mismatch (excluding the variant itself) more than <read length including soft-clipped bases>*<MMfrac> will be treated as non-supporting reads. Default: 0.1')
        printErr('--flankLen : <flankLen> bases on either left or right of the variant position must be cigar M to be considered as ref/alt supporting. Default: 1')
        printErr('only single bam file is allowed.')
        sys.exit(1)
    
    
def main(cmdargs = None):
    helpMessage()
    
    # set some constants 
    add_sample_path = '/home/users/pjh/scripts/vcf_editing/add_empty_sample.awk'
    
    compl_base = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    
    stackDict_path = '/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt'
    stackDict = {k:v for k, v in np.genfromtxt('/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt', dtype=int)}

    # parse commandline arguments
    args = argument_parsing(cmdargs)
    
    # add constants to 'args' dictionary
    args['compl_base'] = compl_base
    args['stackDict'] = stackDict
    args['add_sample_path'] = add_sample_path
    
    write_results(**args)
    
if __name__ == '__main__':
    main()

# args={}
# add_sample_path = '/home/users/pjh/scripts/vcf_editing/add_empty_sample.awk'
# compl_base = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
# stackDict_path = '/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt'
# stackDict = {k:v for k, v in np.genfromtxt('/home/users/pjh/scripts/annotation/short_variants/max_readstack.txt', dtype=int)}
# args['compl_base'] = compl_base
# args['stackDict'] = stackDict
# args['add_sample_path'] = add_sample_path