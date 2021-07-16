# -*- coding: UTF-8 -*-

'''

    MT variant filtering
    - input: .tsv (with POS, REF, ALT을 포함)
    - output: .tsv (with baminfo annotation)

    Version history
    - 2.0: MQ, BQ 를 통과해야만 n_read 증가하도록 수정, duplicated read 제거
    - 3.0: MT 뿐만 아니라 모든 chromosome에 적용할 수 있도록 확장 
    - 4.0: df.loc이 CPU를 너무 많이 점유하는 문제 해결. 모든 계산값을 list에 넣은 후 마지막에 dataframe에 추가하는 방식으로 변경함

'''


import pandas as pd, numpy as np
import pysam
import statistics


def annotate_polyns(chrom, pos, reference, max_repeat_unit_size=5, window_size=5):
    '''
    function used to annotate indel repeat contexts by chnam
    '''
    chrom, pos = str(chrom), int(pos)
    ref = pysam.FastaFile(reference)

    up_base_list=[]
    up_n_list=[]
    down_base_list=[]
    down_n_list=[]
    
    for s in range(1,max_repeat_unit_size+1,1):
        up_base = ref.fetch(chrom,pos-s-1,pos-1)
        up_base_list.append(up_base)
        up_n=1
        while True:
            if ref.fetch(chrom,pos-up_n*s-s-1,pos-up_n*s-1) == up_base: up_n+=1
            else: break
        up_n_list.append(up_n)
        
    if len(set(up_n_list)) == 1:
        max_up_base = up_base_list[0]
        max_up_n = up_n_list[0]
    else:
        max_up_n = max(up_n_list)
        max_up_base = up_base_list[up_n_list.index(max(up_n_list))]
    
    for s in range(1,max_repeat_unit_size+1,1):
        down_base = ref.fetch(chrom,pos,pos+s)
        down_base_list.append(down_base)
        down_n=1
        while True:
            if ref.fetch(chrom,pos+down_n*s,pos+down_n*s+s) == down_base: down_n+=1
            else: break
        down_n_list.append(down_n)
        
    if len(set(down_n_list)) == 1:
        max_down_base = down_base_list[0]
        max_down_n = down_n_list[0]
    else:
        max_down_n = max(down_n_list)
        max_down_base = down_base_list[down_n_list.index(max(down_n_list))]
    
    window           = ref.fetch(chrom,pos-window_size-1,pos+window_size)
    window_base_type = "".join(sorted(list(set(window))))
    window_base_n    = len(window_base_type)
    return f'{max_up_base}:{max_up_n}:{max_down_base}:{max_down_n}:{window_base_type}:{window_base_n}'


#  homopolymer check 위한 rle encoder
def encode_rle(input_string):

    count = 1
    prev = ''
    lst = []

    for character in input_string:

        if character != prev:

            if prev:
                entry = (prev, count)
                lst.append(entry)

            count = 1
            prev = character

        else:
            count += 1

    else:

        try:

            entry = (character, count)
            lst.append(entry)
            return lst

        except Exception as e:

            print("Exception encountered {e}".format(e=e))
            return e


# main function
def annotate_baminfo(df, bam, ref, bq=10, mq=10, rw=5, iw=5, suffix=""):

    # df: POS, REF, ALT column이 있는 panda dataframe
    # bam: bam file (pysam object)
    # ref: reference file (pysam object)
    # bq: base quality threshold
    # mq: mapping quality threshold 
    # rw: window for checking homopolymer 

    list_mean_read_position_from_5p = [None] * len(df)
    list_mean_read_position_from_lt = [None] * len(df)
    list_strandedness               = [None] * len(df)
    list_variant_reads              = [None] * len(df)
    list_variant_reads_forward      = [None] * len(df)
    list_mean_MQ_ref                = [None] * len(df)
    list_mean_MQ_alt                = [None] * len(df)
    list_mean_read_length_ref       = [None] * len(df)
    list_mean_read_length_alt       = [None] * len(df)
    list_reads_with_flanking_indel  = [None] * len(df)
    list_n_seq_changes              = [None] * len(df)
    list_n_base_type_in_rle_window  = [None] * len(df)
    list_n_max_homopolyer_length    = [None] * len(df)
    list_n_reads_in_normal_context  = [None] * len(df)
    list_total_depth                = [None] * len(df)
    list_VAF                        = [None] * len(df)
    list_mean_NM_ref                = [None] * len(df)
    list_mean_NM_alt                = [None] * len(df)
    list_median_NM_ref              = [None] * len(df)
    list_median_NM_alt              = [None] * len(df)
    list_poorly_mapped_reads        = [None] * len(df)

    for i, row in df.iterrows():

        print(row.CHROM, ":", row.POS, "\t", row.REF, ">", row.ALT, sep="")

        if row.REF[0] not in ['A', 'G', 'C', 'T']:
            print("\t- skip")
            continue  # ex) MT:3107 일 경우...

        read_position_from_5p      = list() # position in supporting reads, relative to read length
        read_position_from_lt      = list() # position in supporting reads, relative to lt most position
        n_variant_reads            = 0 # variant read 개수
        n_forward_reads            = 0 # variant read 중 forward read 개수
        n_reads                    = 0 # 전체 read 개수
        mapping_quality_ref_reads  = list()
        mapping_quality_alt_reads  = list()
        read_length_ref_reads      = list()
        read_length_alt_reads      = list()
        flanking_indel             = 0
        n_reads_in_normal_context  = 0
        NM_ref_reads               = list()
        NM_alt_reads               = list()
        n_poorly_mapped_reads      = 0

        # row.CHROM이 int로 되어있으면 강제로 string으로 변환시킴 --> fetch에서 에러가 나지 않도록 하기 위함 
        row.CHROM = str(row.CHROM)

        # flanking 5b에 homopolymer가 있는지 체크
        # reference의 양쪽 끝을 벗어나면 에러가 나므로, 양쪽 끝을 벗어나는 경우 각각 0, chromsome length로 설정
        rle = encode_rle(ref.fetch(row.CHROM, 
            (row.POS - 1 - rw) if (row.POS - 1 - rw) >= 0 else 0, 
            (row.POS + rw) if (row.POS + rw) <= bam.get_reference_length(row.CHROM) else bam.get_reference_length(row.CHROM)))

        n_seq_changes = len(rle)  # base change가 몇번 일어나는가?
        n_base_types_in_rle_window = len(set(x[0] for x in rle))   # base 종류의 개수
        n_max_homopolyer_length = max([x[1] for x in rle])

        # ALT가 *이면 deletion으로 처리하도록 한다. 
        if row.ALT is "*":
            row.POS = row.POS - 1
            row.ALT = ref.fetch(row.CHROM, row.POS - 2, row.POS - 1)  # * 직전 base의 REF 값을 얻음 
            row.REF = row.ALT + row.REF  # deletion으로 처리할 것이므로 REF는 직전 base + * 위치의 base가 됨            

        # snp
        if len(row.REF) == len(row.ALT):  

            for read in bam.fetch(row.CHROM, row.POS - 1, row.POS):

                # duplicated read 이면 건너뜀
                if read.is_duplicate: continue 

                # mapping quality check
                if read.mapping_quality < mq: 
                    n_poorly_mapped_reads = n_poorly_mapped_reads + 1
                    continue

                # 해당 위치의 base를 읽어옴
                aligned_pairs = read.get_aligned_pairs()
                try:
                    element = list(filter(lambda x: x[0] is not None and x[1] is not None and x[1] == (row.POS - 1), aligned_pairs))[0]  
                    index   = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0]         
                except:
                    continue  # 해당하는 위치에 base가 없으면 건너뜀 (deletion, soft clip 등..)

                # MQ 를 통과하고 해당위치에 base가 있으면 n_read 추가 
                n_reads += 1 # MQ 통과해야만 n_read에 추가 

                b = read.query_sequence[aligned_pairs[index][0]]
                if b is not row.ALT: # variant read가 아니면 건너뜀
                    # 건너뛰기 전에 reference read 이면 mapping quality, read length, NM tag 기록
                    if b is row.REF:
                        mapping_quality_ref_reads.append(read.mapping_quality)
                        read_length_ref_reads.append(read.query_length)
                        NM_ref_reads.append(read.get_tag("NM"))
                    continue 

                # base quality check, 만약 threshold 를 넘지 못하면 variant read가 아닌것으로 간주, 건너뛴다. (n_read는 추가됨)
                if read.query_qualities[aligned_pairs[index][0]] < bq: continue                                    

                # read에서의 위치를 percent로 구함
                if not read.is_reverse:
                    read_position_from_5p.append(aligned_pairs[index][0]/read.query_length)
                else:
                    read_position_from_5p.append((read.query_length - aligned_pairs[index][0])/read.query_length)

                read_position_from_lt.append(aligned_pairs[index][0]/read.query_length)

                # variant가 normal context내에 있는지 check
                try:
                    if read.query_sequence[aligned_pairs[index][0]-1] == ref.fetch(row.CHROM, row.POS - 2, row.POS - 1) and read.query_sequence[aligned_pairs[index][0]+1] == ref.fetch(row.CHROM, row.POS, row.POS + 1): 
                        n_reads_in_normal_context += 1
                except IndexError: # variant position이 read의 끝부분이거나 맨 앞 부분이면 list indexing할때 에러가 남. 이경우 한쪽만 확인되는 상태이므로 n_reads_in_normal_context를 증가하지 않음
                    pass

                # strandedness
                n_variant_reads += 1
                if not read.is_reverse: n_forward_reads += 1

                # variant read의 mapping quality 기록
                mapping_quality_alt_reads.append(read.mapping_quality)

                # variant read의 read length 기록
                read_length_alt_reads.append(read.query_length)

                # 5bp 이내에 indel이 있는지 여부 check
                flanking_pairs = aligned_pairs[(index - iw):(index + iw)]
                if None in [x[0] for x in flanking_pairs] or None in [x[1] for x in flanking_pairs]: flanking_indel += 1

                # NM tag
                NM_alt_reads.append(read.get_tag("NM"))


        # insertion
        elif len(row.REF) < len(row.ALT):
            
            for read in bam.fetch(row.CHROM, row.POS - 1, row.POS):

                # duplicated read 이면 건너뜀
                if read.is_duplicate: continue                 

                # mapping quality check
                if read.mapping_quality < mq: 
                    n_poorly_mapped_reads = n_poorly_mapped_reads + 1 
                    continue

                # insertion 직전 위치
                aligned_pairs = read.get_aligned_pairs()
                try:                    
                    element = list(filter(lambda x: x[0] is not None and x[1] is not None and x[1] == (row.POS - 1), aligned_pairs))[0]  
                    index   = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0]                    
                except:
                    continue # base가 없으면 건너뜀

                n_reads += 1 # total read count 하나 증가

                b = read.query_sequence[aligned_pairs[index][0]]
                if b is not row.REF: # insertion 직전 위치의 base는 REF와 같아야 한다.    
                    continue

                # insertion site 직후부터 insertion이 있어야함. 그렇지 않으면 mapping quality 기록하고 건너뜀_pairs):
                if (index + 1) == len(aligned_pairs) or aligned_pairs[index + 1][1] is not None:
                    mapping_quality_ref_reads.append(read.mapping_quality)
                    read_length_ref_reads.append(read.query_length)
                    NM_ref_reads.append(read.get_tag("NM"))
                    continue
                
                # 해당 부위의 query가 reference에 align 되지 않았는지 확인
                if [x[1] for x in aligned_pairs[(index + 1):(index + len(row.ALT))]] != [None] * (len(row.ALT) - 1):                
                    # insertion 길이가 다르다고해서 reference read로 간주하지는 않는다. 
                    continue                
                
                # query의 insertion 된 seq 직후로는 reference에 align 되어야 한다
                try:
                    if aligned_pairs[index + len(row.ALT)][1] is None:
                        continue
                except IndexError: # read의 끝이라면 evidence가 충분하지 않은 것으로 보고 건너뜀
                    continue

                # variant가 normal context내에 있는지 check (variant 직전은 위에서 check 했으므로 variant 이후 base만 check)
                try:
                    if read.query_sequence[aligned_pairs[index][0]+len(row.ALT)] == ref.fetch(row.CHROM, row.POS, row.POS + 1):
                        n_reads_in_normal_context += 1
                except IndexError: 
                    pass

                # query에서의 위치를 percent로 구함
                if not read.is_reverse:                    
                    read_position_from_5p.append(aligned_pairs[index][0]/read.query_length)
                else:
                    read_position_from_5p.append((read.query_length - aligned_pairs[index][0])/read.query_length)

                read_position_from_lt.append(aligned_pairs[index][0]/read.query_length)

                # strandedness
                n_variant_reads += 1
                if not read.is_reverse: n_forward_reads += 1

                # variant read의 mapping quality 기록
                mapping_quality_alt_reads.append(read.mapping_quality)

                # variant read의 read length 기록
                read_length_alt_reads.append(read.query_length)

                # 5bp 이내에 indel이 있는지 여부 check
                flanking_pairs = aligned_pairs[(index-5):index] + aligned_pairs[(index + len(row.ALT)):(index + len(row.ALT) + 5)]
                if None in [x[0] for x in flanking_pairs] or None in [x[1] for x in flanking_pairs]: flanking_indel += 1

                # NM tag
                NM_alt_reads.append(read.get_tag("NM"))

        # deletion
        elif len(row.REF) > len(row.ALT): 

            for read in bam.fetch(row.CHROM, row.POS - 1, row.POS):

                # duplicated read 이면 건너뜀
                if read.is_duplicate: continue                 

                # mapping quality check
                if read.mapping_quality < mq: 
                    n_poorly_mapped_reads = n_poorly_mapped_reads + 1 
                    continue

                # deletion 직전 위치
                aligned_pairs = read.get_aligned_pairs()                
                try:
                    element = list(filter(lambda x: x[0] is not None and x[1] is not None and x[1] == (row.POS - 1), aligned_pairs))[0]  
                    index   = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0]                                        
                except:                    
                    continue  # base가 없으면 건너뜀

                n_reads += 1

                b = read.query_sequence[aligned_pairs[index][0]]
                if b is not row.REF[0]:  # REF가 다르면 건너뜀
                    continue

                # ref 직후부터 deletion이 시작되어야 함(그렇지 않으면 deletion이 없는 read임). 그렇지 않으면 mapping quality, read length 기록하고 건너뜀
                if (index + 1) == len(aligned_pairs) or aligned_pairs[index + 1][0] is not None:
                    mapping_quality_ref_reads.append(read.mapping_quality)
                    read_length_ref_reads.append(read.query_length)
                    NM_ref_reads.append(read.get_tag("NM"))
                    continue

                # 해당 부위의 query가 비어있는지 확인
                if [x[0] for x in aligned_pairs[(index + 1):(index + len(row.REF))]] != [None] * (len(row.REF) - 1):
                    # 비어있는 길이가 다르다고 해서 reference read로 간주하지는 않는다
                    continue

                # query의 deletion 된 seq 직후로는 reference에 align 되어야 한다(그렇지 않으면 더 많이 deletion되어있는 read임)
                try:
                    if aligned_pairs[index + len(row.REF)][0] is None:                        
                        continue
                except IndexError: # deletion 직후의 base를 볼수 없으면 (read의 끝이라면) evidence가 충분하지 않으므로 건너뜀
                    continue

                # variant가 normal context내에 있는지 check (variant 직전은 위에서 check 했으므로 variant 이후 base만 check)                
                try:
                    if read.query_sequence[aligned_pairs[index][0]+1] == ref.fetch(row.CHROM, row.POS + len(row.REF) - 1, row.POS + len(row.REF)):
                        n_reads_in_normal_context += 1
                except IndexError:
                    pass

                # query에서의 위치를 percent로 구함
                if not read.is_reverse:
                    read_position_from_5p.append(aligned_pairs[index][0]/read.query_length)
                else:
                    read_position_from_5p.append((read.query_length - aligned_pairs[index][0])/read.query_length)

                read_position_from_lt.append(aligned_pairs[index][0]/read.query_length)

                # strandedness
                n_variant_reads += 1
                if not read.is_reverse: n_forward_reads += 1

                # variant read의 mapping quality 기록
                mapping_quality_alt_reads.append(read.mapping_quality)

                # variant read의 read length 기록
                read_length_alt_reads.append(read.query_length)

                # 5bp 이내에 indel이 있는지 여부 check
                flanking_pairs = aligned_pairs[(index - 5):index] + aligned_pairs[(index + len(row.REF)):(index + len(row.REF) + 5)]
                if None in [x[0] for x in flanking_pairs] or None in [x[1] for x in flanking_pairs]: flanking_indel += 1

                # NM tag
                NM_alt_reads.append(read.get_tag("NM"))


        # coverage가 전혀 없다면 모두 np.nan으로 처리
        if (n_reads + n_poorly_mapped_reads) is 0:
            mean_read_position_from_5p   = np.nan
            mean_read_position_from_lt   = np.nan
            strandedness                 = np.nan
            mean_ref_MQ                  = np.nan
            mean_alt_MQ                  = np.nan
            mean_ref_read_length         = np.nan 
            mean_alt_read_length         = np.nan
            flanking_indel               = np.nan
            n_seq_changes                = np.nan
            n_base_types_in_rle_window   = np.nan
            n_max_homopolyer_length      = np.nan
            n_variant_reads              = np.nan
            n_forward_reads              = np.nan
            n_reads_in_normal_context    = np.nan
            n_reads                      = np.nan
            VAF                          = np.nan
            mean_NM_ref                  = np.nan
            mean_NM_alt                  = np.nan
            median_NM_ref                = np.nan
            median_NM_alt                = np.nan
            n_poorly_mapped_reads        = np.nan
        else:

            if read_position_from_5p:
                mean_read_position_from_5p = sum(read_position_from_5p) / len(read_position_from_5p)
            else:
                mean_read_position_from_5p = np.nan

            if read_position_from_lt:
                mean_read_position_from_lt = sum(read_position_from_lt) / len(read_position_from_lt)
            else:
                mean_read_position_from_lt = np.nan

            if mapping_quality_ref_reads:
                mean_ref_MQ = sum(mapping_quality_ref_reads) / len(mapping_quality_ref_reads)
            else: # ref read가 하나도 없을 경우..
                mean_ref_MQ = np.nan

            if mapping_quality_alt_reads:
                mean_alt_MQ = sum(mapping_quality_alt_reads) / len(mapping_quality_alt_reads)
            else: # alt read가 하나도 없을 경우..
                mean_alt_MQ = np.nan

            if n_variant_reads:
                strandedness = n_forward_reads / n_variant_reads
            else:
                strandedness = np.nan

            if read_length_ref_reads:
                mean_ref_read_length = sum(read_length_ref_reads) / len(read_length_ref_reads)
            else:
                mean_ref_read_length = np.nan

            if read_length_alt_reads:
                mean_alt_read_length = sum(read_length_alt_reads) / len(read_length_alt_reads)
            else:
                mean_alt_read_length = np.nan

            if NM_alt_reads:
                mean_NM_alt   = sum(NM_alt_reads) / len(NM_alt_reads)
                median_NM_alt = statistics.median(NM_alt_reads)
            else:
                mean_NM_alt   = np.nan
                median_NM_alt = np.nan

            if NM_ref_reads:
                mean_NM_ref   = sum(NM_ref_reads) / len(NM_ref_reads)
                median_NM_ref = statistics.median(NM_ref_reads)
            else:
                mean_NM_ref   = np.nan
                median_NM_ref = np.nan

            if n_reads:
                VAF = n_variant_reads / n_reads * 100
            else:
                VAF = np.nan

        list_mean_read_position_from_5p[i] = mean_read_position_from_5p
        list_mean_read_position_from_lt[i] = mean_read_position_from_lt
        list_strandedness[i]               = strandedness
        list_mean_MQ_ref[i]                = mean_ref_MQ
        list_mean_MQ_alt[i]                = mean_alt_MQ
        list_mean_read_length_ref[i]       = mean_ref_read_length
        list_mean_read_length_alt[i]       = mean_alt_read_length
        list_reads_with_flanking_indel[i]  = flanking_indel
        list_n_seq_changes[i]              = n_seq_changes
        list_n_base_type_in_rle_window[i]  = n_base_types_in_rle_window
        list_n_max_homopolyer_length[i]    = n_max_homopolyer_length
        list_variant_reads[i]              = n_variant_reads
        list_variant_reads_forward[i]      = n_forward_reads
        list_n_reads_in_normal_context[i]  = n_reads_in_normal_context
        list_total_depth[i]                = n_reads
        list_VAF[i]                        = VAF
        list_mean_NM_ref[i]                = mean_NM_ref
        list_mean_NM_alt[i]                = mean_NM_alt
        list_median_NM_ref[i]              = median_NM_ref
        list_median_NM_alt[i]              = median_NM_alt
        list_poorly_mapped_reads[i]        = n_poorly_mapped_reads

        print("\t- VAF:                             : {:.3f}%".format(VAF))
        print("\t- mean position from 5p            : {:.3f}".format(mean_read_position_from_5p))
        print("\t- mean position from lt            : {:.3f}".format(mean_read_position_from_lt))
        print("\t- total depth                      : {}".format(n_reads))
        print("\t- variant reads                    : {}".format(n_variant_reads))
        print("\t- variant reads in normal context  : {}".format(n_reads_in_normal_context))
        print("\t- strandedness                     : {:.3f}".format(strandedness))
        print("\t- MQ difference                    : {:.3f}".format(mean_alt_MQ-mean_ref_MQ))
        print("\t- read length difference           : {:.3f}".format(mean_alt_read_length - mean_ref_read_length))
        print("\t- reads with flanking indel        : {}". format(flanking_indel))
        print("\t- number of sequence changes       : {}".format(n_seq_changes))
        print("\t- number of base types             : {}".format(n_base_types_in_rle_window))
        print("\t- number of max homopolyer length  : {}".format(n_max_homopolyer_length))
        print("\t- mean NM in ref reads             : {}".format(mean_NM_ref))
        print("\t- mean NM in alt reads             : {}".format(mean_NM_alt))
        print("\t- median NM in ref reads           : {}".format(median_NM_ref))
        print("\t- median NM in alt reads           : {}".format(median_NM_alt))
        print("\t- poorly mapped reads              : {}".format(n_poorly_mapped_reads))


    df[f'mean_read_position_from_5p{suffix}']  = list_mean_read_position_from_5p
    df[f'mean_read_position_from_lt{suffix}']  = list_mean_read_position_from_lt
    df[f'strandedness{suffix}']                = list_strandedness
    df[f'mean_MQ_ref{suffix}']                 = list_mean_MQ_ref
    df[f'mean_MQ_alt{suffix}']                 = list_mean_MQ_alt
    df[f'mean_read_length_ref{suffix}']        = list_mean_read_length_ref
    df[f'mean_read_length_alt{suffix}']        = list_mean_read_length_alt
    df[f'reads_with_flanking_indel{suffix}']   = list_reads_with_flanking_indel
    df[f'n_seq_changes{suffix}']               = list_n_seq_changes
    df[f'n_base_type_in_rle_window{suffix}']   = list_n_base_type_in_rle_window
    df[f'n_max_homopolyer_length{suffix}']     = list_n_max_homopolyer_length
    df[f'variant_reads{suffix}']               = list_variant_reads
    df[f'variant_reads_forward{suffix}']       = list_variant_reads_forward
    df[f'n_reads_in_normal_context{suffix}']   = list_n_reads_in_normal_context
    df[f'total_depth{suffix}']                 = list_total_depth
    df[f'VAF{suffix}']                         = list_VAF
    df[f'mean_NM_ref{suffix}']                 = list_mean_NM_ref
    df[f'mean_NM_alt{suffix}']                 = list_mean_NM_alt
    df[f'median_NM_ref{suffix}']               = list_median_NM_ref
    df[f'median_NM_alt{suffix}']               = list_median_NM_alt
    df[f'poorly_mapped_reads{suffix}']         = list_poorly_mapped_reads


    return(df)

if __name__ == '__main__':

    import argparse, os

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help="variant .tsv file (with CHROM, POS, REF, ALT)")
    parser.add_argument('-b', '--bam', help="corresponding bam file", required=True)
    parser.add_argument('-r', '--ref', help="reference", required=True)
    parser.add_argument('-o', '--out', help="output tsv", default=None)
    parser.add_argument('-bq', '--base_quality', help="base quality threshold", default=10)
    parser.add_argument('-mq', '--mapping_quality', help="mapping quality threshold", default=10)
    parser.add_argument('-rw', '--rle_window', help="window for checking homopolyer", default=5)
    parser.add_argument('-iw', '--indel_window', help="window for checking nearby indel", default=5)
    parser.add_argument('-suf', '--suffix', help="suffix to add at column names", default="")
    args = parser.parse_args()

    df  = annotate_baminfo(
        pd.read_csv(args.file, sep='\t', low_memory=False), 
        pysam.AlignmentFile(args.bam, threads=1), 
        pysam.FastaFile(args.ref), 
        int(args.base_quality), 
        int(args.mapping_quality), 
        int(args.rle_window), 
        int(args.indel_window),
        args.suffix)

    if args.out is not None:
        df.to_csv(args.out, index=False, sep='\t')
    else:
        df.to_csv(os.path.splitext(args.file)[0] + '_baminfo.tsv', index=False, sep='\t')