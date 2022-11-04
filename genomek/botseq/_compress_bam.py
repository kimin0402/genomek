import re
import sys
import copy
import array
import pysam 
import numpy as np
import pandas as pd
from collections import defaultdict
from collections import Counter


CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
CIGAR2CODE = dict([y, x] for x, y in enumerate("MIDNSHP=XB"))


def read_to_orientation(v: pysam.libcalignedsegment.AlignedSegment) -> str:
    if v.is_reverse: 
        ori = 'R'
    else: 
        ori = 'F'

    if v.is_read1:
        ori += '1'
    elif v.is_read2:
        ori += '2'
    else: # if read is not paired orientation is set to 0
        ori += '0'
        
    return ori


def cigartuples_to_number(cigar_tuples, start=True):
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
            return 151


def check_reads_overlap(read_1, read_2):
    read_1_start_pos = read_1.reference_start - cigartuples_to_number(read_1.cigartuples, start=True)
    read_1_end_pos = read_1.reference_start + cigartuples_to_number(read_1.cigartuples, start=False)
    read_2_start_pos = read_2.reference_start - cigartuples_to_number(read_2.cigartuples, start=True)
    read_2_end_pos = read_2.reference_start + cigartuples_to_number(read_2.cigartuples, start=False)
    if (read_1_start_pos <= read_2_start_pos <= read_1_end_pos) or (read_2_start_pos <= read_1_start_pos <= read_2_end_pos):
        return True
    else:
        return False


def bam_to_reads(input_bam: pysam.libcalignmentfile.AlignmentFile, chromosome: str=None):
    reads = defaultdict(list)
    for read in input_bam.fetch(contig=chromosome):
        if read.is_secondary:
            continue
        if read.is_unmapped:
            continue
        reads[read.query_name].append(read)
    return reads


def reads_to_umis(reads):
    umis = defaultdict(lambda: defaultdict(list))
    n, m, l = 0, 0, 0
    for k, v in reads.items():
        if len(v) == 1: # if the query name is unpaired (unique)
            l += 1
            read_of_interest = v[0]
            start_pos = read_of_interest.reference_start - cigartuples_to_number(read_of_interest.cigartuples, start=True)
            end_pos = read_of_interest.reference_start + cigartuples_to_number(read_of_interest.cigartuples, start=False)
            try:
                umis[(start_pos, end_pos)][read_to_orientation(v)].append(v)
            except:
                m += 1

        elif len(v) == 2: # if the query name has two reads (paired)
            v_0_start_pos = v[0].reference_start - cigartuples_to_number(v[0].cigartuples, start=True)
            v_0_end_pos = v[0].reference_start + cigartuples_to_number(v[0].cigartuples, start=False)
            v_1_start_pos = v[1].reference_start - cigartuples_to_number(v[1].cigartuples, start=True)
            v_1_end_pos = v[1].reference_start + cigartuples_to_number(v[1].cigartuples, start=False)

            if check_reads_overlap(v[0], v[1]):
                oreintation = ''.join(sorted([read_to_orientation(v[0]), read_to_orientation(v[1])]))

            else:
                if v_0_start_pos < v_1_start_pos:
                    oreintation = ''.join([read_to_orientation(v[0]), read_to_orientation(v[1])])
                elif v_1_start_pos < v_0_start_pos:
                    oreintation = ''.join([read_to_orientation(v[1]), read_to_orientation(v[0])])

            start_pos = min(v_0_start_pos, v_1_start_pos)
            end_pos = max(v_0_end_pos, v_1_end_pos)
            insert_size = abs(start_pos - end_pos)
            if insert_size == 0: continue
            try:
                umis[(start_pos, end_pos)][oreintation].append(v)
            except:
                n += 1

    print(f"Number of total singles: {l}", file=sys.stderr)
    print(f"Number of improperly oriented singles: {m}", file=sys.stderr)
    print(f"Number of improperly oriented pairs: {n}", file=sys.stderr)
    return umis


def seq_to_matrix(sequence, length=151):
    base_to_number = {'A':0, 'C':1, 'G':2, 'T':3}
    matrix = np.zeros((4, length))
    for j, i in enumerate(sequence):
        try:
            matrix[base_to_number[i], j] = 1
        except:
            continue
    return matrix


def matrix_to_seq(matrix):
    number_to_base = {0:'A', 1:'C', 2:'G', 3:'T'}
    seq = ""
    for j in range(matrix.shape[1]):
        try:
            seq += number_to_base[np.argmax(matrix[:,j])]
        except:
            continue
    return seq
            

def reads_to_read(reads: list, read_name: str) -> pysam.libcalignedsegment.AlignedSegment:
    '''
    Given a list of pysam alignment segments, return a compressed pysam alignment segment
    '''
    number_of_reads = len(reads)
    max_query_length = max(x.query_length for x in reads)
    for read in reads:
        if read.query_length == max_query_length:
            read_compressed = copy.deepcopy(read)
            break

    sequence_mtx = np.zeros((4, max_query_length))
    for read in reads:
        qualities = np.zeros(max_query_length)
        qualities[0:len(read.query_qualities)] = read.query_qualities
        sequence = read.query_sequence
        sequence_mtx += seq_to_matrix(sequence, length=max_query_length) * qualities

    qualities = array.array('B', (np.max(sequence_mtx, axis=0) / number_of_reads).astype(int))
    sequence = matrix_to_seq(np.where(sequence_mtx == sequence_mtx.max(axis=0, keepdims=True), 1, 0).astype(int))

    read_compressed.query_name = read_name
    read_compressed.query_sequence = sequence
    read_compressed.query_qualities = qualities

    return read_compressed







def collapse_umis_to_bam(k, v, read_orientation: str, output_bam: pysam.AlignmentFile):
    number_of_reads = len(v[read_orientation])
    if number_of_reads:
        read_name = ':'.join([str(x) for x in k]) + ':'+read_orientation+':' + str(number_of_reads)
        first_query_length = max([x[0].query_length for x in v[read_orientation]]) # v['F1R2'][0][0].query_length
        second_query_length = max([x[1].query_length for x in v[read_orientation]]) # v['F1R2'][0][1].query_length
        r1 = copy.deepcopy(v[read_orientation][0][0])
        r2 = copy.deepcopy(v[read_orientation][0][1])
        r1_sequence_mtx = np.zeros((4, first_query_length))
        r2_sequence_mtx = np.zeros((4, second_query_length))
        
        for reads_f1r2 in v[read_orientation]:
            r1_sequence = reads_f1r2[0].query_sequence
            r2_sequence = reads_f1r2[1].query_sequence
            r1_qualities = np.zeros(first_query_length)
            r2_qualities = np.zeros(second_query_length)
            r1_qualities[0:len(reads_f1r2[0].query_qualities)] = reads_f1r2[0].query_qualities
            r2_qualities[0:len(reads_f1r2[1].query_qualities)] = reads_f1r2[1].query_qualities
            r1_sequence_mtx += seq_to_matrix(r1_sequence, length=first_query_length) * r1_qualities
            r2_sequence_mtx += seq_to_matrix(r2_sequence, length=second_query_length) * r2_qualities
            
        r1_qualities = array.array('B', (np.max(r1_sequence_mtx, axis=0) / number_of_reads).astype(int))
        r2_qualities = array.array('B', (np.max(r2_sequence_mtx, axis=0) / number_of_reads).astype(int))
        r1_sequence = matrix_to_seq(np.where(r1_sequence_mtx == r1_sequence_mtx.max(axis=0, keepdims=True), 1, 0).astype(int))
        r2_sequence = matrix_to_seq(np.where(r2_sequence_mtx == r2_sequence_mtx.max(axis=0, keepdims=True), 1, 0).astype(int))

        r1.query_name = read_name
        r2.query_name = read_name
        r1.query_sequence = r1_sequence
        r2.query_sequence = r2_sequence
        r1.query_qualities = r1_qualities
        r2.query_qualities = r2_qualities
        r1.cigarstring = None
        r2.cigarstring = None

        v_collapse = [r1, r2]
        for v in v_collapse:
            output_bam.write(v)
        return
    else:
        return
            