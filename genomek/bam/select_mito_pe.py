# Select MT reads
# Single end reads

import sys
import pysam
import argparse

chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]


def argument_parsing():
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='Select mito reads')
    parser.add_argument("-i","--input", type=str, required=True, help='(Required) input bam file, or stdin to recieve STDIN')
    parser.add_argument("-o", "--output", type=str, required=True, help='(Required) output bam file, or stdout to output STDOUT')
    parser.add_argument("-c", "--chrom", type=str, default='MT', help='Mito chromosome symbol (MT, chrM ...) Default MT')
    args=parser.parse_args()
    args =vars(args)

    if args['input'] == 'stdin':
        args['input'] = '-'
    if args['output'] == 'stdout':
        args['output'] = '-'
    
    return args


def main(**args):
    input_bam_path = args['input']
    output_bam_path = args['output']
    chrom = args['chrom']

    input_bam = pysam.AlignmentFile(input_bam_path, "r")
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam)

    i = 0
    read_1 = None
    read_2 = None
    for read in input_bam:
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
            
        if not read_1:
            read_1 = read
            continue
        if read_1:
            read_2 = read
            
            try: 
                assert read_1.query_name == read_2.query_name
            except AssertionError:
                print(read_1.query_name, file=sys.stderr)
                print(read_2.query_name, file=sys.stderr)
                exit(1)

            
            if (read_1.reference_name and read_1.reference_name == chrom) and (read_2.reference_name and read_2.reference_name == chrom):
                output_bam.write(read_1)
                output_bam.write(read_2)
                read_1 = None
                read_2 = None
                continue
            
            
            read_1 = None
            read_2 = None
        
        i += 1
            
        if i % 1_000_000 == 0:
            print(i, " reads done", file=sys.stderr)

    print(i, "all done", file=sys.stderr)

    output_bam.close()
    input_bam.close()
    sys.stderr.close()

if __name__ == '__main__':
    args = argument_parsing()
    main(**args)