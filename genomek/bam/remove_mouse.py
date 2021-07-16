import pysam
import argparse


def main(args):
    input_bam_path = args.input
    output_bam_path = args.output
    chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]

    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
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
            
            assert read_1.query_name == read_2.query_name
            
            if (read_1.reference_name and read_1.reference_name.startswith("Mus")) or (read_2.reference_name and read_2.reference_name.startswith("Mus")):
                read_1 = None
                read_2 = None
                continue
            
            output_bam.write(read_1)
            output_bam.write(read_2)
            
            read_1 = None
            read_2 = None
        
        i += 1
            
        if i % 1_000_000 == 0:
            print(i, " reads done")

    print(i, "all done")

    output_bam.close()
    input_bam.close()


if __name__ == '__main__':
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='Human bam editing')
    parser.add_argument("-i","--input", type=str, required = True, help='(Required) input bam file')
    parser.add_argument("-o", "--output", type=str, required = True, help='(Required) output bam file')
    args=parser.parse_args()
    main(args)