"""
March 11 2018 Jongsoo Yoon
filter mutect and strelka mutation calls
outputs three files 
those called in both <sampleName>_1n2.vcf
only in mutect<sampleName>_1-2.vcf
only in strelka <sampleName>_2-1.vcf

"""

import argparse
import cyvcf2
import re
import time
import datetime
import os
import shutil
import sys


# Prepare chromosome conversion dictionary for circos
chrDict = dict()
chromosomeList = [str(i) for i in range(1,23)] + ['X', 'Y']
for i in chromosomeList:
    chrDict.update({i: 'hs'+i})

# make sure this matches with `sv2circos.py` for consistency
svtype2color = dict({'BND': 'black', 'DEL': 'yellow', 'DUP': 'blue', 'INV': 'orange', 'INS': 'green', 'TRA': 'black'})

class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __repr__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)

    def __str__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)

    
    @classmethod
    def fromstring(cls, position_string):
        chromosome = position_string.split(':')[0]
        start = int(position_string.split(':')[1].split('-')[0])
        end = int(position_string.split(':')[1].split('-')[1])
        return Position(chromosome, start, end)

    @staticmethod
    def overlap(position1, position2):
        '''true if position1 and position2 has more than 1 overlapping base'''
        try:
            if isinstance(position1, Position) and isinstance(position2, Position):
                if position1.chromosome == position2.chromosome:
                    if min(position1.end, position2.end) > max (position1.start, position2.start):
                        return True
                    else:
                        return False
                else:
                    return False  # cannot compare if two positions are in different chromosome
            else:
                return None # has to be Posiiton class.
        except:
            Exception
    
    def extend(self, direction, basepairs):
        """extends objects in by specified base pairs, either upstream, downstream, or both"""
        if direction=="up":
            return Position(self.chromosome, max(0, self.start-basepairs), end)
        elif direction=="down":
            return Position(self.chromosome, self.start, self.end + basepairs)
        elif direction=="both":
            return Position(self.chromosome, max(0, self.start - basepairs), self.end + basepairs)
        else:
            print('direction has to be either up, down, or both')
            raise ValueError


def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='separate mutect only calls, strelka only calls, and intersect calls from two vcf files')

    parser.add_argument('-o', '--outputDIR', required=True, default=os.getcwd())
    parser.add_argument('--mutect', required=True, help='Mutect VCF')
    parser.add_argument('--strelka', required=True, help='Strelka VCF')
    parser.add_argument('-s', '--sampleName', required=True, help='Sample Name to be used for output .circosData file')
    args = vars(parser.parse_args())

    outputDIR = args['outputDIR']
    sampleName = args['sampleName']
    mutect = args['mutect']
    strelka = args['strelka']
    return outputDIR, sampleName, mutect, strelka


def vcf2SNVPosition(vcf_file):
    ''' create a generator of Position object of SNVs from a given vcf file
    '''
    variant_list = set()

    for variant in cyvcf2.VCF(vcf_file):
        if variant.FILTER == None:
            variantString= f'{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}'
            variant_list.add(variantString) 
    return variant_list


    

def compare_snvs(vcf1, vcf2, sampleName, outputDIR):
    vcf1_variants = vcf2SNVPosition(vcf1)
    vcf2_variants = vcf2SNVPosition(vcf2)
    in1only = outputDIR + '/' + sampleName + '_mutect_only.vcf'
    in2only = outputDIR + '/' + sampleName + '_strelka_only.vcf'
    inboth = outputDIR + '/' + sampleName + '_mns.vcf'

    
    os.system('bcftools view -h ' + vcf1 + ' > ' + in1only)
    os.system('bcftools view -h ' + vcf2 + ' > ' + in2only)
    os.system('bcftools view -h ' + vcf1 + ' > ' + inboth)

    with open(in1only, 'a') as h:
        for variant in cyvcf2.VCF(vcf1):
            variantString = f'{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}'
            if not variantString in vcf2_variants:
                h.write(str(variant))

    with open(in2only, 'a') as h:
        for variant in cyvcf2.VCF(vcf2):
            variantString = f'{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}'
            if not variantString in vcf1_variants:
                h.write(str(variant))

    with open(inboth, 'a') as h:
        for variant in cyvcf2.VCF(vcf1):
            variantString = f'{variant.CHROM}:{variant.POS}_{variant.REF}>{variant.ALT[0]}'
            if variantString in vcf2_variants:
                h.write(str(variant))


    

    return 0

def main():
    outputDIR, sampleName, mutect, strelka = argument_parser()
    print(outputDIR, sampleName, mutect, strelka)
    compare_snvs(mutect, strelka, sampleName, outputDIR)
    return 0

if __name__=='__main__':
    main()

