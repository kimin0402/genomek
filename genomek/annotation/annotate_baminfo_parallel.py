#!/home/users/changhyunnam/anaconda3/envs/snakemake/bin/python
#211122 Modify to include rasm

import os
import subprocess
import shlex
import re
import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam
from collections import Counter
import statistics
import math
import pyranges as pr

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='tsv file to annotate bam information')
    parser.add_argument('-o', '--output', required=True, help='output tsv file')
    parser.add_argument('-t', '--tumor_bam', required=True, help='tumor bam file')
    parser.add_argument('-n', '--normal_bam', required=True, help='normal bam file')
    parser.add_argument('-r', '--reference', required=True, help='Reference fasta file')
    parser.add_argument('-c', '--chromosome', required=False, help='target chromosome')
    parser.add_argument('-p', '--multicore', required=False, type=int, default=1, help='Number of multiprocessor to parallelize computation')
    args = vars(parser.parse_args())
    return args['input'], args['output'], args['tumor_bam'], args['normal_bam'], args['reference'], args['chromosome'], args['multicore']

def annotate_polyns(chrom, pos, reference, fai_dict, max_repeat_unit_size=5, window_size=5):
    chrom, pos = str(chrom), int(pos)
    ref = pysam.FastaFile(reference)

    if chrom in fai_dict.keys(): chrom_size = fai_dict[chrom]
    else: sys.exit(f"{chrom}:{pos} is not documented in the fai index file")

    up_base_list=[]
    up_n_list=[]
    down_base_list=[]
    down_n_list=[]
    
    for s in range(1,max_repeat_unit_size+1,1):
        up_base = ref.fetch(chrom,pos-s-1,pos-1)
        up_base_list.append(up_base)
        up_n=1
        while pos-up_n*s-s-1 >= 0:
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
        while pos+down_n*s+s <= chrom_size:
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

def comb(n,r):
    """nCR"""
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def fisher_exact(RDF, RDR, ADF, ADR):
    total = RDF + RDR + ADF + ADR
    ref_reads = RDF + RDR
    alt_reads = ADF + ADR
    return comb(ref_reads, RDF) * comb(alt_reads, ADF)/comb(total, RDF + ADF)

def strand_bias(RDF, RDR, ADF, ADR):
    ref_reads = RDF + RDR
    alt_reads = ADF + ADR
    
    pvalues = []
    for i in range(0, min(RDF + RDR, RDF + ADF)+1):
        rdf = i
        adf = RDF + ADF - rdf
        rdr = RDF + RDR - rdf
        adr = ref_reads + alt_reads - rdf - adf - rdr
        if adf >= 0 and rdr >= 0 and adr >= 0:
            pvalue = fisher_exact(rdf, rdr, adf, adr)
            pvalues.append(pvalue)
            if rdf==RDF and adf == ADF and rdr == RDR:
                the_pvalue = pvalue
    fs_value = 0
    for p in pvalues:
        if p <= the_pvalue:
            fs_value += p

    return fs_value


def bgzip_tabix(vcf_file):
    """compress and index vcf file"""
    
    bgzip = "/home/users/changhyunnam/anaconda3/envs/snakemake/bin/bgzip"
    tabix = "/home/users/changhyunnam/anaconda3/envs/snakemake/bin/tabix"

    indexed_vcf = vcf_file + '.gz'
    cmd = f'{bgzip} -f {vcf_file}'
    compress = subprocess.Popen(shlex.split(cmd))
    compress.wait()

    cmd = f'{tabix} -f -p vcf {indexed_vcf}'
    index_run = subprocess.Popen(shlex.split(cmd))
    index_run.wait()
    
    return indexed_vcf

def rerun_mutect_by_pos(chrom,pos,ref,alt,input_tsv,normal_bamFile,tumor_bamFile,reference):
    interval_window = 250

    sampleID = os.path.basename(input_tsv).split(".")[0]
    mutect_dir = os.path.join(os.path.dirname(input_tsv),"mutect",sampleID)
    os.system(f"mkdir -p {mutect_dir}")

    mutect_bamFile = os.path.join(mutect_dir,f"{sampleID}.{chrom}.{pos}.mutect.bam")
    mutect_vcf     = re.sub(".bam$",".vcf",mutect_bamFile)
    mutect_outvcf  = re.sub(".vcf$",".out.vcf.gz",mutect_vcf)
    mutect_log     = re.sub(".vcf$",".log",mutect_vcf)
    
    vcf_header = "/home/users/changhyunnam/scripts/snv_ft/baminfo/vcf_header.vcf"
    with open(mutect_vcf,"w") as f:
        with open(vcf_header,"r") as g: f.write(g.read())
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
                                            
    mutect_vcfgz = bgzip_tabix(mutect_vcf)
    
    JAVA = "/usr/local/bin/java"
    GATK4 = "/home/users/changhyunnam/tools/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar"
    normal_bam = pysam.AlignmentFile(normal_bamFile)
    normal_ID = normal_bam.header['RG'][0]['SM']

    interval = f"{chrom}:{pos-interval_window}-{pos+interval_window}"

    cmd = f"{JAVA} -Xmx32G -jar {GATK4} Mutect2 -I {normal_bamFile} -I {tumor_bamFile} -normal {normal_ID} -L {interval} -alleles {mutect_vcfgz} -R {reference} -O {mutect_outvcf} -bamout {mutect_bamFile}"
    # print(cmd)
    
    log_file = open(mutect_log,"w")
    run_mutect = subprocess.Popen(shlex.split(cmd),stderr=log_file)
    run_mutect.wait()
    log_file.close()

    return mutect_bamFile

def annotate_readinfo_mutect(chrom,pos,ref,alt,input_tsv,normal_bamFile,tumor_bamFile,reference):
    chrom, pos = str(chrom), int(pos)
    ref, alt = str(ref), str(alt)
    # print(f"{chrom}:{pos}:{ref}:{alt}")

    mutect_bamFile = rerun_mutect_by_pos(chrom,pos,ref,alt,input_tsv,normal_bamFile,tumor_bamFile,reference)
    bam = pysam.AlignmentFile(mutect_bamFile)

    ref_reads           = 0
    alt_reads           = 0
    other_reads         = 0

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality<20: continue
        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        try:
            element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == pos-1][0]
            pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0] 
            ridx    = element[0]
        # read_index(ridx) and pair_index(pidx) are not always the same
        # when reads possess deletion, pair_index increases by the size of deletion while read_index is constant
        except:
            continue
        if read.query_qualities[ridx]<20: continue

        b = read.query_sequence[ridx][0]

        if len(ref) == len(alt): # SNV
            if b == ref: ref_reads+=1
            elif b == alt: alt_reads+=1
            else: other_reads+=1
        
        elif len(ref) < len(alt): # Insertion
            if b != ref: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][1] is not None: 
                ref_reads+=1
            elif (pidx+len(alt)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[1] for x in aligned_pairs[(pidx+1):(pidx+len(alt))]] == [None] * (len(alt)-1)) and (aligned_pairs[pidx+len(alt)][1] is not None):
                alt_reads+=1
            else:
                other_reads+=1

        elif len(ref) > len(alt): # Deletion
            if b != ref[0]: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][0] is not None: 
                ref_reads+=1
            elif (pidx+len(ref)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[0] for x in aligned_pairs[(pidx+1):(pidx+len(ref))]] == [None] * (len(ref)-1)) and (aligned_pairs[pidx+len(ref)][0] is not None):
                alt_reads+=1
            else:
                other_reads+=1

    total_depth = ref_reads+alt_reads+other_reads
    if total_depth: vaf = alt_reads/total_depth
    else: vaf = np.nan

    return f'{total_depth}:{ref_reads}:{alt_reads}:{other_reads}:{vaf:0.5f}'

def annotate_readinfo_count(chrom,pos,ref,alt,bamFile):
    chrom, pos = str(chrom), int(pos)
    ref, alt = str(ref), str(alt)
    bam = pysam.AlignmentFile(bamFile)

    print(f"{chrom}:{pos}:{ref}:{alt}")

    ref_reads           = 0
    alt_reads           = 0
    other_reads         = 0

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality<20: continue
        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        try:
            element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == pos-1][0]
            pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0] 
            ridx    = element[0]
        # read_index(ridx) and pair_index(pidx) are not always the same
        # when reads possess deletion, pair_index increases by the size of deletion while read_index is constant
        except:
            continue
        if read.query_qualities[ridx]<20: continue

        b = read.query_sequence[ridx][0]

        if len(ref) == len(alt): # SNV
            if b == ref: ref_reads+=1
            elif b == alt: alt_reads+=1
            else: other_reads+=1
        
        elif len(ref) < len(alt): # Insertion
            if b != ref: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][1] is not None: 
                ref_reads+=1
            elif (pidx+len(alt)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[1] for x in aligned_pairs[(pidx+1):(pidx+len(alt))]] == [None] * (len(alt)-1)) and (aligned_pairs[pidx+len(alt)][1] is not None):
                alt_reads+=1
            else:
                other_reads+=1

        elif len(ref) > len(alt): # Deletion
            if b != ref[0]: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][0] is not None: 
                ref_reads+=1
            elif (pidx+len(ref)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[0] for x in aligned_pairs[(pidx+1):(pidx+len(ref))]] == [None] * (len(ref)-1)) and (aligned_pairs[pidx+len(ref)][0] is not None):
                alt_reads+=1
            else:
                other_reads+=1

    total_depth = ref_reads+alt_reads+other_reads
    if total_depth: vaf = alt_reads/total_depth
    else: vaf = np.nan

    return f'{total_depth}:{ref_reads}:{alt_reads}:{other_reads}:{vaf:0.5f}'

def annotate_readinfo_normal(chrom,pos,ref,alt,bamFile,up_polyn,down_polyn):
    chrom, pos = str(chrom), int(pos)
    ref, alt = str(ref), str(alt)
    bam = pysam.AlignmentFile(bamFile)
    up_polyn, down_polyn = int(up_polyn),int(down_polyn)

    print(f"{chrom}:{pos}:{ref}:{alt}")

    ref_reads           = 0
    alt_reads           = 0
    alt_reads_lowqual   = 0
    total_reads_clipped = 0
    total_reads_discord = 0
    other_reads         = 0
    total_reads_idctx   = 0
    total_reads_idctx_ext = 0 
    total_NM            = []
    isw = 10 # isw(indel search windows)

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality<20: continue
        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        try:
            element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == pos-1][0]
            pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0] 
            ridx    = element[0]
        # read_index(ridx) and pair_index(pidx) are not always the same
        # when reads possess deletion, pair_index increases by the size of deletion while read_index is constant
        except:
            continue

        b = read.query_sequence[ridx][0]

        if read.query_qualities[ridx]<20 and len(ref) == len(alt): 
            if b == alt: alt_reads_lowqual+=1
            continue

        indel_context = 1 if len([aligned_pairs[i] for i in range(max(0,pidx-isw), min(pidx+isw+1,len(aligned_pairs)-1),1) if i != pidx and (aligned_pairs[i][0] is None or aligned_pairs[i][1] is None)]) > 0 else 0 
        indel_context_ext = 1 if len([aligned_pairs[i] for i in range(max(0,pidx-isw-(up_polyn-1)), min(pidx+isw+1+(down_polyn-1),len(aligned_pairs)-1),1) if i != pidx and (aligned_pairs[i][0] is None or aligned_pairs[i][1] is None)]) > 0 else 0 

        if len(ref) == len(alt): # SNV
            if b == ref: ref_reads+=1
            elif b == alt: alt_reads+=1
            else: other_reads+=1
        
        elif len(ref) < len(alt): # Insertion
            if b != ref: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][1] is not None: 
                ref_reads+=1
            elif (pidx+len(alt)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[1] for x in aligned_pairs[(pidx+1):(pidx+len(alt))]] == [None] * (len(alt)-1)) and (aligned_pairs[pidx+len(alt)][1] is not None):
                if read.query_sequence[ridx:ridx+len(alt)] == alt: alt_reads+=1
                else: other_reads+=1
            else:
                other_reads+=1

        elif len(ref) > len(alt): # Deletion
            if b != ref[0]: continue
            if (pidx+1) == len(aligned_pairs): continue
            if aligned_pairs[pidx+1][0] is not None: 
                ref_reads+=1
            elif (pidx+len(ref)) >= len(aligned_pairs):
                other_reads+=1
            elif ([x[0] for x in aligned_pairs[(pidx+1):(pidx+len(ref))]] == [None] * (len(ref)-1)) and (aligned_pairs[pidx+len(ref)][0] is not None):
                alt_reads+=1
            else:
                other_reads+=1

        total_reads_idctx+=indel_context
        total_reads_idctx_ext+=indel_context_ext 
        if re.search(r'S|H', read.cigarstring): total_reads_clipped+=1
        if read.is_paired:
            if not read.mate_is_unmapped:
                if chrom != read.next_reference_name or abs(read.template_length) > 1000: total_reads_discord+=1 
        total_NM.append(read.get_tag("NM"))

    total_depth = ref_reads+alt_reads+other_reads
    if total_depth: vaf = alt_reads/total_depth
    else: vaf = np.nan

    if total_NM:
        mean_NM_total   = sum(total_NM)/len(total_NM)
        median_NM_total = statistics.median(total_NM)
    else: 
        mean_NM_total   = np.nan
        median_NM_total = np.nan

    return f'{total_depth}:{ref_reads}:{alt_reads}:{total_reads_clipped}:{total_reads_discord}:{alt_reads_lowqual}:{other_reads}:{vaf:0.5f}:{total_reads_idctx}:{total_reads_idctx_ext}:{mean_NM_total:0.2f}:{median_NM_total}'


def annotate_readinfo(chrom,pos,ref,alt,bamFile,germline_bamFile,up_polyn,down_polyn):
    chrom, pos = str(chrom), int(pos)
    ref, alt = str(ref), str(alt)
    bam = pysam.AlignmentFile(bamFile)
    germline_bam = pysam.AlignmentFile(germline_bamFile)
    up_polyn, down_polyn = int(up_polyn),int(down_polyn)

    print(f"{chrom}:{pos}:{ref}:{alt}")

    ref_reads         = 0
    ref_reads_clipped = 0
    ref_reads_discord = 0 
    ref_reads_fwd     = 0
    ref_reads_rev     = 0 
    alt_reads         = 0
    alt_reads_clipped = 0
    alt_reads_discord = 0 
    alt_reads_fwd     = 0
    alt_reads_rev     = 0
    other_reads       = 0
    ref_reads_idctx   = 0 # idctx(indel context)
    alt_reads_idctx   = 0
    ref_reads_idctx_ext = 0 # idctx(indel context)
    alt_reads_idctx_ext = 0
    ref_NM            = []
    alt_NM            = []
    alt_RP            = []
    ref_MQ            = []
    alt_MQ            = []
    ref_up_repeat     = []
    ref_down_repeat   = []
    alt_rasm_seq_list = []
    TD_poor = 0
    TD_MQ0  = 0 
    TD_DEL  = 0
    TD_rasm = 0 
    AD_rasm = 0 
    isw = 10 # isw(indel search windows)
    rasm_window = 15 

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality<20: 
            TD_poor+=1
            if read.mapping_quality == 0: TD_MQ0 += 1
            continue
        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        try:
            element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == pos-1][0]
            pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0] 
            ridx    = element[0]
        # read_index(ridx) and pair_index(pidx) are not always the same
        # when reads possess deletion, pair_index increases by the size of deletion while read_index is constant
        except:
            TD_DEL+=1 
            continue
        if read.query_qualities[ridx]<20 and len(ref) == len(alt): continue

        b = read.query_sequence[ridx][0]
        indel_context = 1 if len([aligned_pairs[i] for i in range(max(0,pidx-isw), min(pidx+isw+1,len(aligned_pairs)-1),1) if i != pidx and (aligned_pairs[i][0] is None or aligned_pairs[i][1] is None)]) > 0 else 0 
        indel_context_ext = 1 if len([aligned_pairs[i] for i in range(max(0,pidx-isw-(up_polyn-1)), min(pidx+isw+1+(down_polyn-1),len(aligned_pairs)-1),1) if i != pidx and (aligned_pairs[i][0] is None or aligned_pairs[i][1] is None)]) > 0 else 0 

        if len(ref) == len(alt): # SNV
            if b == ref: vartype = "REF"
            elif b == alt: vartype = "ALT"
            else: vartype = "OTHER"
    
        elif len(ref) < len(alt): # Insertion
            if b != ref: continue
            elif (pidx+1) == len(aligned_pairs): continue
            elif aligned_pairs[pidx+1][1] is not None: vartype = "REF"
            elif (pidx+len(alt)) >= len(aligned_pairs): vartype = "OTHER"
            elif ([x[1] for x in aligned_pairs[(pidx+1):(pidx+len(alt))]] == [None] * (len(alt)-1)) and (aligned_pairs[pidx+len(alt)][1] is not None):
                if read.query_sequence[ridx:ridx+len(alt)] == alt: vartype = "ALT"
                else: vartype = "OTHER"
            else: vartype = "OTHER"
    
        elif len(ref) > len(alt): # Deletion
            if b != ref[0]: continue
            elif (pidx+1) == len(aligned_pairs): continue
            elif aligned_pairs[pidx+1][0] is not None: vartype = "REF"
            elif (pidx+len(ref)) >= len(aligned_pairs): vartype = "OTHER"
            elif ([x[0] for x in aligned_pairs[(pidx+1):(pidx+len(ref))]] == [None] * (len(ref)-1)) and (aligned_pairs[pidx+len(ref)][0] is not None):
                vartype = "ALT"
            else: vartype = "OTHER"

        if vartype == "REF":
            ref_reads+=1
            ref_MQ.append(read.mapping_quality)
            ref_NM.append(read.get_tag("NM"))
            ref_reads_idctx+=indel_context
            ref_reads_idctx_ext+=indel_context_ext

            if re.search(r'S|H', read.cigarstring): ref_reads_clipped+=1
            if read.is_paired:
                if not read.mate_is_unmapped:
                    if chrom != read.next_reference_name or abs(read.template_length) > 1000: ref_reads_discord+=1 

            if read.is_reverse: ref_reads_rev+=1
            else: ref_reads_fwd+=1

            ref_down_repeat.append(len(re.findall(b+"+",read.query_sequence[ridx:min(ridx+10,read.query_length)])[0]))
            ref_up_list = re.findall(alt+"+",read.query_sequence[max(0,ridx-10):ridx])
            ref_up_element = len(ref_up_list[-1]) if len(ref_up_list) else 0
            ref_up_repeat.append(ref_up_element)

        elif vartype == "ALT":
            alt_reads+=1
            alt_MQ.append(read.mapping_quality)
            alt_NM.append(read.get_tag("NM"))
            alt_RP.append(aligned_pairs[pidx][0]/read.query_length)
            alt_reads_idctx+=indel_context
            alt_reads_idctx_ext+=indel_context_ext

            if re.search(r'S|H', read.cigarstring): alt_reads_clipped+=1
            if read.is_paired:
                if not read.mate_is_unmapped:
                    if chrom != read.next_reference_name or abs(read.template_length) > 1000: alt_reads_discord+=1 

            if read.is_reverse: alt_reads_rev+=1
            else: alt_reads_fwd+=1

            ridx_w1_list = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] == pos-1-rasm_window]
            ridx_w1 = ridx_w1_list[0][0] if len(ridx_w1_list) else 0
            
            ridx_w2_list = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] == pos-1+rasm_window+len(ref)+len(alt)-2]
            # To extend sequence length for indel
            ridx_w2 = ridx_w2_list[0][0] if len(ridx_w2_list) else read.query_length-1
            alt_rasm_seq_list.append(read.query_sequence[ridx_w1:ridx_w2])

        else: # vartype == "OTHER"
            other_reads+=1

    if ref_NM:
        mean_NM_ref   = sum(ref_NM)/len(ref_NM)
        median_NM_ref = statistics.median(ref_NM)
    else: 
        mean_NM_ref   = np.nan
        median_NM_ref = np.nan
    if alt_NM: 
        mean_NM_alt   = sum(alt_NM)/len(alt_NM)
        median_NM_alt = statistics.median(alt_NM)
    else: 
        mean_NM_alt   = np.nan
        median_NM_alt = np.nan
    if alt_RP:  
        mean_RP_alt = sum(alt_RP)/len(alt_RP)
        if len(alt_RP)>1: std_RP_alt = statistics.stdev(alt_RP)
        else: std_RP_alt = 0
    else: 
        mean_RP_alt = np.nan
        std_RP_alt = np.nan
    if ref_MQ: mean_MQ_ref = sum(ref_MQ)/len(ref_MQ)
    else: mean_MQ_ref = np.nan
    if alt_MQ: mean_MQ_alt = sum(alt_MQ)/len(alt_MQ)
    else: mean_MQ_alt = np.nan

    total_depth = ref_reads+alt_reads+other_reads
    if total_depth:
        vaf = alt_reads/total_depth
        sb = strand_bias(ref_reads_fwd, ref_reads_rev, alt_reads_fwd, alt_reads_rev)
    else:
        vaf = np.nan
        sb = np.nan

    if ref_up_repeat: 
        try:
            ref_up_polyn = max([i for i, j in Counter(ref_up_repeat).most_common(2) if j>ref_reads*0.4])
        except:
            ref_up_polyn = Counter(ref_up_repeat).most_common(1)[0][0]
    else: ref_up_polyn = np.nan

    if ref_down_repeat: 
        try:
            ref_down_polyn = max([i for i, j in Counter(ref_down_repeat).most_common(2) if j>ref_reads*0.4])
        except:
            ref_down_polyn = Counter(ref_down_repeat).most_common(1)[0][0]
    else: ref_down_polyn = np.nan

    alt_rasm_seq = "N"
    alt_rasm_seq_N = 0
    for k, v in Counter(alt_rasm_seq_list).items():
        if v > alt_rasm_seq_N: alt_rasm_seq, alt_rasm_seq_N = k, v
        elif v == alt_rasm_seq_N and len(k) > len(alt_rasm_seq): alt_rasm_seq, alt_rasm_seq_N = k, v

    for read in germline_bam.fetch(chrom,pos-1,pos):
        # if read.is_duplicate: continue
        # if read.mapping_quality<20: continue    
        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        ridx_w1_list = [pair for pair in aligned_pairs if pair[1] == pos-1-rasm_window]
        if not ridx_w1_list: continue
        if ridx_w1_list[0][0] is None: continue
        rasm_seq = read.query_sequence[ridx_w1_list[0][0]:(ridx_w1_list[0][0]+len(alt_rasm_seq))]        

        TD_rasm += 1 
        if rasm_seq == alt_rasm_seq: AD_rasm += 1
    
    if TD_rasm: VAF_rasm = AD_rasm/TD_rasm
    else: VAF_rasm = np.nan

    return f'{total_depth}:{ref_reads}:{alt_reads}:{ref_reads_clipped}:{ref_reads_discord}:{alt_reads_clipped}:{alt_reads_discord}:{other_reads}:{vaf:0.4f}:{ref_reads_fwd}:{ref_reads_rev}:{alt_reads_fwd}:{alt_reads_rev}:{sb:0.4f}:{mean_NM_ref:0.2f}:{mean_NM_alt:0.2f}:{median_NM_ref}:{median_NM_alt}:{mean_RP_alt:0.2f}:{std_RP_alt:0.2f}:{mean_MQ_ref:0.2f}:{mean_MQ_alt:0.2f}:{ref_reads_idctx}:{ref_reads_idctx_ext}:{alt_reads_idctx}:{alt_reads_idctx_ext}:{ref_up_polyn}:{ref_down_polyn}:{TD_poor}:{TD_MQ0}:{TD_DEL}:{TD_rasm}:{AD_rasm}:{VAF_rasm:0.4f}'

def find_mate(bamFile,read):
    """finds the mate of read"""
    mate_bam = pysam.AlignmentFile(bamFile)
    for mate_read in mate_bam.fetch(read.next_reference_name, read.next_reference_start, read.next_reference_start + 1):
        if mate_read.query_name == read.query_name:return mate_read
    return None

def find_read_range(bamFile, chrom, pos):
    bam = pysam.AlignmentFile(bamFile)
    read_range1, read_range2 = pos-1, pos-1
    mate_read_range1, mate_read_range2 = pos-1, pos-1

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality<20: continue
        if not read.is_paired: continue

        if read.reference_start<read_range1: read_range1 = read.reference_start
        if read.reference_end>read_range2: read_range2 = read.reference_end

        mate_read = find_mate(bamFile,read)
        if mate_read.reference_name != chrom: continue
        if mate_read.reference_start<mate_read_range1: mate_read_range1 = mate_read.reference_start
        if mate_read.reference_end>mate_read_range2: mate_read_range2 = mate_read.reference_end
        
    return read_range1, read_range2, mate_read_range1, mate_read_range2

def find_nearby_germline(bamFile, reference, chrom, pos, mate_read_range1, mate_read_range2, MQ_cutoff = 20, BQ_cutoff = 20, depth_cutoff = 10):
    bam = pysam.AlignmentFile(bamFile)
    ref = pysam.FastaFile(reference)

    germline_list = []
    indel_pat = re.compile('([ACGTN])([\+|-])([0-9]+)([ACGTN]+)')
    # Insertion : 'C+4TCTG'
    # Deletion : 'T-39NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

    for pos_i in range(mate_read_range1, mate_read_range2+1):
        if pos_i == pos-1: continue

        pup = bam.pileup(chrom,pos_i,pos_i+1,min_mapping_quality=MQ_cutoff,min_base_quality=BQ_cutoff,truncate=True)
        pupc = next(pup)

        base_list = [base.upper() for base in pupc.get_query_sequences(add_indels=True)]
        total_depth = len(base_list)
        if total_depth < depth_cutoff: continue # if depth < 10, go to next position ()

        ref_base = ref.fetch(chrom,pos_i,pos_i+1)
        if set(ref_base) == set(base_list): continue # if all bases are ref base, go to next position

        base_count = Counter(base_list)
        if base_count[ref_base]/total_depth > 0.75: continue 
        # if vaf of ref base > 0.75 ( = vaf of total alt base < 0.25), go to next position

        alt_base_list = [base for base in base_count.keys() if base != ref_base and base != "*"] # base "*" for deleted sites
        alt_vaf_list = [base_count[base]/total_depth for base in alt_base_list]

        for alt_base, alt_vaf in zip(alt_base_list,alt_vaf_list):
            if alt_vaf > 0.25:
                indel_re = indel_pat.match(alt_base)
                if indel_re is None: 
                    germline_list.append(f"{chrom}_{pos_i+1}_{ref_base}>{alt_base}")
                elif indel_re.group(2) == "+":
                    germline_list.append(f"{chrom}_{pos_i+1}_{indel_re.group(1)}>{indel_re.group(1)+indel_re.group(4)}")
                elif indel_re.group(2) == "-":
                    del_bases = ref.fetch(chrom,pos_i+1,pos_i+int(indel_re.group(3))+1) # Replace N (deleted base) to reference sequence
                    germline_list.append(f"{chrom}_{pos_i+1}_{indel_re.group(1)+del_bases}>{indel_re.group(1)}")
                    
    return(germline_list) # pos is 1-based


def find_vartype(read,ridx,ref,alt):
    b = read.query_sequence[ridx][0]

    if len(ref) == len(alt): # SNV
        if b == ref: vartype = "REF"
        elif b == alt: vartype = "ALT"
        else: vartype = "OTHER"

    elif len(ref) < len(alt): # Insertion
        if b != ref: vartype = "OTHER"
        elif (pidx+1) == len(aligned_pairs): vartype = "OTHER"
        elif aligned_pairs[pidx+1][1] is not None: vartype = "REF"
        elif (pidx+len(alt)) >= len(aligned_pairs): vartype = "OTHER"
        elif ([x[1] for x in aligned_pairs[(pidx+1):(pidx+len(alt))]] == [None] * (len(alt)-1)) and (aligned_pairs[pidx+len(alt)][1] is not None):
            vartype = "ALT"
        else: vartype = "OTHER"

    elif len(ref) > len(alt): # Deletion
        if b != ref[0]: vartype = "OTHER"
        elif (pidx+1) == len(aligned_pairs): vartype = "OTHER"
        elif aligned_pairs[pidx+1][0] is not None: vartype = "REF"
        elif (pidx+len(ref)) >= len(aligned_pairs): vartype = "OTHER"
        elif ([x[0] for x in aligned_pairs[(pidx+1):(pidx+len(ref))]] == [None] * (len(ref)-1)) and (aligned_pairs[pidx+len(ref)][0] is not None):
            vartype = "ALT"
        else: vartype = "OTHER"
    
    return vartype


def phase_germline_correct_NM(bamFile,chrom,pos,ref,alt,germline_list,NM_list,MQ_cutoff=20,BQ_cutoff=20):    
    if len(germline_list) == 0: return "None", "None", NM_list
    
    bam = pysam.AlignmentFile(bamFile)
    phase_dict = {germline_query:{"WT":0,"VAR":0} for germline_query in germline_list}
    NM_dict = {"REF":[],"ALT":[]}

    for read in bam.fetch(chrom,pos-1,pos):
        if read.is_duplicate: continue
        if read.mapping_quality < MQ_cutoff: continue
        if not read.is_paired: continue        

        aligned_pairs = read.get_aligned_pairs() # list of (read index(0-based), position(0-based))
        try:
            element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == pos-1][0]
            pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == element][0] 
            ridx    = element[0]
        # read_index(ridx) and pair_index(pidx) are not always the same
        # when reads possess deletion, pair_index increases by the size of deletion while read_index is constant
        except:
            continue
        if read.query_qualities[ridx] < BQ_cutoff: continue

        vartype = find_vartype(read,ridx,ref,alt)
        if vartype == "OTHER": continue

        read_range = list(range(read.reference_start, read.reference_end))

        mate_read = find_mate(bamFile, read)
        if mate_read.reference_name == chrom:
            mate_range = list(range(mate_read.reference_start, mate_read.reference_end))
            mate_aligned_pairs = mate_read.get_aligned_pairs()
        else:
            mate_range, mate_aligned_pairs = [], []

        NM_gap = 0
        for germline_query in germline_list:
            germ_pos, germ_ref, germ_alt = re.split(":|_|>",germline_query)[1:]
            germ_pos = int(germ_pos)-1 # 1-based string to 0-based integer

            if germ_pos not in read_range + mate_range: continue
            if germ_pos in read_range:
                try: 
                    germ_element = [pair for pair in aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == germ_pos][0]
                    germ_pidx    = [i for i in range(len(aligned_pairs)) if aligned_pairs[i] == germ_element][0] 
                    germ_ridx    = germ_element[0]
                    germ_vartype = find_read_vartype(read,germ_ridx,germ_ref,germ_alt)
                except:
                    continue

                if vartype == "ALT":
                    if germ_vartype == "REF": phase_dict[germline_query]["WT"]+=1
                    elif germ_vartype == "ALT": phase_dict[germline_query]["VAR"]+=1

                if germ_vartype == "ALT":
                    germ_NM = 1 if len(germ_ref) == 1 and len(germ_alt) == 1 else abs(len(germ_ref)-len(germ_alt))
                    NM_gap+=germ_NM

            if germ_pos in mate_range:
                try: 
                    mate_germ_element = [pair for pair in mate_aligned_pairs if pair[0] is not None and pair[1] is not None and pair[1] == germ_pos][0]
                    mate_germ_pidx    = [i for i in range(len(mate_aligned_pairs)) if mate_aligned_pairs[i] == mate_germ_element][0] 
                    mate_germ_ridx    = mate_germ_element[0]
                    germ_vartype = find_read_vartype(mate_read,mate_germ_ridx,germ_ref,germ_alt)
                except:
                    continue

                if vartype == "ALT":
                    if germ_vartype == "REF": phase_dict[germline_query]["WT"]+=1
                    elif germ_vartype == "ALT": phase_dict[germline_query]["VAR"]+=1

        corrected_NM = read.get_tag("NM") - NM_gap
        NM_dict[vartype].append(corrected_NM)

    for k, v in phase_dict.items():
        if v["WT"] + v["VAR"] == 0: v["phase"] = "unknown"
        elif v["WT"] * v["VAR"] == 0 : v["phase"] = "phased"
        else: v["phase"] = "unphased"

    phase_string = ";".join([f"{k}_WT{v['WT']}_VAR{v['VAR']}_{v['phase']}" for k, v in phase_dict.items()])
    phase_list = [v["phase"] for k,v in phase_dict.items()]
    if "unphased" in phase_list: phase_status = "unphased"
    elif "phased" in phase_list: phase_status = "phased"
    else: phase_status = "unknown"

    ref_NM, alt_NM = NM_dict["REF"], NM_dict["ALT"]
    if ref_NM: mean_NM_ref, median_NM_ref = sum(ref_NM)/len(ref_NM), statistics.median(ref_NM)
    else: mean_NM_ref, median_NM_ref = np.nan, np.nan
    if alt_NM: mean_NM_alt, median_NM_alt = sum(alt_NM)/len(alt_NM), statistics.median(alt_NM)
    else: mean_NM_alt, median_NM_alt = np.nan, np.nan
                             
    return phase_string, phase_status, [f"{mean_NM_ref:0.2f}",f"{mean_NM_alt:0.2f}", median_NM_ref, median_NM_alt]

def annotate_readinfo_nearby(chrom,pos,ref,alt,mean_NM_ref,mean_NM_alt,median_NM_ref,median_NM_alt,bamFile,reference):
    chrom, pos = str(chrom), int(pos)
    ref, alt = str(ref), str(alt)

    #1. Find read range and mate read range
    read_range1, read_range2, mate_read_range1, mate_read_range2 = find_read_range(bamFile,chrom,pos)
    # print(read_range1, read_range2, mate_read_range1, mate_read_range2)

    #2. Find nearby germline variants in mate read range
    germline_list = find_nearby_germline(bamFile,reference,chrom,pos,mate_read_range1,mate_read_range2,20,20,10)
    # print(germline_list)

    #3. Phase nearby germline variants and correct NM values
    NM_list = [mean_NM_ref,mean_NM_alt,median_NM_ref,median_NM_alt]
    phase_string, phase_status, correct_NM_list = phase_germline_correct_NM(bamFile,chrom,pos,ref,alt,germline_list,NM_list,20,20)
    mean_NM_ref_c, median_NM_ref_c, mean_NM_alt_c, median_NM_alt_c = correct_NM_list
    # print(phase_string)
    
    return f'{phase_string}:{phase_status}:{mean_NM_ref_c}:{mean_NM_alt_c}:{median_NM_ref_c}:{median_NM_alt_c}'

def filter_satellite(chrom,pos,pr_satellite):
        chrom, pos = str(chrom), int(pos)
        print(chrom, pos)
        pr_variant = pr.PyRanges(chromosomes = chrom, starts = [pos-1], ends = [pos])
        return "NSAT" if pr_satellite.intersect(pr_variant).empty else "SAT"

def annotate_baminfo(df,tumor_bamFile,normal_bamFile,reference,fai_dict):
    #0. remove variant on satellite region
    # df["satellite"] = df.apply(lambda x: filter_satellite(x['CHROM'],x['POS'],pr_satellite), axis=1)
    # df = df[df["satellite"] == "NSAT"]
    # df = df.drop("satellite",axis=1)
    # if len(df) == 0: return df

    #0. annotate repeat context
    df["polyns"] = df.apply(lambda x: annotate_polyns(x['CHROM'],x['POS'], reference, fai_dict, 5, 5), axis=1)
    df[["up_base","up_polyn","down_base","down_polyn","window_base_type","window_base_n"]] = df["polyns"].str.split(":",expand=True,)
    df = df.drop("polyns", axis=1)
    # print("#0. annotate repeat context : Finished")

    #1. annotate normal readinfo
    df["normal_readinfo"] = df.apply(lambda x: annotate_readinfo_normal(x['CHROM'],x['POS'],x['REF'],x['ALT'],normal_bamFile,x['up_polyn'],x['down_polyn']), axis=1)
    df[["TD_normal","RD_normal","AD_normal","CLIPPED_normal","DISC_normal","ADLQ_normal","OTHER_normal","VAF_normal","ID_context_normal","ID_context_extend_normal","mean_NM_normal","median_NM_normal"]] = df["normal_readinfo"].str.split(":",expand=True)
    df = df.drop("normal_readinfo", axis=1)
    # print("#1. annotate normal readinfo : Finished")

    #2. annotate tumor readinfo
    df["tumor_readinfo"] = df.apply(lambda x: annotate_readinfo(x['CHROM'],x['POS'],x['REF'],x['ALT'],tumor_bamFile,normal_bamFile,x['up_polyn'],x['down_polyn']), axis=1)
    df[["TD","RD","AD","CLIPPED_RD","DISC_RD","CLIPPED_AD","DISC_AD","OTHER","VAF","RDF","RDR","ADF","ADR","SB","mean_NM_ref","mean_NM_alt","median_NM_ref","median_NM_alt","mean_RP_alt","std_RP_alt","mean_MQ_ref","mean_MQ_alt","ID_context_ref","ID_context_extend_ref","ID_context_alt","ID_context_extend_alt","ref_up_polyn","ref_down_polyn","TD_poor","TD_MQ0","TD_DEL","TD_rasm","AD_rasm","VAF_rasm"]] = df["tumor_readinfo"].str.split(":",expand=True)
    df = df.drop("tumor_readinfo", axis=1)
    # print("#2. annotate tumor readinfo : Finished")        



    #3. annotate mutect readinfo
    # df["mutect_readinfo"] = df.apply(lambda x: annotate_readinfo_count(x['CHROM'],x['POS'],x['REF'],x['ALT'],mutect_bamFile), axis=1)
    # df[["TD_mutect","RD_mutect","AD_mutect","OTHER_mutect","VAF_mutect"]] = df["mutect_readinfo"].str.split(":",expand=True)
    # df = df.drop("mutect_readinfo", axis=1)
    # print("#3. annotate mutect readinfo : Finished")

    #4. annotate repeat context
    # df["polyns"] = df.apply(lambda x: annotate_polyns(x['CHROM'],x['POS'], reference, fai_dict, 5, 5), axis=1)
    # df[["up_base","up_polyn","down_base","down_polyn","window_base_type","window_base_n"]] = df["polyns"].str.split(":",expand=True,)
    # df = df.drop("polyns", axis=1)
    # print("#4. annotate repeat context : Finished")

    return df

# def annotate_baminfo_normal(df,normal_bamFile):
#     #1. annotate normal readinfo
#     df["normal_readinfo"] = df.apply(lambda x: annotate_readinfo_normal(x['CHROM'],x['POS'],x['REF'],x['ALT'],normal_bamFile), axis=1)
#     df[["TD_normal","RD_normal","AD_normal","CLIPPED_normal","OTHER_normal","VAF_normal","ID_context_normal","mean_NM_normal","median_NM_normal"]] = df["normal_readinfo"].str.split(":",expand=True)
#     df = df.drop("normal_readinfo", axis=1)
#     print("#1. annotate normal readinfo : Finished")
#     return df

# def annotate_baminfo_tumor(df,tumor_bamFile,reference,fai_dict):
#     #2. annotate tumor readinfo
#     df["tumor_readinfo"] = df.apply(lambda x: annotate_readinfo(x['CHROM'],x['POS'],x['REF'],x['ALT'],tumor_bamFile), axis=1)
#     df[["TD","RD","AD","CLIPPED_AD","OTHER","VAF","RDF","RDR","ADF","ADR","SB","mean_NM_ref","mean_NM_alt","median_NM_ref","median_NM_alt","mean_RP_alt","ID_context_ref","ID_context_alt","ref_up_polyn","ref_down_polyn"]] = df["tumor_readinfo"].str.split(":",expand=True)
#     df = df.drop("tumor_readinfo", axis=1)
#     print("#2. annotate tumor readinfo : Finished")        

#     #2b. annotate repeat context
#     df["polyns"] = df.apply(lambda x: annotate_polyns(x['CHROM'],x['POS'], reference, fai_dict, 5, 5), axis=1)
#     df[["up_base","up_polyn","down_base","down_polyn","window_base_type","window_base_n"]] = df["polyns"].str.split(":",expand=True,)
#     df = df.drop("polyns", axis=1)
#     print("#2b. annotate repeat context : Finished")
#     return df

# def annotate_baminfo_mutect(df,input_tsv,normal_bamFile,tumor_bamFile,reference):
#     #3. annotate mutect readinfo
#     df["mutect_readinfo"] = df.apply(lambda x: annotate_readinfo_mutect(x['CHROM'],x['POS'],x['REF'],x['ALT'],input_tsv,normal_bamFile,tumor_bamFile,reference), axis=1)
#     df[["TD_mutect","RD_mutect","AD_mutect","OTHER_mutect","VAF_mutect"]] = df["mutect_readinfo"].str.split(":",expand=True)
#     df = df.drop("mutect_readinfo", axis=1)
#     print("#3. annotate mutect readinfo : Finished")
#     return df

# def annotate_phase(df,tumor_bamFile,reference):
#     df["nearby_readinfo"] = df.apply(lambda x: annotate_readinfo_nearby(x['CHROM'],x['POS'],x['REF'],x['ALT'],x['mean_NM_ref'],x['mean_NM_alt'],x['median_NM_ref'],x['median_NM_alt'],tumor_bamFile,reference), axis=1)
#     df[["phase_string","phase_status","mean_NM_ref_corrected","mean_NM_alt_corrected","median_NM_ref_corrected","median_NM_alt_corrected"]] = df["nearby_readinfo"].str.split(":",expand=True)
#     df = df.drop("nearby_readinfo", axis=1)
#     print("4. annotate nearby readinfo : Finished")
#     return df

def annotate_as_split(df,annotate_fxn,ncore,*args):
    if (len(df) < ncore) or (ncore == 1) :
        df_annotated = annotate_fxn(df,*args)
    else: 
        df_split = np.array_split(df,ncore)
        arg_list = [] 
        for df in df_split:
            arg_list.append((df,)+args)
        with mp.Pool(ncore) as pool:
            result = pool.starmap(annotate_fxn,arg_list)
        df_annotated = pd.concat(result)
    return df_annotated

def main():
    input_tsv, output_tsv, tumor_bamFile, normal_bamFile, reference, target_chrom, ncore = argument_parser()

    if os.path.isfile(reference+".fai"):
        fai = pd.read_csv(reference+".fai",sep="\t",names=["NAME","LENGTH","A","B","C"])
        fai_dict = {chrom:size for chrom, size in zip(fai["NAME"], fai["LENGTH"]) if not re.search("GL",chrom)}
    else: sys.exit(f"fai index file was not found for given reference file")

    df = pd.read_csv(input_tsv, sep='\t',dtype={'CHROM':'str'})
    if target_chrom is not None: df = df[df["CHROM"] == str(target_chrom)]
    if len(df) == 0:
        df.to_csv(output_tsv,sep="\t",index=False)
    else:
        df_bam = annotate_as_split(df,annotate_baminfo,ncore,tumor_bamFile,normal_bamFile,reference,fai_dict)
        df_bam.to_csv(output_tsv,sep="\t",index=False)

if __name__=='__main__':
    main()
