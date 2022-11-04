import sys
import pysam
import subprocess
import tempfile
import pyranges as pr
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from cachetools import cached
from cachetools.keys import hashkey
from ..tools import print_err


def alnlist_to_seqdict(aln_list, seq_idt: str) -> dict:
    result = {}
    for i, aln in enumerate(aln_list):
        result[f"{seq_idt}_{i}"] = aln.query_sequence
    return result


def seqdict_to_fasta(seq_dict: dict, fasta_path: str):
    records = []
    for k, v in seq_dict.items():
        records.append(SeqRecord(Seq(v), k, "", ""))
    SeqIO.write(records, fasta_path, 'fasta')


def psl_to_df(psl_path: str) -> pd.DataFrame:
    columns = ['CHROM', 'Range', 'Score', 'Pct', 'Strand']
    df = pd.DataFrame(columns=columns)
    blat_result = SearchIO.parse(psl_path, 'blat-psl')
    
    n = 0
    for query in blat_result:
        n += 1
        for hit in query:
            for hsp in hit:
                row = {}
                row['CHROM'] = hsp.hit_id.replace('chr', '').replace('M', 'MT')
                row['Range'] = hsp.hit_range
                row['Score'] = hsp.score
                row['Pct'] = hsp.ident_pct
                df = df.append(row, ignore_index=True)
    
    df = df.astype({'CHROM': str, 'Score': int, 'Pct': float})
    df = df.groupby(['CHROM', 'Range']).agg({'Score': 'mean', 'Pct': 'mean', 'CHROM': 'count'}).dropna().rename(columns={'CHROM': 'Count'})
    df = df.reset_index().sort_values(by=['Score', 'Pct'], ascending=False)
    return df


def check_interval_isin_pr(interval: pr.PyRanges, bed: pr.PyRanges) -> bool:
    if len(interval) != 1:
        raise ValueError("interval should be pr.PyRanges object with length 1")
    return interval.count_overlaps(bed).as_df().loc[0, 'NumberOverlaps'] > 0


def blat_a_sequence(sequence: str, sequence_name: str,
                    blat_path: str='/home/users/kimin/tools/blat/gfClient', 
                    reference_dir_path: str="/home/users/kimin/projects/00_Reference/",
                    port: int=2882, minScore: int=20, minIdentity: int=20) -> list:
    NUMT_tsv_path = "/home/users/kimin/projects/12_MT/00_reference/GRCh37_NumtS.tsv"
    NUMT_pr = pr.PyRanges(df=pd.read_csv(NUMT_tsv_path, sep='\t', index_col=None, usecols=[1,2,3,6]))
    input_dict = {sequence_name: sequence}
    result = []

    with tempfile.TemporaryDirectory() as temp_dir:
        query_fasta_path = f"{temp_dir}/temp.fa"
        output_psl_path = f"{temp_dir}/temp.psl"
        seqdict_to_fasta(input_dict, query_fasta_path)
        cmd = f"{blat_path} -minScore={minScore} -minIdentity={minIdentity} 127.0.0.1 {port} {reference_dir_path} {query_fasta_path} {output_psl_path}"
        subprocess.run(cmd, shell=True)

        df = psl_to_df(output_psl_path)
        if len(df.query("CHROM == 'MT'")):
            MT_top_entry = df.query("CHROM == 'MT'").iloc[0]
        else:
            MT_top_entry = {'CHROM': 'MT', 'Range': (0,0), 'Score': None, 'Pct': None}
        if len(df.query("CHROM != 'MT'")):
            NU_top_entry = df.query("CHROM != 'MT'").iloc[0]
            NU_top_pr = pr.PyRanges(chromosomes=[NU_top_entry['CHROM']], starts=[NU_top_entry['Range'][0]], ends=[NU_top_entry['Range'][1]])
            NUMT = check_interval_isin_pr(NU_top_pr, NUMT_pr)
        else:
            NU_top_entry = {'CHROM': None, 'Range': (0,0), 'Score': None, 'Pct': None}
            NUMT = False
        Top_three_string = ','.join([f"{getattr(row, 'CHROM')}:{getattr(row, 'Range')[0]}-{getattr(row, 'Range')[1]} {getattr(row, 'Score'):.2f} {getattr(row, 'Pct'):.2f}" for row in df.iloc[0:3].itertuples()])

        result.append(NUMT)
        result.append(f"{NU_top_entry['CHROM']}:{NU_top_entry['Range'][0]}-{NU_top_entry['Range'][1]}")
        result.append(NU_top_entry['Score'])
        result.append(NU_top_entry['Pct'])
        result.append(f"{MT_top_entry['CHROM']}:{MT_top_entry['Range'][0]}-{MT_top_entry['Range'][1]}")
        result.append(MT_top_entry['Score'])
        result.append(MT_top_entry['Pct'])
        result.append(Top_three_string)

    return result

@cached(cache={}, key=lambda umis_key, umis_set: hashkey(umis_key))

def blat_umis_reads(umis_key: str, umis_set: dict , 
                    blat_path: str='/home/users/kimin/tools/blat/gfClient', 
                    reference_dir_path: str="/home/users/kimin/projects/00_Reference/",
                    port: int=2882, minScore: int=20, minIdentity: int=50) -> list:

    result = []
    ids = umis_key.split(":")
    bam_id = ids[0]
    umis_id = (int(ids[1]), int(ids[2]))
    orientations = list(umis_set[bam_id][umis_id].keys())
    num_of_reads = len(umis_set[bam_id][umis_id][orientations[0]][0])
    if num_of_reads > 2:
        raise ValueError(f"{umis_id} {orientations[0]} first fragment is neither paired nor unique")

    NUMT_tsv_path = "/home/users/kimin/projects/12_MT/00_reference/GRCh37_NumtS.tsv"
    NUMT_pr = pr.PyRanges(df=pd.read_csv(NUMT_tsv_path, sep='\t', index_col=None, usecols=[1,2,3,6]))

    with tempfile.TemporaryDirectory() as temp_dir:

        for i in range(num_of_reads):
            input_dict = {}
            for ori in orientations:
                input_dict = {**input_dict,
                              **alnlist_to_seqdict(list(zip(*umis_set[bam_id][umis_id][ori]))[i], f"{umis_key}_{i}_{ori}")}

            query_fasta_path = f"{temp_dir}/{i}.fa"
            output_psl_path = f"{temp_dir}/{i}.psl"
            seqdict_to_fasta(input_dict, query_fasta_path)
            cmd = f"{blat_path} -minScore={minScore} -minIdentity={minIdentity} 127.0.0.1 {port} {reference_dir_path} {query_fasta_path} {output_psl_path}"
            subprocess.run(cmd, shell=True)

            df = psl_to_df(output_psl_path)
            df = df.query(f"Count > {len(input_dict) * 0.1}")
            if len(df.query("CHROM == 'MT'")):
                MT_top_entry = df.query("CHROM == 'MT'").iloc[0]
            else:
                MT_top_entry = {'CHROM': 'MT', 'Range': (0,0), 'Score': None, 'Pct': None}
            if len(df.query("CHROM != 'MT'")):
                NU_top_entry = df.query("CHROM != 'MT'").iloc[0]
                NU_top_pr = pr.PyRanges(chromosomes=[NU_top_entry['CHROM']], starts=[NU_top_entry['Range'][0]], ends=[NU_top_entry['Range'][1]])
                NUMT = check_interval_isin_pr(NU_top_pr, NUMT_pr)
            else:
                NU_top_entry = {'CHROM': None, 'Range': (0,0), 'Score': None, 'Pct': None}
                NUMT = False
            Top_three_string = ','.join([f"{getattr(row, 'CHROM')}:{getattr(row, 'Range')[0]}-{getattr(row, 'Range')[1]} {getattr(row, 'Score'):.2f} {getattr(row, 'Pct'):.2f}" for row in df.iloc[0:3].itertuples()])

            result.append(NUMT)
            result.append(f"{NU_top_entry['CHROM']}:{NU_top_entry['Range'][0]}-{NU_top_entry['Range'][1]}")
            result.append(NU_top_entry['Score'])
            result.append(NU_top_entry['Pct'])
            result.append(f"{MT_top_entry['CHROM']}:{MT_top_entry['Range'][0]}-{MT_top_entry['Range'][1]}")
            result.append(MT_top_entry['Score'])
            result.append(MT_top_entry['Pct'])
            result.append(Top_three_string)

    if num_of_reads == 1:
        result = result * 2

    print_err(f"{umis_key} done")

    return result