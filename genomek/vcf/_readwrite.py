import os
import sys
import re
import io
import gzip
import tempfile
import subprocess
import pysam
import cyvcf2
import pandas as pd
from pandas.api.types import CategoricalDtype

CSQ_columns = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|NEAREST|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE"

def read_vcf(vcf_path:str, hg38:bool = False) -> pd.DataFrame:
    '''
    Input: vcf_path in string. supports both .gz and .vcf
    Output: pandas dataframe.
    hg38: default false. if set True, different chromosome strings will be used
    
    Easiest way to load vcf files. Does not support format parsing.
    CHROM column is converted to categorical dtype, and NAs in CHROM columns are removed.
    
    adopted from https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
    '''
    if hg38:
        chromosomes = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y', 'M']]
    else:
        chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]

    chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)
    
    if vcf_path.endswith('.gz'):
        with gzip.open(vcf_path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
        
    else:
        with open(vcf_path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

    df = pd.read_csv(io.StringIO(''.join(lines)),
                     dtype={'#CHROM': chrom_cat_type, 'POS': int, 
                            'ID': str, 'REF': str, 'ALT': str,
                            'QUAL': str, 'FILTER': str, 'INFO': str},
                     sep='\t').rename(columns={'#CHROM': 'CHROM'})

    df = df.loc[~df['CHROM'].isna()]
    return df




def _vcftotsv(vcf_path:str = "", info=None, format=None, sample=None, human=False, na="NA") -> pd.DataFrame:
    '''
    vcf: Input vcf file path. May be vcf, vcf.gz, bcf, or bcf.gz.
    info: Comma-separated list of INFO fields to include. If not set, all fields are included.
    format: Comma-separated list of FORMAT fields to include. If not set, all fields are included.
    sample: Comma-separated list of sample names to include. If not set, all samples are included.
    human: Field name is prepended to each field content.
    na: String corresponding to missing value. Default: NA
    '''
    # check vcf path and load vcf
    if vcf_path == "" or not os.path.isfile(vcf_path):
        print('vcf file does not exist')
        sys.exit(1)
    else:
        vcf = cyvcf2.VCF(vcf_path)

    # missing value character
    NAchar = na

    with tempfile.NamedTemporaryFile(mode='w') as ofile:
    # with open('temp.vcf', 'w') as ofile:
        # get default INFO and FORMAT field names and sample names
        INFOfields = []
        FORMATfields = []
        for i in vcf.header_iter():
            if i['HeaderType'] == 'INFO':
                INFOfields.append(i['ID'])
            if i['HeaderType'] == 'FORMAT':
                FORMATfields.append(i['ID'])

        samples_all = vcf.samples
        samples = samples_all

        # change default values if arguments are set
        if info != None:
            INFOfields = info.split(',')
        if format != None:
            FORMATfields = format.split(',')
        if sample != None:
            samples = sample.split(',')

        # get indices of sample names
        samples_idx = [ samples_all.index(x) for x in samples ]
        
        # print meta lines
        raw_header = vcf.raw_header.strip().split('\n')
        for line in raw_header[:-1]:
            ofile.write(line + '\n')

        # print header line
        raw_header_split = raw_header[-1].split('\t')
        header = list()
        header.extend(raw_header_split[:7])
        for key in INFOfields:
            header.append(key)
        if len(raw_header_split) > 8:
            for key in FORMATfields:
                for sample in samples:
                    header.append(f'{key}:::{sample}')

        # Added var_type by kimin 20201228
        header.append('VAR_TYPE')
        ofile.write('\t'.join(header) + '\n')

        # set INFO parsing regex pattern
        pat_INFO = re.compile('^([^=]+)(=(.+))?$')

        # iterate over variant lines
        for v in vcf:
            # print(v)
            raw_str = str(v).strip().split('\t')

            tmp = list()
            fixed = [ x if x != '.' else NAchar for x in raw_str[:7] ]
            tmp.extend( fixed )

            raw_INFOdict = dict()
            for x in raw_str[7].split(';'):
                mat = pat_INFO.match(x)
                if mat.group(2) == None:
                    raw_INFOdict[x] = 'TRUE'
                else:
                    raw_INFOdict[mat.group(1)] = mat.group(3)

            for key in INFOfields:
                if key in raw_INFOdict.keys():
                    val = raw_INFOdict[key]
                    if val == '.':
                        val = NAchar

                    if human:
                        tmp.append(f'{key}={val}')
                    else:
                        tmp.append(f'{val}')
                else:
                    if human:
                        tmp.append(f'{key}={NAchar}')
                    else:
                        tmp.append(f'{NAchar}')

            if len(raw_str) > 8: # FORMAT columns are included

                # add FORMAT fields
                raw_FORMATfields = raw_str[8].split(':')
                raw_FORMATvalues = [ x.split(':') for x in raw_str[9:] ]

                for key in FORMATfields:
                    if key in raw_FORMATfields:
                        key_idx = raw_FORMATfields.index(key)
                        for idx, sample in zip(samples_idx, samples):
                            val = raw_FORMATvalues[idx][key_idx]
                            if val == '.':
                                val = NAchar

                            if human:
                                tmp.append(f'{key}:::{sample}={val}')
                            else:
                                tmp.append(f'{val}')
                    else:
                        for sample in samples:
                            if human:
                                tmp.append(f'{key}:::{sample}={NAchar}')
                            else:
                                tmp.append(f'{NAchar}')

            # Added var_type by kimin 20201228
            tmp.append(v.var_type)

            ofile.write('\t'.join(tmp) + '\n')
        
        ofile.flush()    

        df = read_vcf(ofile.name)

    return df 


def _tsvtovcf(tsv, output, vcf=None, f=None, O='z', na='.', d=':::'):
    '''
    Core function for tidyvcf.tsvtovcf method
    tsv: input tsv file path
    output: output file path
    vcf: template vcf path from which meta lines are copied
    O: output file format as in bcftools common output. Default is "z"
    na: na character. Default is "."
    d: delimiter between format field name and sample name. Default is ":::"
    '''

    infile = open(tsv, 'r')
    ofile_path = output
    NAchar = na
    delim = d
    ofile_format = O
    template_vcf = cyvcf2.VCF(vcf) if vcf != None else None
    fasta = pysam.FastaFile(f) if f != None else None
    vcfVer = 'VCFv4.3'

    # parse input tsv file format
    for line in infile:
        if line.startswith('##'):
            continue

        if line.startswith('#'):
            linesp = line.replace('\n', '').split('\t')
            if len(linesp) >= 7 and linesp[0:7] == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']:
                mode = 1
                input_hdr = linesp[7:]
                break
            elif len(linesp) >= 5 and linesp[0:5] == ['#CHROM', 'POS', 'ID', 'REF', 'ALT']:
                mode = 2
                input_hdr = linesp[5:]
                break
            else:
                printErr('Input tsv file has an invalid format.')
                sys.exit(1)

    # parse INFO, FORMAT, and sample names from input_hdr
    delim_pat = re.compile(delim)
    # these can be empty lists
    INFO_fields = [ x for x in input_hdr if delim_pat.search(x) == None ]
    FORMAT_fields_pre = [ x.split(delim) for x in input_hdr if delim_pat.search(x) != None ]
    FORMAT_fields = list( set( [ x[0] for x in FORMAT_fields_pre ] ) )
    samples = list( set( [ x[1] for x in FORMAT_fields_pre ] ) )
    INFO_fields_idx = [ (idx, x) for (idx, x) in enumerate(input_hdr) if delim_pat.search(x) == None ]
    FORMAT_fields_idx = [ (idx, x.split(delim)[0], x.split(delim)[1]) for (idx, x) in enumerate(input_hdr) if delim_pat.search(x) != None ]

    tmpfile_dir = os.path.dirname(ofile_path)
    tmp_vcf = tempfile.NamedTemporaryFile(mode = 'w+', dir = tmpfile_dir)

    # write meta lines
    if template_vcf != None:
        raw_hdr = template_vcf.raw_header.strip().split('\n')[:-1] # only meta lines
        for line in raw_hdr:
            tmp_vcf.write(line + '\n')
    else:
        tmp_vcf.write(f'##fileformat={vcfVer}\n')
        if f != None:
            for chrom, length in zip(fasta.references, fasta.lengths):
                tmp_vcf.write(f'##contig=<ID={chrom},length={length}>\n')
        for field in INFO_fields:
            tmp_vcf.write(f'##INFO=<ID={field},Number=.,Type=String,Description="">\n')
        for field in FORMAT_fields:
            tmp_vcf.write(f'##FORMAT=<ID={field},Number=.,Type=String,Description="">\n')


    # write header line
    output_hdr = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if len(FORMAT_fields) > 0:
        output_hdr.append('FORMAT')
        for sample in samples:
            output_hdr.append(sample)
    tmp_vcf.write('\t'.join(output_hdr) + '\n')


    # loop over tsv variant lines
    for line in infile:
        linesp = line.replace('\n', '').split('\t')
        if mode == 1:
            fixed = linesp[:7]
            input_cols = linesp[7:]
        elif mode == 2:
            fixed = linesp[:5] + ['.', '.']
            input_cols = linesp[5:]
        fixed = [ '.' if x == NAchar else x for x in fixed ]

        # set initial values
        INFO_pre = dict()
        FORMAT_pre = dict()
        for field in FORMAT_fields:
            FORMAT_pre[field] = dict()

        # get INFO contents
        for idx, field in INFO_fields_idx:
            INFO_pre[field] = input_cols[idx]
        INFO = { k:v for k,v in INFO_pre.items() if v != NAchar }
        
        # get FORMAT contents
        for idx, field, sample in FORMAT_fields_idx:
            FORMAT_pre[field][sample] = input_cols[idx]

        FORMAT_cols = list()
        for field in FORMAT_fields:
            for x in FORMAT_pre[field].values():
                if x != NAchar:
                    FORMAT_cols.append(field)
                    break
        if 'GT' in FORMAT_cols:
            FORMAT_cols.remove('GT')
            FORMAT_cols = ['GT'] + FORMAT_cols

        FORMAT = dict()
        for sample in samples:
            FORMAT[sample] = list()
        for field in FORMAT_cols:
            for sample in samples:
                val = '.' if FORMAT_pre[field][sample] == NAchar else FORMAT_pre[field][sample]
                FORMAT[sample].append(val)

        if not 'GT' in FORMAT_cols:
            FORMAT_cols = ['GT'] + FORMAT_cols
            for sample in samples:
                FORMAT[sample] = ['./.'] + FORMAT[sample]

        # write results
        result = fixed
        if len(INFO_fields) == 0:
            result.append('.')
        else:
            result.append( ';'.join([ f'{k}={v}' for k,v in INFO.items() ]) )
            if len(FORMAT_fields) > 0:
                result.append( ':'.join(FORMAT_cols) )
                for sample in samples:
                    result.append( ':'.join(FORMAT[sample]) )

        tmp_vcf.write('\t'.join(result) + '\n')

    tmp_vcf.seek(0)

    # convert tmp_vcf with bcftools
    subprocess.run(f'bcftools view {tmp_vcf.name} -O {ofile_format} -o {ofile_path}', shell=True)

    # subprocess.run(['bcftools', 'view', tmp_vcf.name, '-O', ofile_format, '-o', ofile_path])
    if ofile_format in 'zb':
        subprocess.run(['bcftools', 'index', ofile_path])



class tidydf:
    def __init__(self, vcf_path, info=None, format=None, sample=None, human=False, na="NA"):
        self.df = _vcftotsv(vcf_path=vcf_path, info=info, format=format, sample=sample, human=human, na=na)
        self.vcf_path = os.path.abspath(vcf_path)
    def tsvtovcf(self, output_path="", output_format="z"):
        print("Warning: the new vcf does not have proper data type specification at the header")

        # Index original vcf file just in case
        cmd = f'bcftools index {self.vcf_path}'    
        subprocess.run(cmd, shell=True)

        # Write df to tsv at a temporary file
        with tempfile.TemporaryDirectory() as tempdir:
            tmp_tsv_path = f"{tempdir}/tmp.tsv"

            # Write df to tsv
            try: 
                df_to_save = self.df.drop('VAR_TYPE', axis=1).copy()
            except: 
                df_to_save = self.df.copy()
            df_to_save.rename(columns={'CHROM': '#CHROM'}).to_csv(tmp_tsv_path, sep='\t', na_rep='.', index=False)

            # Write tsv to vcf
            _tsvtovcf(tsv=tmp_tsv_path, output=output_path, vcf=self.vcf_path, O=output_format, na='.', d=':::')



#########################
## Deprecated Functions##
#########################
def df2csv(df,fname,formats=[],sep='\t'):
    """
    adopted from sigprofiler 진짜 빠른진 잘 모르겠음
    # function is faster than to_csv
    # 7 times faster for numbers if formats are specified, 
    # 2 times faster for strings.
    # Note - be careful. It doesn't add quotes and doesn't check
    # for quotes or separators inside elements
    # We've seen output time going down from 45 min to 6 min 
    # on a simple numeric 4-col dataframe with 45 million rows.
    """
    if len(df.columns) <= 0:
        return
    Nd = len(df.columns)
    Nd_1 = Nd 
    #formats = myformats[:] # take a copy to modify it
    Nf = 0
    formats.append('%s')
    # make sure we have formats for all columns
    if Nf < Nd:
        for ii in range(Nf,Nd, 1):
            coltype = df[df.columns[ii]].dtype
            ff = '%s'
            if coltype == np.int64:
                ff = '%d'
            elif coltype == np.float64:
                ff = '%f'
            formats.append(ff)
    fh=open(fname,'w', buffering=200000)
    header = ['MutationType'] + list(df.columns)
    fh.write('\t'.join(header) + '\n')
    for row in df.itertuples(index=True):
        ss = ''
        for ii in range(0,Nd+1,1):
            ss += formats[ii] % row[ii]
            if ii < Nd_1:
                ss += sep
        fh.write(ss+'\n')
    fh.close()


    
def get_arrays(vcf_path: str, *args, only_filter=False, **kwargs) -> dict:
    '''
    Prototype work. Uses cyvcf to parse a vcf file and returns a dictionary of arrays.

    Example usage
    get_arrays(vcf_path, "CHROM", "POS, "REF", DP_int32_1="INFO", AF_float64_1="FORMAT")
    results = scripts.get_arrays(v, "CHROM", "POS", "REF", "ALT", "var_type", only_filter=True, AF_float64_1="FORMAT", DP_int32_1="FORMAT", AD_int32_2="FORMAT")
    '''
    vcf = cyvcf2.VCF(vcf_path)
    sample_list = [x.replace('-', '_') for x in vcf.samples]
    len_samples = len(sample_list)
    if only_filter: len_vcf = sum(1 for variant in cyvcf2.VCF(vcf_path) if variant.FILTER is None)
    else: len_vcf = sum(1 for variant in cyvcf2.VCF(vcf_path))
    result_dict_1 = {}
    result_dict_2 = {}
    result_dict_3 = {}
    result_dict_4 = {}
    result_dict_5 = {}
    result_dict_6 = {}

    if "CHROM" in args: result_dict_1["CHROM"] = np.empty(len_vcf, dtype=np.dtype('U10'))
    if "POS" in args: result_dict_1["POS"] = np.empty(len_vcf, dtype=np.int64)
    if "REF" in args: result_dict_1["REF"] = np.empty(len_vcf, dtype=np.dtype(('U', 32)))
    if "var_type" in args: result_dict_1["var_type"] = np.empty(len_vcf, dtype=np.dtype(('U', 5)))
    if "ALT" in args: result_dict_2["ALT"] = np.empty(len_vcf, dtype=np.dtype(('U', 32)))
    if "FILTER" in args: result_dict_3["FILTER"] = np.empty(len_vcf, dtype=np.dtype(bool))

    for k, v in kwargs.items():
        if v == "INFO":
            key, dtype, shape = k.split('_')
            result_dict_4[v+'_'+key] = np.empty(len_vcf, dtype=np.dtype((dtype, int(shape))))
        elif v == "FORMAT":
            key, dtype, shape = k.split('_')
            result_dict_5[v+'_'+key] = np.empty(len_vcf, dtype=np.dtype((dtype, (len_samples, int(shape)))))

    i = 0
    for variant in vcf:
        if only_filter and variant.FILTER is not None: 
            continue
        for k, v in result_dict_1.items():
            v[i] = getattr(variant, k)
        for k, v in result_dict_2.items():
            v[i] = ','.join(getattr(variant, k))
        for k, v in result_dict_3.items():
            v[i] = getattr(variant, k) is None
        for k, v in result_dict_4.items():
            INFO, KEY = k.split('_')
            DICT = getattr(variant, INFO)
            v[i] = DICT.get(KEY, 0)
        for k, v in result_dict_5.items():
            METHOD, KEY = k.split('_')
            METHOD = getattr(variant, METHOD.lower())
            v[i] = METHOD(KEY)
        i += 1
    
    for k, v in result_dict_5.items():
        for i, sample in enumerate(sample_list):
            if v.shape[2] == 1:
                new_k = k + '_' + sample    
                result_dict_6[new_k] = v[:,i,:].ravel()
            else: 
                for j in range(v.shape[2]):
                    new_k = k + str(j) + '_' + sample
                    result_dict_6[new_k] = v[:,i,j].ravel()

    result = {**result_dict_1, **result_dict_2, **result_dict_3, **result_dict_4, **result_dict_6}

    return result



def get_arrays_specific(vcf_path: str, *args, **kwargs) -> dict:
    '''
    Prototype work. Uses cyvcf to parse a vcf file and returns a dictionary of arrays.
    Unfinished work.
    '''
    if "CHROM" in args: 
        chrom_array = np.empty(len_vcf, dtype=np.dtype('U4'))
    if "POS" in args:
        pos_array = np.empty(len_vcf, dtype=np.int64)
    if "REF" in args:
        ref_array = np.empty(len_vcf, dtype=np.dtype('U32'))
    if "ALT" in args:
        alt_array = np.empty(len_vcf, dtype=np.dtype('U32'))
    if "VAF1" in args:
        vaf_array1 = np.empty(len_vcf, dtype=np.float32)
    if "VAF2" in args:
        vaf_array2 = np.empty(len_vcf, dtype=np.float32)

        if "CHROM" in args:
            chrom_array[i] = variant.CHROM
        if "POS" in args:
            pos_array[i] = variant.start
        if "REF" in args:
            ref_array[i] = variant.REF
        if "ALT" in args:
            alt_array[i] = variant.ALT[0]
        if "VAF1" in args:
            vaf_array1[i] = variant.format("AF")[0,0]
        if "VAF2" in args:
            vaf_array2[i] = variant.format("AF")[1,0]
    
    result['CHROM'] = chrom_array
    result['POS'] = pos_array
    result['REF'] = ref_array
    result['ALT'] = alt_array
    result['VAF1'] = vaf_array1
    result['VAF2'] = vaf_array2
    
    return result