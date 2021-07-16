import pysam


## Small functions from sigprofiler
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','[':'[',']':']','>':'>'}[B] for B in x][::-1])
revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N'}[B] for B in x][::-1])



def context(fasta, CHROM, POS, ALT, cont=3):
    context_start = int(cont/2) 
    ref_context = fasta.fetch(CHROM, POS - context_start - 1, POS - context_start - 1 + cont)
    if (ref_context[context_start] == 'A') or (ref_context[context_start] == 'G'):
        ref_context = revcompl(ref_context)
        ALT = revcompl(ALT)
    subsitution = ref_context[context_start] + ">" + ALT
    return ref_context, subsitution


def context_3(fasta, CHROM, POS, ALT, cont=3):
    '''
    fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")
    df_tumor_blood = df_tumor_blood.astype({"CHROM":'str'})
    df_tumor_blood[['context_3', 'substitution']] = df_tumor_blood.query("var_type == 'snp'").apply(lambda x: scripts.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')
    pandas_count = df_tumor_blood[['context_3', 'substitution']].groupby('context_3').count()
    pandas_count = sorted(df_tumor_blood['context_3'].value_counts().index, key=lambda x: (x[1], x[0], x[2]))
    pandas_count['96'] = pandas_count.index
    '''
    context_start = int(cont/2)
    CHROM = str(CHROM)
    ref_context = fasta.fetch(CHROM, POS - context_start - 1, POS - context_start - 1 + cont)
    if (ref_context[context_start] == 'A') or (ref_context[context_start] == 'G'):
        ref_context = revcompl(ref_context)
        ALT = revcompl(ALT)
    subsitution = ref_context[context_start] + ">" + ALT
    ref_context = ref_context[0] + "[" + subsitution + "]" + ref_context[2]
    return ref_context, subsitution


def context_3_strand(fasta, CHROM, POS, ALT, cont=3):
    '''
    fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")
    df_tumor_blood = df_tumor_blood.astype({"CHROM":'str'})
    df_tumor_blood[['context_3', 'substitution']] = df_tumor_blood.query("var_type == 'snp'").apply(lambda x: scripts.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')
    pandas_count = df_tumor_blood[['context_3', 'substitution']].groupby('context_3').count()
    pandas_count = sorted(df_tumor_blood['context_3'].value_counts().index, key=lambda x: (x[1], x[0], x[2]))
    pandas_count['96'] = pandas_count.index
    '''
    context_start = int(cont/2)
    CHROM = str(CHROM)
    ref_context = fasta.fetch(CHROM, POS - context_start - 1, POS - context_start - 1 + cont)
    if (ref_context[context_start] == 'A') or (ref_context[context_start] == 'G'):
        ref_context = revcompl(ref_context)
        ALT = revcompl(ALT)
        strand = 'Heavy'
    else:
        strand = 'Light'
    subsitution = ref_context[context_start] + ">" + ALT
    ref_context = ref_context[0] + "[" + subsitution + "]" + ref_context[2]
    return ref_context, subsitution, strand

