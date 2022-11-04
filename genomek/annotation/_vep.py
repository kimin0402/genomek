import io
import re
import pandas as pd
from ..tools._gadgets import chrom_cat_type_37 as chrom_cat_type
from ..tools._gadgets import fasta_37 as fasta

def df_to_vep_df(df: pd.DataFrame) -> pd.DataFrame:
    def row_operation(row):
        '''
        https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html
        '''
        chrom = row['CHROM']
        strand = '+'
        ref = row['REF']
        alt = row['ALT']
        pos = row['POS']
        if len(ref) == len(alt):
            start = pos
            end = pos + len(ref) - 1
            allele = f"{ref}/{alt}"
        elif len(ref) < len(alt):
            start = pos + 1
            end = pos
            allele = f"-/{alt[len(ref):]}"
        elif len(ref) > len(alt):
            start = pos + len(alt)
            end = pos + len(ref) - 1
            allele = f"{ref[len(alt):]}/-"
        return chrom, start, end, allele, strand
    
    result_df = pd.DataFrame()
    result_df[['chrom', 'start', 'end', 'allele', 'strand']] = df.apply(lambda x: row_operation(x), axis=1, result_type='expand')
    return result_df

def vep_df_to_df(df: pd.DataFrame, fasta) -> pd.DataFrame:
    '''
    vep_df: Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2 ->
    df: CHROM, POS, REF, ALT
    fasta: pysam fastafile
    '''
    def row_operation(row):
        chrom = row['Chromosome']
        if row['Reference_Allele'] != '-' and row['Tumor_Seq_Allele2'] != '-':
            pos = row['Start_Position']
            ref = row['Reference_Allele']
            alt = row['Tumor_Seq_Allele2']
        elif row['Reference_Allele'] == '-':
            pos = row['Start_Position']
            ref = fasta.fetch(reference=row['Chromosome'], start=row['Start_Position']-1, end=row['Start_Position'])
            alt = ref + row['Tumor_Seq_Allele2']
        elif row['Tumor_Seq_Allele2'] == '-':
            pos = row['Start_Position'] - 1
            ref = fasta.fetch(reference=row['Chromosome'], start=row['Start_Position']-2, end=row['End_Position'])
            alt = ref[0]
        return chrom, pos, ref, alt
    result_df = pd.DataFrame()
    result_df[['CHROM', 'POS', 'REF', 'ALT']] = df.apply(lambda x: row_operation(x), axis=1, result_type='expand')
    return result_df



def vep_file_to_df(vep_path: str) -> pd.DataFrame:
    '''
    path_vep = f"02_vcf/vep/{group}/{group}.output.0416.vep"
    '''
    with open(vep_path) as file:
        lines = [l for l in file if not l.startswith("##")]
    df_vep = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    df_vep.columns = [x.replace("#", "") for x in df_vep.columns]
    df_vep_dict = df_vep.to_dict(orient='index')
    for index, row in df_vep_dict.items():
        for data in row['Extra'].split(';'):
            k, v = data.split('=')
            row[k] = v
        del row['Extra']
        df_vep_dict[index] = row
    df_vep = pd.DataFrame.from_dict(df_vep_dict, orient='index')
    def row_operation(row):
        '''
        recover_chrom_pos_ref_alt
        '''
        regex = re.search(r"^([\w]+)_(\d+)_([ACTGN-]+)/([ACTGN-]+)$", row['Uploaded_variation'])
        chrom, start, ref, alt = regex.group(1), int(regex.group(2)), regex.group(3), regex.group(4)
        if row['VARIANT_CLASS'] == 'SNV' or row['VARIANT_CLASS'] == 'substitution':
            pos = start
        elif row['VARIANT_CLASS'] == 'deletion':
            pos = start - len(alt)  # actually, len(alt) is always 1 if previous vcfs were correctly normalized
            hidden_base = fasta.fetch(reference=chrom, start=pos-1, end=pos)
            ref = hidden_base+ref
            alt = hidden_base
        elif row['VARIANT_CLASS'] == 'insertion':
            pos = start - 1
            hidden_base = fasta.fetch(reference=chrom, start=pos-1, end=pos)
            ref = hidden_base
            alt = hidden_base+alt
        else:
            raise Exception("VARIANT_CLASS not defined! (something seriously wrong with the vep file)")
        return chrom, pos, ref, alt

    df_vep[['CHROM', 'POS', 'REF', 'ALT']] = df_vep.apply(lambda x: row_operation(x), axis=1, result_type='expand')
    # df_vep.drop(columns=['Uploaded_variation', 'Location', 'Allele', 'REF_ALLELE'], inplace=True)
    df_vep = df_vep[['CHROM', 'POS', 'REF', 'ALT'] + [x for x in df_vep.columns if not x in ['CHROM', 'POS', 'REF', 'ALT']]]
    df_vep.columns = ['CHROM', 'POS', 'REF', 'ALT'] + [f"VEP_{x}" for x in df_vep.columns if not x in ['CHROM', 'POS', 'REF', 'ALT']]
    df_vep = df_vep.astype({'CHROM': chrom_cat_type})
    df_vep.query("VEP_PICK == '1'", inplace=True)
    df_vep.reset_index(drop=True, inplace=True)
    return df_vep



# Prioritize Sequence Ontology terms in order of severity, as estimated by Ensembl:
# https://ensembl.org/info/genome/variation/prediction/predicted_data.html
# https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl perl script 에서 가져옴
effectPriority = {
        'transcript_ablation' : 1, # A feature ablation whereby the deleted region includes a transcript feature
        'exon_loss_variant' : 1, # A sequence variant whereby an exon is lost from the transcript
        'splice_donor_variant' : 2, # A splice variant that changes the 2 base region at the 5' end of an intron
        'splice_acceptor_variant' : 2, # A splice variant that changes the 2 base region at the 3' end of an intron
        'stop_gained' : 3, # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
        'frameshift_variant' : 3, # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
        'stop_lost' : 3, # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
        'start_lost' : 4, # A codon variant that changes at least one base of the canonical start codon
        'initiator_codon_variant' : 4, # A codon variant that changes at least one base of the first codon of a transcript
        'disruptive_inframe_insertion' : 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
        'disruptive_inframe_deletion' : 5, # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
        'conservative_inframe_insertion' : 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence between existing codons
        'conservative_inframe_deletion' : 5, # An inframe decrease in cds length that deletes one or more entire codons from the coding sequence but does not change any remaining codons
        'inframe_insertion' : 5, # An inframe non synonymous variant that inserts bases into the coding sequence
        'inframe_deletion' : 5, # An inframe non synonymous variant that deletes bases from the coding sequence
        'protein_altering_variant' : 5, # A sequence variant which is predicted to change the protein encoded in the coding sequence
        'missense_variant' : 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
        'conservative_missense_variant' : 6, # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
        'rare_amino_acid_variant' : 6, # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
        'transcript_amplification' : 7, # A feature amplification of a region containing a transcript
        'splice_region_variant' : 8, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
        'start_retained_variant' : 9, # A sequence variant where at least one base in the start codon is changed, but the start remains
        'stop_retained_variant' : 9, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
        'synonymous_variant' : 9, # A sequence variant where there is no resulting change to the encoded amino acid
        'incomplete_terminal_codon_variant' : 10, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
        'coding_sequence_variant' : 11, # A sequence variant that changes the coding sequence
        'mature_miRNA_variant' : 11, # A transcript variant located with the sequence of the mature miRNA
        'exon_variant' : 11, # A sequence variant that changes exon sequence
        '5_prime_UTR_variant' : 12, # A UTR variant of the 5' UTR
        '5_prime_UTR_premature_start_codon_gain_variant' : 12, # snpEff-specific effect, creating a start codon in 5' UTR
        '3_prime_UTR_variant' : 12, # A UTR variant of the 3' UTR
        'non_coding_exon_variant' : 13, # A sequence variant that changes non-coding exon sequence
        'non_coding_transcript_exon_variant' : 13, # snpEff-specific synonym for non_coding_exon_variant
        'non_coding_transcript_variant' : 14, # A transcript variant of a non coding RNA gene
        'nc_transcript_variant' : 14, # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
        'intron_variant' : 14, # A transcript variant occurring within an intron
        'intragenic_variant' : 14, # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
        'INTRAGENIC' : 14, # snpEff-specific synonym of intragenic_variant
        'NMD_transcript_variant' : 15, # A variant in a transcript that is the target of NMD
        'upstream_gene_variant' : 16, # A sequence variant located 5' of a gene
        'downstream_gene_variant' : 16, # A sequence variant located 3' of a gene
        'TFBS_ablation' : 17, # A feature ablation whereby the deleted region includes a transcription factor binding site
        'TFBS_amplification' : 17, # A feature amplification of a region containing a transcription factor binding site
        'TF_binding_site_variant' : 17, # A sequence variant located within a transcription factor binding site
        'regulatory_region_ablation' : 17, # A feature ablation whereby the deleted region includes a regulatory region
        'regulatory_region_amplification' : 17, # A feature amplification of a region containing a regulatory region
        'regulatory_region_variant' : 17, # A sequence variant located within a regulatory region
        'regulatory_region' :17, # snpEff-specific effect that should really be regulatory_region_variant
        'feature_elongation' : 18, # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
        'feature_truncation' : 18, # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
        'intergenic_variant' : 19, # A sequence variant located in the intergenic region, between genes
        'intergenic_region' : 19, # snpEff-specific effect that should really be intergenic_variant
        '' : 20
    }

# self made list
effectsCheck = ['splice_acceptor_variant','splice_donor_variant','transcript_ablation','exon_loss_variant'] + \
['stop_gained', 'frameshift_variant', 'protein_altering_variant', 'stop_lost'] + \
['initiator_codon_variant', 'start_lost', 'splice_region_variant'] + \
['inframe_insertion', 'inframe_deletion'] + \
['missense_variant','coding_sequence_variant','conservative_missense_variant','rare_amino_acid_variant']

def AA3to1(string):
    result = string
    aa3to1 = {'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Asx':'B', 
              'Cys':'C', 'Glu':'E', 'Gln':'Q', 'Glx':'Z', 'Gly':'G',
              'His':'H', 'Ile':'I', 'Leu':'L', 'Lys':'K', 'Met':'M',
              'Phe':'F', 'Pro':'P', 'Ser':'S', 'Thr':'T', 'Trp':'W',
              'Tyr':'Y', 'Val':'V', 'Xxx':'X', 'Ter':'*', '%3D':'='}
    for k,v in aa3to1.items():
        result = result.replace(k, v)
    return result


def GetVariantClassification(effect, var_type, inframe):
    '''
    MSKCC vcf2maf script 에서 가져옴. perl 로 돼 있어서 내가 조금 손을 봄
    https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl
    '''
    if not effect: # In case VEP was skipped
        return "Targeted_Region" 
    elif( effect in ['splice_acceptor_variant','splice_donor_variant','transcript_ablation','exon_loss_variant']):
        return "Splice_Site" 
    elif( effect == 'stop_gained'):
        return "Nonsense_Mutation" 
    elif(( effect == 'frameshift_variant' or ( effect == 'protein_altering_variant' and not inframe )) and var_type == 'DEL' ):
        return "Frame_Shift_Del" 
    elif(( effect == 'frameshift_variant' or ( effect == 'protein_altering_variant' and not inframe )) and var_type == 'INS' ):
        return "Frame_Shift_Ins" 
    elif( effect == 'stop_lost' ):
        return "Nonstop_Mutation" 
    elif( effect in ['initiator_codon_variant', 'start_lost']):
        return "Translation_Start_Site" 
    elif( effect in ['inframe_insertion'] or ( effect == 'protein_altering_variant' and inframe and var_type == 'INS' )):
        return "In_Frame_Ins" 
    elif( effect in ['inframe_deletion'] or ( effect == 'protein_altering_variant' and inframe and var_type == 'DEL' )):
        return "In_Frame_Del" 
    elif( effect in ['missense_variant','coding_sequence_variant','conservative_missense_variant','rare_amino_acid_variant']):
        return "Missense_Mutation" 
    if ( effect in ['transcript_amplification','intron_variant','INTRAGENIC','intragenic_variant']):
        return "Intron" 
    elif( effect in ['splice_region_variant', 'splice_donor_5th_base_variant', 'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant'] ):
        return "Splice_Region" 
    elif( effect in ['incomplete_terminal_codon_variant','synonymous_variant','start_retained_variant','stop_retained_variant','NMD_transcript_variant'] ):
        return "Silent" 
    elif( effect in ['mature_miRNA_variant','exon_variant','non_coding_exon_variant','non_coding_transcript_exon_variant','non_coding_transcript_variant','nc_transcript_variant'] ):
        return "RNA" 
    elif( effect in ['5_prime_UTR_variant','5_prime_UTR_premature_start_codon_gain_variant'] ):
        return "5'UTR" 
    elif( effect == '3_prime_UTR_variant' ):
        return "3'UTR" 
    elif( effect in ['TF_binding_site_variant','regulatory_region_variant','regulatory_region','intergenic_variant','intergenic_region'] ):
        return "IGR" 
    elif( effect == 'upstream_gene_variant' ):
        return "5'Flank" 
    elif( effect == 'downstream_gene_variant' ):
        return "3'Flank" 
    else:
        # Annotate everything else simply as a targeted region
        # TFBS_ablation, TFBS_amplification,regulatory_region_ablation, regulatory_region_amplification,
        # feature_elongation, feature_truncation
        return "Targeted_Region";

def df_to_maf_df(df, prefix='VEP_'):
    def row_operation(row, prefix='VEP_'):
        '''
        https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/
        https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
        '''
        Hugo_Symbol = row[f'{prefix}SYMBOL']
        Chromosome = row['CHROM']
        Reference_Allele = row[f'{prefix}REF_ALLELE']
        Tumor_Seq_Allele2 = row[f'{prefix}Allele']
        pos_list = [int(x) for x in re.findall(r"\d+", row[f'{prefix}Location'].split(":")[1])]
        if len(pos_list) == 1:
            Start_Position, End_Position = pos_list[0], pos_list[0]
        elif len(pos_list) == 2:
            pos_list = sorted(pos_list)
            Start_Position, End_Position = pos_list[0], pos_list[1]
        
        # Variant_Type
        if row[f'{prefix}VARIANT_CLASS'] == 'SNV':
            Variant_Type = 'SNP'
        elif row[f'{prefix}VARIANT_CLASS'] == 'substitution':
            if len(Reference_Allele) == len(Tumor_Seq_Allele2):
                if len(Reference_Allele) == 2:
                    Variant_Type = 'DNP'
                elif len(Reference_Allele) == 3:
                    Variant_Type = 'TNP'
                elif len(Reference_Allele) >= 4:
                    Variant_Type = 'ONP'
            else:
                raise Exception("some substitution variants have different ref, alt length")
        elif row[f'{prefix}VARIANT_CLASS'] == 'insertion':
            Variant_Type = 'INS'
        elif row[f'{prefix}VARIANT_CLASS'] == 'deletion':
            Variant_Type = 'DEL'
        
        # Variant_Classification    
        inframe = abs(len(Reference_Allele.replace("-",""))-len(Tumor_Seq_Allele2.replace("-","")))
        inframe = True if inframe % 3 == 0 else False
        effect = row[f'{prefix}Consequence']
        if effect:
            effect = sorted(effect.split(','), key=lambda x: effectPriority.get(x, 21))[0]
        else:
            effect = 'intergenic_variant'
        Variant_Classification = GetVariantClassification(effect, Variant_Type, inframe)

        # HGVSp_Short
        if row[f'{prefix}HGVSp']:
            HGVSp_Short = row[f'{prefix}HGVSp'].split(':')[1]
            HGVSp_Short = AA3to1(HGVSp_Short)
        elif effect in ['splice_donor_variant', 'splice_acceptor_variant']:
            HGVSp_Short = None
            # Something has to be done here. Check github page later https://github.com/mskcc/vcf2maf/issues/312
        else:
            HGVSp_Short = None

        return Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type, HGVSp_Short

    result_df = pd.DataFrame()
    result_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'HGVSp_Short']] = df.apply(lambda x: row_operation(x, prefix=prefix), axis=1, result_type='expand')
    result_df['IMPACT'] = df[f'{prefix}IMPACT']
    result_df['Gene_Id'] = df[f'{prefix}Gene']
    result_df['i_transcript_name'] = df[f'{prefix}Feature']
    return result_df
