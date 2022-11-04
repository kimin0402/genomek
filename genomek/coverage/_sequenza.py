import collections

def get_sequenza_seqz(chrom, pos, ref, alt, nbam, tbam):
    npup = nbam.pileup(chrom, pos-1, pos, truncate = True, min_base_quality=20, min_mapping_quality=20)
    try:
        npupcol = next(npup)
        ndp = npupcol.nsegments
    except StopIteration:
        ndp = 0

    tpup = tbam.pileup(chrom, pos-1, pos, truncate = True, min_base_quality=20, min_mapping_quality=20)
    try:
        tpupcol = next(tpup)
        tdp = tpupcol.nsegments
        tdp_edit = tpupcol.get_num_aligned()
    except StopIteration:
        tdp = 0
        tdp_edit = 0
    
    if ndp == 0 or tdp_edit == 0:
        return [None] * 13

    dr = round(tdp/ndp, 3)
    bases = collections.Counter( tpupcol.get_query_sequences() )
    good = tdp

    zyg = 'het'
    AB_tumor = '.'
    tumor_strand = 0
    vaf_ref = round( ( bases[ref.upper()] + bases[ref.lower()] ) / tdp_edit, 3 )
    vaf_alt = round( ( bases[alt.upper()] + bases[alt.lower()] ) / tdp_edit, 3 )
    if max(vaf_ref, vaf_alt) == vaf_ref:
        Af = vaf_ref
        Bf = vaf_alt
        AB_normal = ref + alt
    else:
        Af = vaf_alt
        Bf = vaf_ref
        AB_normal = alt + ref

    # print(chrom, pos, ref, ndp, tdp, dr, Af, Bf, zyg, gc, good, AB_normal, AB_tumor, tumor_strand)
    return chrom, pos, ref, ndp, tdp, dr, Af, Bf, zyg, good, AB_normal, AB_tumor, tumor_strand