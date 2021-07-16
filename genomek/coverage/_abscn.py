
def absolute_copy(cov1, cov2, C_cov, N_cov, ploidy, acf):
    '''
    cov1: cancer sample's average coverage of the bin
    cov2: normal sample's average coverage of the bin
    C_cov: cancer sample's total average depth
    N_cov: normal sample's total average depth
    ploidy: estimated ploidy from sequenza
    acf: aberrant cell fraction estimated from sequenza
    '''
    if cov2 > 0 :
        return (cov1/C_cov/cov2*N_cov-(1-acf))*ploidy/acf
    else:
        return None


    
    