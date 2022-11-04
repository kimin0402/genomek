
def lsc(s: str, words: str='ACGT') -> float:
    '''
    Input
    s: sequence string
    words: set of bases for calculation
    
    Output
    lingustic sequence complexity in float
    
    Reference
    PMID: 12050064 
    '''
    
    l = len(words) # word size
    m = len(s) # sequence size
    max_vocab_list = [min(l**k, m-k+1) for k in range(1, m+1)]
    actual_vocab_list = [len({s[i:i+j] for i in range(0, m-j+1)}) for j in range(1,m+1)]
        
    return sum(actual_vocab_list)/sum(max_vocab_list)


def lsc10(s: str, words: str='ACGT') -> float:
    '''
    Input
    s: sequence string
    words: set of bases for calculation
    
    Output
    lingustic sequence complexity in float
    
    Reference
    PMID: 12050064 
    '''
    
    l = len(words) # word size
    m = len(s) # sequence size
    w = 10 # window size
    max_vocab_list = [min(l**k, m-k+1) for k in range(1, w+1)]
    actual_vocab_list = [len({s[i:i+j] for i in range(0, m-j+1)}) for j in range(1, w+1)]
        
    return sum(actual_vocab_list)/sum(max_vocab_list)



def sdust_score(s: str) -> float:
    def occurrences(string, sub):
        count = start = 0
        while True:
            start = string.find(sub, start) + 1
            if start > 0:
                count+=1
            else:
                return count
    d = {}
    for k in [f"{x}{y}{z}" for x in 'ACTG' for y in 'ACTG' for z in 'ACTG']:
        d[k] = occurrences(s, k)
    return sum(v*(v-1)/2 for v in d.values())/(len(s)-3) * 10