import sys

def h_plus_m_mod_h(h:int, m:int):
    return h+(m%h)

def gh_left(h:int, m:int):
    floor = int(m/h)
    pof2 = 1<< floor
    return pof2 * h
    
def gh_right(h:int, m:int):
    floor = int(m/h)
    pof2 = 1<< floor
    return pof2 * (m%h)

if __name__ == "__main__" :
    M = int(sys.argv[1])
    argc = len(sys.argv)
    print(" h + (M mod h)")
    for h in range(2, M+1):
        ghm1 = h_plus_m_mod_h(h-1, M)
        gh   = h_plus_m_mod_h(h, M)
        print(gh <= ghm1)

    print(" 2^{M/h} h")
    for h in range(2, M+1):
        ghm1 = gh_left(h-1, M)
        gh   = gh_left(h, M)
        print(gh <= ghm1)

    print(" 2^{M/h} h")
    for h in range(2, M+1):
        ghm1 = gh_right(h-1, M)
        gh   = gh_right(h, M)
        print(gh <= ghm1)
