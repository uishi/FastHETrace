import math
import numpy as np
import sys

def g(h, M):
    fl = int(M / h)
    pof2 =  1 << fl
    r = M % h
    return pof2 * (h + r)

if __name__ == "__main__" :
    argc = len(sys.argv)
    if argc != 2:
        print("python3 ./xxxx M ")
        exit(1)

    logN = 15
    M = int(sys.argv[1])
    print("M = {}".format(M))
    ghm = g(M, M)
    print("g(M) = {}".format(ghm))
    # upto M-1
    for h in range(1, M):
        gh  = g(h,   M)
        bound = 1 + gh - ghm
        print(" C({}) < C({}) when CH/CIP > {}".format(h, M, bound))
