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
    for h in range(2, M+1):
        gh  = g(h,   M)
        gh1 = g(h-1, M)
        bound = gh1 - gh + 1
        #print(" C({}) <= C({}) when CH/CIP<= {}".format(h, h-1, bound))
        print("h={}: C({}) < C({}) when CH/CIP > {}".format(h-1, h-1, h, bound))
    print("h={}".format(M))
