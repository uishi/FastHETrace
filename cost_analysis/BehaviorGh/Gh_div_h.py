import math
import numpy as np
import sys
import matplotlib.pyplot as plt

def g(h, M):
    fl = int(M / h)
    pof2 =  1 << fl
    r = M % h
    return pof2 * (h + r)

def G(h, M):
    return g(h, M) - h

if __name__ == "__main__" :
    argc = len(sys.argv)
    if argc != 2:
        print("python3 ./xxxx M ")
        exit(1)

    logN = 15
    M = int(sys.argv[1])
    print("M = {}".format(M))
    prev_gh  = g(1,   M)
    for h in range(2, M+1):
        print("G({})/{} = {}".format(h, h, G(h, M)/h))

    print()
