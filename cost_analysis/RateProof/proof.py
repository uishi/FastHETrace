import math
import numpy as np
import sys

if __name__ == "__main__" :
    argc = len(sys.argv)
    logN = 15
    Ms = list(range(2, logN + 1))
    for M in Ms:
        #for M in range(2, logN+1):
        print("M = {}".format(M))
        for h in range(1, M):
            pow2 = 1 << math.floor(M/h)
            e = M - h
            rem = M % h
            c = pow2 * (h + rem) - (2*M)
            c /= e
            c += 1
            print(" \tC(M={}) win C(M-e={})  e = {}: if CH/CIP <= x = {:.3f},  CH/CIP > x = {:.2f} means unroll at e win!".format(M, h, e, c, c) )
