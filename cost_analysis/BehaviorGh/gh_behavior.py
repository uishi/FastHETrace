import math
import numpy as np
import sys
import matplotlib.pyplot as plt

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
    prev_gh  = g(1,   M)
    for h in range(2, M+1):
        gh  = g(h,   M)
        diff = prev_gh - gh
        print("\tg({}) - g({}): {}".format(h-1, h, diff))
        if prev_gh < gh:
            print("non-decreasing!!!")
            break
        #print(" C({}) <= C({}) when CH/CIP<= {}".format(h, h-1, bound))
        #print(" g({}) = {}".format( h, gh))
        prev_gh = gh
    print("monotonical decrease : )")

    print("\n\ng(h) values:")
    for h in range(1, M+1):
        gh  = g(h,   M)
        print("\tg({}) = {}".format(h, gh))

    hs = list( range(1, M+1) )
    ghs = [g(h, M) for h in hs]
    plt.plot(hs, ghs)
    plt.show()
