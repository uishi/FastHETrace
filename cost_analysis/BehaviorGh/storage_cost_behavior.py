import math
import numpy as np
import sys
import matplotlib.pyplot as plt

def g_storage(h, M):
    fl = int(M / h)
    pof2 =  1 << fl
    r = M % h
    return pof2 * (h + r) - h

if __name__ == "__main__" :
    argc = len(sys.argv)
    if argc != 2:
        print("python3 ./xxxx M ")
        exit(1)

    logN = 15
    M = int(sys.argv[1])
    print("M = {}".format(M))
    prev_gh  = g_storage(1,   M)
    for h in range(2, M+1):
        gh  = g_storage(h,   M)
        diff = gh - prev_gh
        print("\tg_storage({}) - g_storage({}): {}".format(h-1, h, diff))
        prev_gh = gh

    hs = list( range(1, M+1) )
    ghs = [g_storage(h, M) for h in hs]
    plt.plot(hs, ghs)
    plt.show()
