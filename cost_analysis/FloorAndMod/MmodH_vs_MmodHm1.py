import math
import numpy as np
import sys


if __name__ == "__main__" :
    argc = len(sys.argv)
    if argc != 2:
        print("python3 ./xxxx M ")
        exit(1)

    logN = 15
    M = int(sys.argv[1])
    print("M = {}".format(M))
    for h in range(2, M+1):
        print("h ={}".format(h), end=",")
        mmh = int(M / h)
        mmhm1 = int(M / (h-1))
        desc=""
        if mmh > mmhm1:
            desc = ">"
        elif mmh < mmhm1:
            desc = "<"
        else:
            desc = "="

        print("(FL(M / h), prd, FL(M / h-1)) = ({}  {}  {})".format(mmh, desc, mmhm1), end=" ")

        mmh = M % h
        mmhm1 = M % (h-1)
        desc=""
        if mmh > mmhm1:
            desc = ">"
        elif mmh < mmhm1:
            desc = "<"
        else:
            desc = "="

        print("(M mod h, prd, M mod h-1) = ({}  {}  {})".format(mmh, desc, mmhm1))

