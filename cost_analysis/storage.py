import math, sys
import numpy as np
import matplotlib.pyplot as plt
import collections
import math
import timeit
import os

import HKcost as HK
import Util as u
import cost_analysis_helper as CostHelper

colors = ["k", "r", "b", "g"]

def plot_storage(logN, M, Ms, res, prefixname):
    filedir = "/tmp/result/numkey"
    os.system("mkdir -p {}".format(filedir))

# log scale
    cnt = 0
    for M in Ms:
        print("M = {}, costs = {}".format(M, res[M].costs))
        hs = list(range(1, M + 1))
        plt.plot(hs, res[M].num_ips_unroll,      label="M = {}".format( M) ,color=colors[cnt])
        #plt.plot(hs, res[M].num_ips_unroll,      label="M = {}".format( M))
        cnt += 1

    plt.legend()
    plt.xlabel("h"   , fontsize=22)
    plt.yscale("log")
    #plt.ylabel("Cost (log-scaled)", fontsize=22)
    plt.ylabel("# Evaluation Keys (log-scaled)", fontsize=22)
    plt.xticks(np.arange(1, M+1, step=1))  
    plt.savefig("{}/{}_numkey_logN{}_logscale.eps".format(filedir, prefixname, logN))
    #plt.show()
    plt.clf()

# non-log scale (if M > 7, h=1, 2 cutoff)
    cnt = 0
    for M in Ms:
        print("M = {}, costs = {}".format(M, res[M].costs))
        if M > 7:
            hs = list(range(3, M + 1))
            plt.plot(hs, res[M].num_ips_unroll[2:],      label="M = {}".format( M) ,color=colors[cnt])
        else:
            hs = list(range(1, M + 1))
            plt.plot(hs, res[M].num_ips_unroll,      label="M = {}".format( M) ,color=colors[cnt])
        #plt.plot(hs, res[M].num_ips_unroll,      label="M = {}".format( M))
        cnt += 1

    plt.legend()
    plt.xlabel("h"   , fontsize=22)
    #plt.ylabel("Cost (log-scaled)", fontsize=22)
    plt.ylabel("# Evaluation Keys", fontsize=22)
    plt.xticks(np.arange(1, M+1, step=1))  
    plt.savefig("{}/{}_numkey_logN{}.eps".format(filedir, prefixname, logN))
    #plt.show()
    plt.clf()

if __name__ == "__main__" :
    # max
    argc = len(sys.argv)
    if argc != 3:
        print("python3 ./xxx [logN] [PrefixName]")
        exit(1)

    logN   = int(sys.argv[1])
    prefixname = sys.argv[2]
    # Depth, k, kscost choice do not matter for storage cost as long as compressed evks are not supported
    depth  = 19
    k      = 20
    ks_id  = 1
    ell = depth
    num_rns_ctxt =  ell + 1

    is_coeff   = 0
    if is_coeff:
        print("Coefficient Rep")
    else:
        print("Eval Rep")

    if num_rns_ctxt % k != 0:
        print("wrong")
        exit(1)

    if logN < 10 or logN > 16:
        print("logN must be 13 14 15")
        exit(1)

    print("NOTE: This script assumes Eval-Eval (Input-Output) nalysis")


    if (ell + 1 ) %k != 0:
        print("# mod + 1 must be divisible by #Special mod")
        exit(1)

    #print("**** We use a fixed (d, k) = ({}, {}) for M = {}".format(d, k, M))

    m_ks = HK.KSMethod.LATGOHK20
    if ks_id == 0:
        m_ks = HK.KSMethod.LATGOHK20
        print("KS = LATIGO")
    elif ks_id == 1:
        m_ks = HK.KSMethod.PALIHK20
        print("KS = PALISADE")
    else:
        print("KS id must be 0 or 1")
        exit(1)


    res = {}

    #Ms = list(range(2, logN + 1))
    Ms = [2, logN >> 1,  logN]
    for M in Ms:
        hs = list(range(1, M + 1))
        ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, ell, k, M, hs, m_ks)
        res[M] = ca.cost()
        print("M = {}, opth={}, max_speedup = {}".format(M, res[M].optimal_h, res[M].max_speedup))
    #ca.print()
    for M in Ms:
        print("M = {}, costs = {}".format(M, res[M].costs))

    plot_storage(logN, M, Ms, res, prefixname)
