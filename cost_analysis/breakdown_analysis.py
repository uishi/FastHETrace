import math, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import collections
import math
import timeit
import os

import HKcost as HK
import Util as u
import cost_analysis_helper as CostHelper

colors = ["k", "r", "b", "g"]
is_debug = False
is_png = False
eta = 0

def plot_cost(logN, M, Ms, res, k, ell, prefixname):
    filedir = "/tmp/result/cost_plot"
    os.system("mkdir -p {}".format(filedir))

    cnt = 0
    for M in Ms[logN]:
        print("M = {}, costs = {}".format(M, res[M].costs))
        hs = list(range(1, M + 1))
        plt.plot(hs, res[M].costs,      label=r"M = {}".format(M)                    ,color=colors[cnt])
        #plt.plot(hs, res[M].costs,     label=r"({}, {})".format(logN , M)                    ,color=colors[cnt])
        cnt += 1

    plt.legend(fontsize=22)
    plt.xlabel("h"   , fontsize=22)
    plt.yscale("log")
    plt.ylabel("Cost (log-scaled)", fontsize=22)
    plt.xticks(np.arange(1, M+1, step=1))  
    plt.tick_params(axis='y', which='major', labelsize=12)
    savename = "{}/{}_k{}_cost_logN{}_ell{}".format(filedir, prefixname, k, logN, ell)
    if is_png:
        plt.savefig("{}.png".format(savename))
    else:
        plt.savefig("{}.eps".format(savename))
    #plt.show()
    plt.clf()

def plot_cost_breakdown(logN, M, Ms, res,  k,ell,  prefixname):
    WIDTH = 0.9

    filedir = "/tmp/result/breakdown"
    for M in Ms[logN]:
        filedir = "/tmp/result/breakdown/logN{}M{}".format(logN, M)
        os.system("mkdir -p {}".format(filedir))
        plt.xlabel("h"   , fontsize=22)
        plt.ylabel("Cost ", fontsize=22)
        #plt.yscale("log")
        #hs = list(range(1, M+1))
        #plt.xticks(np.arange(1, M+1, step=1))  
        #data_hoist1   =  np.array(res[M].total_c_hoist1)
        #data_ip       =  np.array(res[M].total_c_ips   )
        #data_hoist2   =  np.array(res[M].total_c_hoist2)
        #plt.xticks(np.arange(2, M+1, step=1))  

        beg=2
        if M == 2:
            beg=1
        if M > 2:
            beg=3
        hs = list(range(beg, M+1))
        plt.xticks(np.arange(beg, M+1, step=1))  
        data_hoist1   =  np.array(res[M].total_c_hoist1[beg-1:])
        data_ip       =  np.array(res[M].total_c_ips   [beg-1:])
        data_hoist2   =  np.array(res[M].total_c_hoist2[beg-1:])
        print("\n\n\n M = {}".format(M))
        print("C_H1 = {}".format(data_hoist1))
        print("C_IP = {}".format(data_ip))
        print("C_H2 = {}".format(data_hoist2))
        p_hoist2   = plt.bar(hs, data_hoist2,   color="dodgerblue",  width=WIDTH)
        p_ip       = plt.bar(hs, data_ip,     bottom=data_hoist2, color="yellow",  width=WIDTH)
        p_hoist1   = plt.bar(hs, data_hoist1, bottom=data_hoist2 + data_ip, color="lime",  width=WIDTH)
        if is_debug:
            plt.title(r"$\ell = {}, k ={}, \beta={} $".format(ell, k, math.ceil((ell +1)/k)))

        #plt.legend((p_hoist1[0], p_ip[0], p_hoist2[0]), (r'$h \cdot C_{H1}$', r"$(g(h)-h)\cdot C_{IP}$", r"$h \cdot C_{H2}$"))
        if ell == 9:
            plt.ylim(pow(10, 7), pow(10, 9)*1.4)
        elif ell == 1:
            plt.ylim(pow(10, 7), pow(10, 8)*1.5)

        # https://stackoverflow.com/questions/41296313/stacked-bar-chart-with-centered-labels
        for rh2, rip, rh1, in zip(p_hoist2, p_ip, p_hoist1):
            h1  = rh1.get_height()
            hip = rip.get_height()
            h2  = rh2.get_height()

            total_cost =  h1 + hip + h2
            percent_h2 =  100 * h2 / total_cost
            percent_h1 =  100 * h1 / total_cost
            percent_ip =  100 * hip/ total_cost

            # Use curly braces for format inside of latext symbol that also uses curly braces:
            # Solution: Use double curly braces
            #   https://tex.stackexchange.com/questions/139643/string-formatting-when-using-pythontex
            color_str = "black"
            font_s = 8
            if ell == 9:
                #plt.text(rh2.get_x() + rh2.get_width() / 2., h2/3,            r"$\mathbb{{{:.1f}}}$".format(percent_h2), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                #plt.text(rip.get_x() + rip.get_width() / 2., h2 + hip/10,      r"$\mathbb{{{:.1f}}}$".format(percent_ip), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                #plt.text(rh1.get_x() + rh1.get_width() / 2., h2 + hip + h1/2., r"$\mathbb{{{:.1f}}}$".format(percent_h1), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                plt.text(rh2.get_x() + rh2.get_width() / 2., h2/3,             "{:.1f}".format(percent_h2), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)
                plt.text(rip.get_x() + rip.get_width() / 2., h2 + hip/10,      "{:.1f}".format(percent_ip), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)
                plt.text(rh1.get_x() + rh1.get_width() / 2., h2 + hip + h1/2., "{:.1f}".format(percent_h1), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)
            elif ell == 1:
                #plt.text(rh2.get_x() + rh2.get_width() / 2., h2/1.5,            r"$\mathbb{{{:.1f}}}$".format(percent_h2), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                #plt.text(rip.get_x() + rip.get_width() / 2., h2 + hip/4,      r"$\mathbb{{{:.1f}}}$".format(percent_ip), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                #plt.text(rh1.get_x() + rh1.get_width() / 2., h2 + hip + h1/2., r"$\mathbb{{{:.1f}}}$".format(percent_h1), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=10)
                plt.text(rh2.get_x() + rh2.get_width() / 2., h2/1.5,            "{:.1f}".format(percent_h2), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)
                plt.text(rip.get_x() + rip.get_width() / 2., h2 + hip/4,        "{:.1f}".format(percent_ip), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)
                plt.text(rh1.get_x() + rh1.get_width() / 2., h2 + hip + h1/2.,  "{:.1f}".format(percent_h1), color=color_str, ha="center", va="bottom", fontweight="bold", fontsize=font_s)

        # legend outside
        kw = dict(ncol=3, loc="lower center", frameon=False)    
        plt.legend((p_hoist1[0], p_ip[0], p_hoist2[0]), (r'$h \cdot C_{H1}$', r"$(g(h)-h)\cdot C_{IP}$", r"$h \cdot C_{H2}$"), bbox_to_anchor=[0.5,1.00],**kw, fontsize=14)

        savename = "{}/{}_k{}_lgN{}_M{}_ell{}".format(filedir, prefixname, k, logN, M, ell)
        if is_png:
            plt.savefig("{}.png".format(savename))
        else:
            plt.savefig("{}.eps".format(savename))
        plt.clf()

if __name__ == "__main__" :
    # max
    argc = len(sys.argv)
    if argc != 9:
        print("python3 ./xxx [logN] [ell] [k] [KS:0(lattigo) 1(PALISADE)] [IsLastDigit: 0 (NO) 1 (YES)] [PrefixName] [ISDEBUG] [ISPNG]")
        exit(1)

    logN               = int(sys.argv[1])
    ell                = int(sys.argv[2])
    k                  = int(sys.argv[3])
    ks_id              = int(sys.argv[4])
    is_last_digit = bool(int(sys.argv[5]))
    prefixname         =     sys.argv[6]
    is_debug           = int(sys.argv[7])
    is_png             = int(sys.argv[8])
    num_rns_ctxt =  ell + 1

    is_coeff   = 0
    if is_coeff:
        print("Coefficient Rep")
    else:
        print("Eval Rep")

    if logN < 10 or logN > 16:
        print("logN must be 13 14 15")
        exit(1)

    print("NOTE: This script assumes Eval-Eval (Input-Output) nalysis")



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

    Ms = {}
    Ms[logN] = [2, logN >> 1, logN]
    for M in Ms[logN]:
        hs = list(range(1, M + 1))
        ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, ell, k, M, hs, m_ks, eta, is_last_digit)
        res[M] = ca.cost()
        print("M = {}, opth={}, max_speedup = {}".format(M, res[M].optimal_h, res[M].max_speedup))
    #ca.print()
    for M in Ms[logN]:
        print("M = {}, costs = {}".format(M, res[M].costs))

    plot_cost(logN, M, Ms, res, k, ell, prefixname)
    plot_cost_breakdown(logN, M, Ms, res, k, ell, prefixname)

    for M in Ms[logN]:
        result = res[M]
        n = len(result.numkeys_blowup)
        print("M={} Speed::Mem ".format(M) )
        for i in range(0, n):
            if result.speedups[i] > result.numkeys_blowup[i]:
                print("SPEED h={} {}::{}".format(i+1, result.speedups[i], result.numkeys_blowup[i]))
            else:
                print("MEMOR h={} {}::{}".format(i+1, result.speedups[i], result.numkeys_blowup[i]))
        print("Rate= {}".format(result.unit_ch_div_cip))
