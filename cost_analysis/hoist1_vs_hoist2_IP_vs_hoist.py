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
is_debug = False
is_png = False
eta = 0

if __name__ == "__main__" :
    # max
    argc = len(sys.argv)
    if argc != 8:
        print("python3 ./xxx [logN] [ell]  [KS:0(lattigo) 1(PALISADE)] [IsLastDigit:0 (NO) 1 (YES)] [PrefixName] [ISDEBUG] [ISPNG]")
        exit(1)

    logN          =  int(sys.argv[1])
    ell           =  int(sys.argv[2])
    ks_id         =  int(sys.argv[3])
    is_last_digit = bool(int(sys.argv[4]))
    prefixname    =      sys.argv[5]
    is_debug      =  int(sys.argv[6])
    is_png        =  int(sys.argv[7])
    L = ell

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
    M =  15
    hs = list(range(1, M + 1))
    ks = []
    for k in range(1, L+2):
        if (L + 1) % k ==0:
            ks.append(k)
    print("Fixing k anc see variable ell? Yes 1: No (Fixed=ell, var k): 0")
    is_k_fix = int(input())
    print("M={}",format(M))
    if is_k_fix > 0:
        ells = list(range(0,ell+1))
        for k in ks:
            print("****k ={}".format(k))
            rate = []
            ell_h1s = []
            ell_h2s = []
            cms = []
            for l in ells:
                beta = math.ceil((l + 1)/k)
                ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, k, M, hs, m_ks, eta, is_last_digit)
                res = ca.cost()
                h1 =  res.unit_ch1
                h2 =  res.unit_ch2
                ch =  h1 + h2
                cip=  res.unit_cip
                cm = res.costs[-1]
                cms.append(cm)
        #        if h1 > h2:
        #            print("\tH1BIG H1 = {}, H2 ={}".format(h1, h2))
        #        else:
        #            print("\tH2BIG H1 = {}, H2 ={}".format(h1, h2))
        #        print("\tCH! / CH2 = {}".format(h1/h2))
                winner = "H1"
                if h2 < h1:
                    winner = "H2"

                print("(ell, beta) = ({}, {}) , opth={}, max_speedup = {:.2f}, Win(H1 vs H2) = {}, C(M)={} CH/CIP={:.3f}".format(l,beta, res.optimal_h, res.max_speedup, winner, math.log(res.costs[-1], 2) ,ch/cip))
            plt.plot(ells, cms, linestyle="-", label="k={}".format(k))
            plt.xlabel("ell")
            plt.ylabel("cost")
            plt.xticks(np.arange(0, max(ells)+1, 1), fontsize=5)
            plt.yscale("log")
            plt.legend()
            dirnam="/tmp/result/no_last_digit_costM/"
            if is_last_digit:
                dirnam="/tmp/result/last_digit_costM/"
            os.system("mkdir -p {}".format(dirnam))
        plt.savefig("{}{}costM_logN{}_L{}.eps".format(dirnam, prefixname, logN, L))
        plt.clf()

# RATE
        for k in ks:
            print()
            rate = []
            for l in ells:
                beta = math.ceil((l + 1)/k)
                ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, k, M, hs, m_ks, eta, is_last_digit)
                res = ca.cost()
                h1 =  res.unit_ch1
                h2 =  res.unit_ch2
                ch =  h1 + h2
                cip=  res.unit_cip
                rate.append(ch/cip)
            plt.plot(ells, rate, linestyle="-", marker="o", label="k={}".format(k))
        plt.xlabel(r"$\ell$", fontsize=22)
        plt.ylabel(r"Ratio ($C_H/C_{IP}$)", fontsize=22)
        plt.xticks(np.arange(0, max(ells)+1, 1))
        kw = dict(ncol=4, loc="lower center", frameon=False)
        plt.legend(bbox_to_anchor=[0.5,1.00],**kw, fontsize=14)
        dirnam="/tmp/result/no_last_digit_rate/"
        if is_last_digit:
            dirnam="/tmp/result/last_digit_rate/"
        os.system("mkdir -p {}".format(dirnam))
        plt.savefig("{}{}ch_to_cip_logN{}_L{}.eps".format(dirnam, prefixname, logN, L))
        plt.clf()

# Seperate H1 and IP
        for k in ks:
            print("k={}".format(k))
            h1_plus_h2s = []

            l = ells[0]
            ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, k, M, hs, m_ks, eta, is_last_digit)
            res = ca.cost()
            prev_h1 =  res.unit_ch1
            prev_h2 =  res.unit_ch2
            prev_ch =  prev_h1 + prev_h2
            prev_cip=  res.unit_cip

            
            beta = math.ceil((l + 1)/k)
            cips = [prev_cip]
            ch1s = [prev_h1]
            ch2s = [prev_h2]
            chs  = [prev_ch]
            betas = [beta]
            for l in ells[1:]:
                beta = math.ceil((l + 1)/k)
                ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, k, M, hs, m_ks, eta, is_last_digit)
                res = ca.cost()
                h1 =  res.unit_ch1
                h2 =  res.unit_ch2
                ch =  h1 + h2
                cip=  res.unit_cip

                ch_increase = ch / prev_ch
                cip_increase = cip / prev_cip

                prev_h1 = h1
                prev_h2 = h2
                prev_cip = cip
                cips.append(cip)
                ch1s.append(h1)
                ch2s.append(h2)
                chs.append(ch)
                betas.append(beta)

            plt.plot(ells, cips, linestyle="-", marker="o", label="k={}CIP".format(k))
            #plt.plot(ells, chs, linestyle="-", marker="o", label="k={}CH".format(k))
            #plt.plot(ells, ch2s, linestyle="-", marker="o", label="k={}CH1".format(k))
            #plt.plot(ells, ch2s, linestyle="-", marker="o", label="k={}CH2".format(k))
            for i in range(1, len(chs)):
                ch_rate  = chs[i]/chs[i-1]
                cip_rate = cips[i]/cips[i-1]
                beta_str = "Same"
                if betas[i] > betas[i-1]:
                    beta_str = "UP"
                print(" l-1 ({}) -> l({}) b_l-1 ({})-> b ({}) _lCH1 slope ={:.2f}, CH2 slope ={:.2f}, CH slope = {:.2f}, CIPSlope = {:.9f}, beta: {}\t CH > CIP: {}".format(i-1, i, betas[i-1], betas[i], ch1s[i]/ch1s[i-1], ch2s[i]/ch2s[i-1], ch_rate, cip_rate, beta_str, ch_rate > cip_rate))
        plt.xlabel(r"$\ell$", fontsize=22)
        plt.ylabel("costs", fontsize=22)
        plt.yscale("log")
        plt.xticks(np.arange(0, max(ells)+1, 1))
        kw = dict(ncol=4, loc="lower center", frameon=False)
        plt.legend(bbox_to_anchor=[0.5,1.00],**kw, fontsize=14)
        dirnam = "/tmp/result/no_last_digit_up/"
        if is_last_digit:
            dirnam = "/tmp/result/last_digit_up/"

        os.system("mkdir -p {}".format(dirnam))
        plt.savefig("{}{}ch_to_cip_logN{}_L{}.eps".format(dirnam, prefixname, logN, L))
        plt.clf()


        #    plt.plot(ells, rate, label="k={}".format(k), marker="o")
        #plt.xlabel("ell")
        #plt.ylabel("Ratio")
        #plt.xticks(np.arange(0, max(ells)+1, 1), fontsize=5)
        #plt.legend()
        #plt.savefig("/tmp/result/rate/logN{}_L{}.eps".format(logN, L))
    else:
        for l in range(1, ell+1):
            print("*** ell = {}".format(l))
            k = ks[0]
            beta = math.ceil((l + 1)/k)
            ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, ks[0], M, hs, m_ks, eta, is_last_digit)
            res = ca.cost()
            h1 =  res.unit_ch1
            h2 =  res.unit_ch2
            ch =  h1 + h2
            cip=  res.unit_cip
            winner = "H1"
            if h2 < h1:
                winner = "H2"
            print("M = {} (k, beta) = ({},{}) , opth={}, max_speedup = {:.2f}, Win(H1 vs H2) = {} C(M)={:.2f} CH/CIP={:.2f}".format(M, k, beta, res.optimal_h, res.max_speedup, winner, math.log(res.costs[-1],2), ch/cip))

            prev_res = res
            for k in ks[1:]:
                beta = math.ceil((l + 1)/k)
                ca = CostHelper.CostAnalysisAdvanced(is_coeff, logN, l, k, M, hs, m_ks, eta, is_last_digit)
                res = ca.cost()
                h1 =  res.unit_ch1
                h2 =  res.unit_ch2
                ch =  h1 + h2
                cip=  res.unit_cip
        #        if h1 > h2:
        #            print("\tH1BIG H1 = {}, H2 ={}".format(h1, h2))
        #        else:
        #            print("\tH2BIG H1 = {}, H2 ={}".format(h1, h2))
        #        print("\tCH! / CH2 = {}".format(h1/h2))
                winner = "H1"
                if h2 < h1:
                    winner = "H2"
                updown = "down"
                if prev_res.costs[-1] < res.costs[-1]:
                    updown = "up"
                print("M = {} (k, beta) = ({},{}) , opth={}, max_speedup = {:.2f}, Win(H1 vs H2) = {} C(M)={:.2f} CH/CIP={:.2f}    C(M)...\t{}".format(M, k, beta, res.optimal_h, res.max_speedup, winner, math.log(res.costs[-1],2), ch/cip, updown))
                prev_res = res
