import math, sys
import numpy as np
import matplotlib.pyplot as plt
import collections
import math
import timeit

import HKcost as HK
import Util as u

from cost_analysis_base import Results


# Accept the following 2 cases
#  1. ell < L for a fixed alpha an
#  2. k does not divide ell + 1
class CostAnalysisAdvanced:
    def __init__ (self, is_coeff:bool, logN: int, ell:int, k:int,  M:int, hs:list, m_ks: HK.KSMethod, eta = 0, is_last_digit = False):
        self.is_coeff = is_coeff
        self.logN = logN
        self.ell  = ell
        self.k    = k
        self.d    = math.ceil((ell+1)/k)
        self.M    = M
        self.hs   = hs
        self.m_ks  = m_ks
        self.eta  = eta
        self.is_last_digit = is_last_digit
        self.res  = Results()

    def cost(self):
        [self.res.unit_ch1, self.res.unit_cip, self.res.unit_ch2] = HK.compute_HKcost_ee(self.logN, self.ell, self.k, self.d, self.m_ks)
        if self.is_last_digit == False:
            [self.res.unit_ch1, self.res.unit_cip, self.res.unit_ch2] = HK.compute_HKcost_ee_advanced(self.logN, self.ell, self.k, self.d, self.m_ks)
            
        self.res.unit_ch_div_cip = (self.res.unit_ch1 + self.res.unit_ch2)/self.res.unit_cip
        for h in self.hs:
            # 1. Base number of ctxts (in bits) for each unroll
            logN_div_h_floor = int(self.M / h)
            #  Base number of ctxts  for each unroll
            num_ctxt_base = 1 << logN_div_h_floor
            # 2. The number of Ceil = Floor + 1unrolling
            rem = self.M % h
            # Total Number of key-switchings in all unrolled iteration
            # Since each unrolled iteration consists of exactly no-keyswitching, we subtract it (i.e., h in total).
            #num_ctxt_rem = 1 << (logN_div_h_floor + 1)
            #total_num_ks = (num_ctxt_base - 1) * (h - rem)  + (num_ctxt_rem - 1) * rem
            total_num_ks = (num_ctxt_base * (h + rem)) - h
            self.res.num_ips_unroll.append(total_num_ks)
    
            # Compute unrolled cost
            hoist1_cost = 0
            ip_cost     = 0
            hoist2_cost = 0
            if self.is_coeff:
                raise Exception("Costs for advanced cc is not defined.")
            else:
                [hoist1_cost, ip_cost, hoist2_cost] = HK.unrolled_ee_cost_advanced(h, total_num_ks, self.logN, self.ell, self.k, self.d,  self.m_ks)
            if hoist1_cost == 0 or ip_cost == 0 or hoist2_cost == 0:
                raise Exception("Costs should not be zero.")
    
            total_cost  = hoist1_cost + ip_cost + hoist2_cost
    
            self.res.total_c_hoist1.append(hoist1_cost)
            self.res.total_c_ips   .append(ip_cost)
            self.res.total_c_hoist2.append(hoist2_cost)
            self.res.total_c_hoist .append(hoist1_cost + hoist2_cost)
    
            self.res.percentage_hoist1.append(100 * hoist1_cost / total_cost)
            self.res.percentage_ip    .append(100 * ip_cost / total_cost)
            self.res.percentage_hoist2.append(100 * hoist2_cost / total_cost)
            self.res.percentage_all_hoist.append(100 * (hoist1_cost + hoist2_cost) / total_cost)
            self.res.costs.append(total_cost)

        self.res.optimal_h = np.argmin(self.res.costs) + 1
        self.res.speedups  = [self.res.costs[-1]/x for x in self.res.costs]
        self.res.max_speedup  = self.res.costs[-1]/min(self.res.costs)
        self.res.numkeys_blowup  = [x/self.res.num_ips_unroll[-1] for x in self.res.num_ips_unroll]
        return self.res
    
    def cost_w_automorph(self):
        [self.res.unit_ch1, self.res.unit_cip, self.res.unit_ch2] = HK.compute_HKcost_ee(self.logN, self.ell, self.k, self.d, self.m_ks)
        if self.is_last_digit == False:
            [self.res.unit_ch1, self.res.unit_cip, self.res.unit_ch2] = HK.compute_HKcost_ee_advanced(self.logN, self.ell, self.k, self.d, self.m_ks)

        self.res.unit_ch_div_cip = (self.res.unit_ch1 + self.res.unit_ch2)/self.res.unit_cip
        for h in self.hs:
            # 1. Base number of ctxts (in bits) for each unroll
            logN_div_h_floor = int(self.M / h)
            #  Base number of ctxts  for each unroll
            num_ctxt_base = 1 << logN_div_h_floor
            # 2. The number of Ceil = Floor + 1unrolling
            rem = self.M % h
            # Total Number of key-switchings in all unrolled iteration
            # Since each unrolled iteration consists of exactly no-keyswitching, we subtract it (i.e., h in total).
            #num_ctxt_rem = 1 << (logN_div_h_floor + 1)
            #total_num_ks = (num_ctxt_base - 1) * (h - rem)  + (num_ctxt_rem - 1) * rem
            total_num_ks = (num_ctxt_base * (h + rem)) - h
            self.res.num_ips_unroll.append(total_num_ks)
    
            # Compute unrolled cost
            hoist1_cost = 0
            ip_cost     = 0
            hoist2_cost = 0
            if self.is_coeff:
                raise Exception("Costs for advanced cc is not defined.")
            else:
                [hoist1_cost, ip_cost, hoist2_cost] = HK.unrolled_ee_cost_advanced(h, total_num_ks, self.logN, self.ell, self.k,  self.m_ks)
            if hoist1_cost == 0 or ip_cost == 0 or hoist2_cost == 0:
                raise Exception("Costs should not be zero.")
    

            total_cost  = hoist1_cost + ip_cost + hoist2_cost

            cauto = self.eta * (1 << self.logN)
            # No hoist
            if h == self.hs[-1]:
                total_cost  += h * cauto *(2 * self.ell + 2) # only on R_{Q}^2
            # Hoist
            else:
                total_cost  += total_num_ks * cauto *(self.ell + 1 + 2*(self.ell + self.k + 1))
    
            self.res.total_c_hoist1.append(hoist1_cost)
            self.res.total_c_ips   .append(ip_cost)
            self.res.total_c_hoist2.append(hoist2_cost)
            self.res.total_c_hoist .append(hoist1_cost + hoist2_cost)
    
            self.res.percentage_hoist1.append(100 * hoist1_cost / total_cost)
            self.res.percentage_ip    .append(100 * ip_cost / total_cost)
            self.res.percentage_hoist2.append(100 * hoist2_cost / total_cost)
            self.res.percentage_all_hoist.append(100 * (hoist1_cost + hoist2_cost) / total_cost)
            self.res.costs.append(total_cost)

        self.res.optimal_h = np.argmin(self.res.costs) + 1
        self.res.speedups  = [self.res.costs[-1]/x for x in self.res.costs]
        self.res.max_speedup  = self.res.costs[-1]/min(self.res.costs)
        return self.res

    ##### For experiment each among many h's
#    optimial_hs = []
#    optimial_costs_unroll = []
#
#    optimial_costs_unroll.append(min(costs))
    def print(self):
        print("\t Unroll speedups:\n\t[",end="")
        for x in self.res.costs:
           print(" {:.2f} ".format(self.res.costs[-1]/x), end ="")
        print("]")
           
        print("\t Storage Blowups:\n\t [",end="")
        for x in self.res.num_ips_unroll:
           print(" {:.2f} ".format(x/self.res.num_ips_unroll[-1]), end ="")
        print("]")

        print("\t Unroll Max-speedup:\t {}".format(self.res.costs[-1] / min(self.res.costs)))
        print("\n\n\t Costs {}".format(self.res.costs))
        print("        hs = {}".format(list(self.hs)))
        print("      #IPS = {}".format(self.res.num_ips_unroll))
        #print("  OptCosts = {}".format(res.optimial_costs_unroll))
        #print("  Opt    H = {}".format(res.optimial_hs))

        print("BreakDown Percentage:\n")
        print(" Hoist1: ",end="")
        for x in self.res.percentage_hoist1:
           print(" {:.2f} \t".format(x), end ="")
        print()

        print("     IP: ",end="")
        for x in self.res.percentage_ip:
           print(" {:.2f} \t".format(x), end ="")
        print()

        print(" Hoist2: ",end="")
        for x in self.res.percentage_hoist2:
           print(" {:.2f} \t".format(x), end ="")
        print()

        print(" All Hoist: ",end="")
        for x in self.res.percentage_all_hoist:
           print(" {:.2f} \t".format(x), end ="")
        print()


        print("Magnitude Analysis Alpha terms vs Beta terms")
        print(" Hoist1 + IP (Beta Term): ",end="")
        for i in range(0, len(self.res.total_c_hoist1)):
           print(" {} \t".format(self.res.total_c_hoist1[i] + self.res.total_c_ips[i]), end ="")
        print()

        print(" Hoist2                 : ",end="")
        for x in self.res.total_c_hoist2:
           print(" {} \t".format(x), end ="")
        print()

        print(" Which one is smaller?                 : ",end="")
        for i in range(0, len(self.res.total_c_hoist1)):
           if self.res.total_c_hoist2[i]  > (self.res.total_c_hoist1[i] + self.res.total_c_ips[i]):
               print(" Beta \t", end ="")
           elif self.res.total_c_hoist2[i]  < (self.res.total_c_hoist1[i] + self.res.total_c_ips[i]):
               print(" Alph \t", end ="")
        print()
