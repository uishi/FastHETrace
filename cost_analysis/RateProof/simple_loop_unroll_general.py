import math
import numpy as np
import sys

if __name__ == "__main__" :
    argc = len(sys.argv)
    if argc != 2:
        print("python3 ./xxxx")
        exit(1)

    logN = 15
    M = int(sys.argv[1])
    for cip in [1]:
        #for M in range(2, logN+1):
        print("M = {}".format(M))
        for ch in [0.5, 1, 2, 3, 10, 100, 1000, 1000000000000000]:
            print("\t(CIP, CH) = {}, {}".format(cip, ch))
            cs = []
            for h in range(1, M+1):
                logN_div_h_floor = int(M / h)
                #  Base number of ctxts  for each unroll
                num_ctxt_base = 1 << logN_div_h_floor
                # 2. The number of Ceil = Floor + 1unrolling
                rem = M % h
                # Total Number of key-switchings in all unrolled iteration
                # Since each unrolled iteration consists of exactly no-keyswitching, we subtract it (i.e., h in total).                #num_ctxt_rem = 1 << (logN_div_h_floor + 1)
                #total_num_ks = (num_ctxt_base - 1) * (h - rem)  + (num_ctxt_rem - 1) * rem
                total_num_ks = (num_ctxt_base * (h + rem)) - h

                unroll_cost_ip = total_num_ks * cip
                unroll_cost_hoist = h * ch
                cs.append(unroll_cost_ip + unroll_cost_hoist)
            #    if c_ip > c_hoist:
            #        print(" \t c_ip    dominant at h = {} (c)={}".format(h, (cs[-1])))
            #    else:
            #        print(" \t c_hoist dominant at h = {} (c)={}".format(h, (cs[-1])))
            argh = np.argmin(cs) + 1
            conventional = cs[-1]
            is_unroll_all_win = True

            for h in range(0, len(cs)-1): # skip the last
                if cs[h] > cs[-1]:
                    is_unroll_all_win = False
                    break
               
            if argh < M and cs[argh-1] != cs[-1]:
                print("\t\tUnroll well deserved! h = {} Speedup: {}".format(argh, conventional/cs[argh-1]))
                if is_unroll_all_win:
                    print("\t\t\tUnroll ALLWIN!!!!")
            else:
                print("\t\tUnroll does not work")
            #print("\t\tARGH={}\n\n".format(argh))
