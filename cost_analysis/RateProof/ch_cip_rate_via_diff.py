import math
import sys

logN = 14
N =  1 << logN

def cip(ell, beta, k):
    return 2*beta *(ell + 1 + k) 

# Original formula
def ch(ell, beta,  k):
    nlogn = logN *(ell * (beta + 5) + beta + 2*k + 3)
    n     = ell**2 + ell * (6 + 2*k) + 4*k + 5
    return nlogn + n

def ch_nlogn(ell, beta,  k):
    nlogn = logN *(ell * (beta + 5) + beta + 2*k + 3)
    return nlogn

def ch_n(ell, beta,  k):
    n     = ell**2 + ell * (6 + 2*k) + 4*k + 5
    return n

# Correct
def ch_diff_beta_same(ell, beta, k):
    return logN*(beta + 5) + (7 + 2*k+2*ell)

def ch_diff_beta_same_nlogn(ell, beta, k):
    return logN*(beta + 5) 

def ch_diff_beta_same_n(ell, beta, k):
    return  (7 + 2*k+2*ell)

def cip_diff_beta_same(ell, beta, k):
    return 2*beta

def ch_diff_beta_up(ell, beta, k):
    return logN * (beta + ell + 7) + 7 + 2*k + 2*ell

def ch_diff_beta_up_nlogn(ell, beta, k):
    return logN * (beta + ell + 7) 

def ch_diff_beta_up_n(ell, beta, k):
    return 7 + 2*k + 2*ell

# 2*b*(l + 1 +k)
# b->b+1: (2b+2)*(l + 2 + k) = 2*b*(l + 1 +k) + 2*b + 2(l + 2+ k)
def cip_diff_beta_up(ell, beta, k):
    return 2*beta + 2*(ell + 2 + k)


def ch_grow(l, k):
    ell = l-1
    beta_old = math.ceil((ell+1) /k)
    old = ch(ell, beta_old, k)

    ell = l
    beta_new =  math.ceil((ell+1) /k)
    new = ch(ell, beta_new, k)

    updown = "Down"
    if new > old:
        updown = "Up"

    return updown, old, new


def cip_grow(l, k):
    ell = l - 1
    beta_old = math.ceil((ell+1) /k)
    old = cip(ell, beta_old, k)

    ell = l
    beta_new = math.ceil((ell+1) /k)
    new = cip(ell,  beta_new, k)

#    print("\tIPNEW = {}".format(new))
    updown = "Down"
    if new > old:
        updown = "Up"
    return updown, old, new
    

if __name__ == "__main__" :
   argc = len(sys.argv)
   if argc != 3:
       print("p ./xxxx ell OutLoopK?(0:ell, 1:k)")
       exit(1)
   ell = int(sys.argv[1])
   is_k = int(sys.argv[2])
   ks = []
   for k in range(1, ell + 2):
       if (ell+1) % k == 0:
           ks.append(k)
   print(ks)
   ells =  range(2, ell+1)

   loop1 = ells
   loop2 = ks
   if is_k:
       loop1 = ks
       loop2 = ells

# NOTE: We are interested in the change from the previous ell to current ell where prev_ell + 1 + cur_ell
   for a in loop1:
       print()
       for b in loop2:
           l = a
           k = b
           if is_k:
               l = b
               k = a

           ell = l - 1
           beta_old = math.ceil((ell+1)/k)

           ell = l
           beta_new = math.ceil((ell+1)/k)

           BetaStr="Stay"
           if beta_old != beta_new:
               BetaStr="Up"


           ell = l - 1
           beta = beta_old
           ch_diff  = ch_diff_beta_same (ell, beta ,k)
           cip_diff = cip_diff_beta_same(ell, beta, k)
# Seperate thins into n term and nlogn term
           ch_diff_nlogn  = ch_diff_beta_same_nlogn (ell, beta ,k)
           ch_diff_n      = ch_diff_beta_same_n     (ell, beta ,k)
           if beta_old != beta_new:
               ch_diff  =  ch_diff_beta_up(ell, beta ,k)
               cip_diff = cip_diff_beta_up(ell, beta, k)
               ch_diff_nlogn  = ch_diff_beta_up_nlogn (ell, beta ,k)
               ch_diff_n      = ch_diff_beta_up_n     (ell, beta ,k)

           cip_old = cip(ell, beta, k)
           ch_old  =  ch(ell, beta, k)

           ch_greater_than_cip_est = cip_old * ch_diff > ch_old * cip_diff
   #        print("\t CIP*CHDIFF, CH*CIPDIFF = {}, {}".format(cip_old *ch_diff, ch_old * cip_diff))

# Seperate thins into n term and nlogn term
           ch_old_nlogn  =  ch_nlogn(ell, beta, k)
           ch_old_n      =  ch_n    (ell, beta, k)

           ch_greater_than_cip_est_nlogn = cip_old * ch_diff_nlogn > ch_old_nlogn * cip_diff
           ch_greater_than_cip_est_n     = cip_old * ch_diff_n     > ch_old_n     * cip_diff

           # l -1 -> l
           a_str = "l"
           b_str = "k"
           if is_k:
               a_str = "k"
               b_str = "l"

           print("  {}, {}, bo, bn=({}, {}, {}->{})\t Beta:{}\tCHCIP: est (total, nlogn, n) = ({}::::::::::\t {}\t {})".\
           format(a_str, \
                  b_str, \
                  a,\
                  b,\
                  beta_old,\
                  beta_new,\
                  BetaStr, \
                  ch_greater_than_cip_est, ch_greater_than_cip_est_nlogn, ch_greater_than_cip_est_n))

   print(" NOTE all up down w.r.t ell")
