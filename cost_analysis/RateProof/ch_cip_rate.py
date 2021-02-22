import math
import sys

logN = 15
N =  1 << logN

def cip(ell, beta, k):
    return 2*beta *(ell + 1 + k) 

# Original formula
def ch(ell, beta,  k):
    nlogn = logN *(ell * (beta + 5) + beta + 2*k + 3)
    n     = ell**2 + ell * (6 + 2*k) + 4*k + 5
    return nlogn + n

#def ch(ell, beta,  k):
#    nlogn = (ell * (beta + 5) + beta + 2*k + 3)
#    n     = 0
#    return nlogn + n

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

   for a in loop1:
       print()
       for b in loop2:
           l = a
           k = b
           if is_k:
               l = b
               k = a

           beta_old = math.ceil(l/k)
           beta_new = math.ceil((l+1)/k)
           BetaStr="Stay"
           if beta_old != beta_new:
               BetaStr="Up"

           CHupdown ,  ch_old,  ch_new  =  ch_grow(l, k)
           CIPupdown, cip_old, cip_new  = cip_grow(l, k)
           ch_change  =  ch_new / ch_old
           cip_change = cip_new / cip_old
           ch_plus_cip_new = ch_new + cip_new

           ch_diff  = ch_new   - ch_old
           cip_diff = cip_new  - cip_old

           ch_greater_than_cip_est = cip_old * ch_diff > ch_old * cip_diff
           #ch_greater_than_cip_est = cip_old > ch_old

           chcip="CH/CIP DOWN"
           if ch_change > cip_change:
               chcip="CH/CIP UP"
           
           # l -1 -> l
           a_str = "l"
           b_str = "k"
           if is_k:
               a_str = "k"
               b_str = "l"

           print("  {}, {}, b=({}, {}, {})\t Beta:{}\t CH={}, \tCIP={}, \tCH +CIP={}\t CCH:{},\tCCIP:{}\t(CCH, \tCCIP) = ({:.2f}, {:.2f}) -> CH/CIP {} est = {}".format(a_str, b_str, a, b, beta_new, BetaStr, ch_new, cip_new, ch_plus_cip_new, CHupdown, CIPupdown, ch_change, cip_change, chcip, ch_greater_than_cip_est))

   print(" NOTE all up down w.r.t ell")
