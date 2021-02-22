import math
import sys

logN = 14
N =  1 << logN

def f(a, b, l):
    return -2*(a**2)*(b**2) -6*(a**2)*b + 4*(a**2)*l + 2*a*(b**2) + 4*a*(l**2) + 8*(a**2) + 22*a + 18*l + 12 

def g(a, b, l):
    return -2*a*(b**2) + 4*(a**2) - 8*a*b + 14*a*l + 10*(l**2) + 18*a + 30*l + 20

def gg(a, b, l):
    return 4*a*l - 2*a*(b**2) - 6*a*b + 2*(b**2) + 4 *(l**2) + 8*a + 22

if __name__ == "__main__" :
   argc = len(sys.argv)
   if argc != 2:
       print("p ./xxxx ell")
       exit(1)
   ell = int(sys.argv[1])
   ks = []
   for k in range(1, ell + 2):
       if (ell+1) % k == 0:
           ks.append(k)
   #print(ks)
   ells =  range(0, ell+1)
   for k in ks:
       for ell in ells:
           beta = math.ceil((ell+1)/k)
           if (f(k,beta, ell)) <0:
               print("Negative")
           if (g(k,beta, ell)) <0:
               print("Negative")
           if (gg(k,beta, ell)) <0:
               print("Negative")
