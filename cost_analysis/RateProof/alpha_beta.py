import math
import numpy as np

ks = [1, 2, 5, 10, 20]
ells = np.arange(1, 39, 1)

for ell in ells:
    #print("ell = {}".format(ell))
    prods = []
    keys=[]
    for k in ks:
        beta = math.ceil((ell+1)/k)
        prod = k*beta
        prods.append(prod)
        keys.append("k = {}, beta={}".format(k, beta))
        #print(" k * beta = {} * {} = {}".format(k,beta, prod))
    i = np.argmax(prods)
    print("ell = {} -> {}".format(ell, keys[i]))
