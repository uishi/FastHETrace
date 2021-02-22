import math
import sys

base = 2
M = int(sys.argv[1])
print("\n\nM={}".format(M))
divfl = int(M/1)
mod = M % 1
p = math.pow(base, divfl)
c = (p *1)
d = (p *mod)
X = c + d
print("h=1, c ={} d= {} X={}".format( c, d, X))
prevdivfl = divfl
for h in range(2, M+1):
    divfl = int(M/h)
    p = math.pow(base, divfl)
    if divfl != prevdivfl:
        print("=================")
    mod = M % h
    c = (p *h)
    d = (p *mod)
    X = c + d
    print("h={} fl={} 2^fl*h ={} r={} hr= {} X={}".format(h, divfl, c, mod, d, X))
    prevdivfl = divfl
