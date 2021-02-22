import os

logN=15

# Lattigo
ksid = 0
os.system("python3 breakdown_analysis.py {} 9 1  {} 0 lattigo 0 0".format(logN, ksid))
os.system("python3 breakdown_analysis.py {} 9 10 {} 0 lattigo 0 0".format(logN, ksid))
os.system("python3 breakdown_analysis.py {} 1 1  {} 0 lattigo 0 0".format(logN, ksid))
os.system("python3 breakdown_analysis.py {} 1 2  {} 0 lattigo 0 0".format(logN, ksid))

# PALISADE
ksid = 1
islastdigit=0
os.system("python3 breakdown_analysis.py {} 9 1  {} {} pali 0 0".format(logN, ksid, islastdigit))
os.system("python3 breakdown_analysis.py {} 9 10 {} {} pali 0 0".format(logN, ksid, islastdigit))
os.system("python3 breakdown_analysis.py {} 1 1  {} {} pali 0 0".format(logN, ksid, islastdigit))
os.system("python3 breakdown_analysis.py {} 1 2  {} {} pali 0 0".format(logN, ksid, islastdigit))
