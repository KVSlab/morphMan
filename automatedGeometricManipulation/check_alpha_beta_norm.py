from os import path, listdir
import numpy as np
pathvals = "plots_angle/optimal_alphabeta.txt"

info = np.loadtxt(pathvals, dtype="str")
print info
norms = []
for i in range(len(info)-20,len(info)):
    alpha = float(info[i][2].split("=")[1])
    beta = float(info[i][3].split("=")[1])
    norms.append(np.sqrt(alpha**2 + beta**2))

print norms
print np.mean(norms)
