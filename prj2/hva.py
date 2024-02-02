###########################################################
#             Harmonic Vibrational Analysis               #
###########################################################

import numpy as np
from math import *
geometry = open('geom.dat','r')
lines = geometry.readlines()
atnum = []
x = []
y = []
z = []
for i in range(len(lines)):
    toks = lines[i].split()
    if i == 0:
       total_atoms = int(toks[0])
    else:
       atnum.append(int(float(toks[0])))
       x.append(float(toks[1]))
       y.append(float(toks[2]))
       z.append(float(toks[3]))
#
# ******* Reading Hessian matrix *******
#
hessian = open('hess.dat','r')
total_atoms = int(hessian.readline())
lines = hessian.readlines()
row = col = total_atoms*3
F = np.zeros((row,col))
m = 0
for i in range(row):
    k = 0
    for j in range(m,m+3):
        toks = lines[j].split()
        F[i,3*k] = toks[0]
        F[i,3*k+1] = toks[1]
        F[i,3*k+2] = toks[2]
        k += 1
    m += 3
print(F)

#
# ******* Mass-weighted Hessian matrix *******
#
mass = {1: 1.007825, 2: 4.002603, 3: 6.015122, 4: 9.012183,
        5: 10.012936, 6: 12.000000, 7: 14.003074, 8: 15.994914}
F_m = np.zeros((row,col))
def mat_conversion(li,ui,lj,uj,a,b):
    for i in range(row):
        for j in range(col):
            if li <= i < ui and lj <= j < uj:
                F_m[i,j] = F[i,j]/(sqrt(mass[a]*(mass[b])))
    return F_m[i,j]
mat_conversion(0,3,0,3,atnum[0],atnum[0])
mat_conversion(0,3,3,6,atnum[0],atnum[1])
mat_conversion(0,3,6,9,atnum[0],atnum[2])
mat_conversion(3,6,0,3,atnum[1],atnum[0])
mat_conversion(3,6,3,6,atnum[1],atnum[1])
mat_conversion(3,6,6,9,atnum[1],atnum[2])
mat_conversion(6,9,0,3,atnum[2],atnum[0])
mat_conversion(6,9,3,6,atnum[2],atnum[1])
mat_conversion(6,9,6,9,atnum[2],atnum[2])

print(F_m)
#
# ******* Diagonalization of F_m *******
#
w,l = np.linalg.eig(F_m)
w = sorted(w)
print (' Hessian eigenvalues (hartree/amu-bohr^2):')
for i in range(total_atoms*3):
    print ('%12.10f'%w[i])
#
# ******* Harmonic vibrational frequencies *******
#
conv = (4.359744650e-18)/(((5.29e-11)**2)*(1.660539066e-27)) # energy conversion from au
c = 2.997924e10                                              # to cm^-1 
omega = ['']*(total_atoms*3)
print (' Harmonic vibrational frequencies (cm^-1):')
print(pi)
for i in range(9):
    omega[i] = sqrt((abs(w[i])*conv)/(4*((pi)**2)*(c)**2))
    print ('%9.4f'%omega[i])

# ******************** X ***********************






