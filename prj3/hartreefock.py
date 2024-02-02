######################################################
#              Hartree-Fock SCF code                 #
######################################################
from math import *
import numpy as np




N = 10  # total number of electrons
result = open('result.dat','w')
#
# ******* Reading nuclear repulsion energy ********
#
nu_rep = open('enuc.dat','r')
nu_rep = float(nu_rep.read())
result.write('Nuclear repulsion energy = ' + str(nu_rep) + '\n')
#
#Print Matrix
#
def print_mat(mat,output_statement):
    nrow = len(mat)
    result.write('\n' + str(output_statement) + '\n')
    for i in range(nrow):
        line = ''
        for j in range(nrow):
            line = line + '%11.7f'%mat[i,j]+ '  '
        result.write(str(line) + '\n')
    return mat        
#
# ******* Reading precomputed one-electron integrals ********
#
def matrix(filename,m):
    values = open(filename,'r')
    lines = values.readlines()
    for i in range(len(lines)):
        toks = lines[i].split()
        m[int(toks[0])-1,int(toks[1])-1] = float(toks[2])
    return
S = np.zeros((7,7))
T = np.zeros((7,7))
V = np.zeros((7,7))
matrix('overlap.dat',S)
matrix('kin.dat',T)
matrix('nuc_att.dat',V)
for i in range(7):
    for j in range(7):
        S[i,j] = S[j,i]
        T[i,j] = T[j,i]
        V[i,j] = V[j,i]
print_mat(S,'Overlap Integrals:')
print_mat(T,'Kinetic-Energy Integrals:')
print_mat(V,'Nuclear Attraction Integrals:')
#
# ******** Core Hamiltonian construction ********
#
H = np.zeros((7,7))
for i in range(7):
    for j in range(7):
        H[i,j] =  T[i,j] + V[i,j]
print_mat(H,'Core Hamiltonian:')
#
# ******* Reading precomputed two-electron integrals ********
#
TEI = np.zeros((7,7,7,7))
eri = open('eri.dat','r')
lines = eri.readlines()
for i in range(len(lines)):
    toks = lines[i].split()
    TEI[int(toks[0])-1,int(toks[1])-1,int(toks[2])-1,int(toks[3])-1] = float(toks[4])
    p = int(toks[0])-1
    q = int(toks[1])-1
    r = int(toks[2])-1
    s = int(toks[3])-1
    TEI[q,p,r,s] = float(toks[4])
    TEI[p,q,s,r] = float(toks[4])
    TEI[q,p,s,r] = float(toks[4])
    TEI[r,s,p,q] = float(toks[4])
    TEI[s,r,p,q] = float(toks[4])
    TEI[r,s,q,p] = float(toks[4])
    TEI[s,r,q,p] = float(toks[4])
#
# ******* Symmetric Orthogonalization *******
#
w,L = np.linalg.eig(S)
lamb = np.diag(w)
lamb_sym = np.zeros((7,7))
for i in range(7):
    for j in range(7):
        if i == j:
            lamb_sym[i,j] = (lamb[i,j])**(-0.5) # matrix /\^(-1/2)
S_sym = np.zeros((7,7))
S_sym = np.dot(L,np.dot(lamb_sym,np.transpose(L)))
print_mat(S_sym,'S^-1/2 Matrix:')
#
# ******* Initial guess density *******	
#
F0 = np.dot(np.transpose(S_sym),np.dot(H,S_sym)) # initial Fock matrix
F0_evals,C0prime = np.linalg.eig(F0)
idx = np.argsort(F0_evals)
F0_evals = F0_evals[idx]
C0prime = C0prime[:,idx]
C0 = np.zeros((7,7))
C0 = np.dot(S_sym,C0prime)    # initial coefficient matrix
D = np.zeros((7,7))
for i in range(7):
    for j in range(7):
        for a in range(int(N/2)):
            D[i,j] += C0[i,a]*C0[j,a] # initial density matrix
print_mat(F0,'Initial F Matrix:')
print_mat(C0,'Initial C matrix:')
print_mat(D,'Initial density matrix:')
#
# ******* Initial SCF energy *******
#
E_ele = 0.0
for i in range(7):
    for j in range(7):
        E_ele += D[i,j]*(H[i,j]+H[i,j])
E_total = E_ele + nu_rep 
result.write('\nIter = 0' + ', E(total) = ' + str(E_total))
E = []
D_list = []
E.append(E_total)
D_list.append(D)
#
# ******* New Fock matrix *******
#
delta_E = 1.0
rms_D = 1.0
iteration = 1 

# iteration starts here ---->

while delta_E > 10**(-12) or rms_D > 10**(-12) and iteration > 0:
      G = np.zeros((7,7))
      F = np.zeros((7,7))
      for i in range(7):
          for j in range(7):
              for k in range(7): 
                  for l in range(7):
                      G[i,j] += D[k,l]*(2*TEI[i,j,k,l] - TEI[i,k,j,l])
      F = H + G   # new Fock matrix
      print_mat(F,'New Fock matrix:')
#
# ******* New density *******
#
      Fprime = np.dot(np.transpose(S_sym),np.dot(F,S_sym)) # orthogonalization of new Fock matrix
      Fprime_evals,Cprime = np.linalg.eig(Fprime)
      idx = np.argsort(Fprime_evals)
      Fprime_evals = Fprime_evals[idx]
      Cprime = Cprime[:,idx]
      C = np.dot(S_sym,Cprime)    # new coefficient matrix
      D = np.zeros((7,7))
      for i in range(7):
          for j in range(7):
              for a in range(int(N/2)):
                  D[i,j] += C[i,a]*C[j,a] # new density matrix
      D_list.append(D)
#
# ******* New SCF energy *******
#
      E_ele = 0.0
      for i in range(7):
          for j in range(7):
              E_ele += D[i,j]*(H[i,j]+F[i,j])
      E_total = E_ele + nu_rep
      E.append(E_total)
      delta_E = abs(E[iteration] - E[iteration-1])
      rms_D = 0.0
      for i in range(7):
          for j in range(7):
              rms_D += (D_list[iteration][i,j]-D_list[iteration-1][i,j])**2
      rms_D = sqrt(rms_D)
      result.write('\nIter = ' + str(iteration) + ', E(total) = ' + str(E_total) + \
      ', Delta(E) = ' + str(delta_E) + ', RMS(D) = ' + str(rms_D))
      iteration += 1

# **********************  X  ************************





