from math import *
import numpy as np
indata = open('geom.dat','r')
outdata = open('result.dat','w+')
natoms = indata.readline() 
outdata.write('Number of atoms:'+natoms)
geom = indata.read()
outdata.write('\nInput Cartesian coordinates:\n'+geom)
result = geom.split()
atnum = ['']*7
x = ['']*7
y = ['']*7
z = ['']*7
i = 0
while i <= 6:
	atnum[i] = result[i*4]
	x[i] = result[i*4+1] 
	y[i] = result[i*4+2]
	z[i] = result[i*4+3]
	i += 1 

# *********** Bond length calculation *********

outdata.write('\n'+'Interatomic distances (bohr):')
R = np.zeros((7,7))
for i in range(7):
	for j in range(7):
		if i != j:
			R[i][j] = sqrt(((float(x[i])-float(x[j]))**2)+((float(y[i])-float(y[j]))**2)+((float(z[i])-float(z[j]))**2))
		if i > j:
			outdata.write('\n'+str(i)+' '+str(j)+3*' '+'%.5f'%R[i][j])

# ********* Bond angle calculation ***********

outdata.write(2*'\n'+'Bond angles:')
ex = np.zeros((7,7))
ey = np.zeros((7,7))
ez = np.zeros((7,7))
phi = np.zeros((7,7,7)) 

for i in range(7):
	for j in range(7):
		if i != j:
			ex[i][j] = -(float(x[i])-float(x[j]))/R[i][j]
			ey[i][j] = -(float(y[i])-float(y[j]))/R[i][j]
			ez[i][j] = -(float(z[i])-float(z[j]))/R[i][j]
for i in range(7):
	for j in range(7):
		for k in range(7):
			if i != j and j != k and k != i:
				phi[i][j][k] = degrees(acos(((ex[i][j])*(ex[k][j]))+((ey[i][j])*(ey[k][j]))+((ez[i][j])*(ez[k][j]))))
			if i > j and j > k and R[i][j] < 4.0 and R[k][j] < 4.0:
				outdata.write('\n'+str(i)+'-'+str(j)+'-'+str(k)+3*' '+'%.6f'%(phi[i][j][k]))

# ********** Out-of-Plane angles ************

outdata.write(2*'\n'+'Out-of-plane angles:')
#e_jkl_x = [[[0 for m in range(7)] for n in range(7)] for o in range(7)]
#e_jkl_y = [[[0 for m in range(7)] for n in range(7)] for o in range(7)]
#e_jkl_z = [[[0 for m in range(7)] for n in range(7)] for o in range(7)]
#e_xx = [[[[0 for m in range(7)] for n in range(7)] for o in range(7)] for p in range(7)]
#e_yy = [[[[0 for m in range(7)] for n in range(7)] for o in range(7)] for p in range(7)]
#e_zz = [[[[0 for m in range(7)] for n in range(7)] for o in range(7)] for p in range(7)]
#theta = [[[[0 for m in range(7)] for n in range(7)] for o in range(7)] for p in range(7)]
for i in range(7):
	for j in range(7):
		for k in range(7):
			for l in range(7):
				if i != j and i != k and i != l and  j != k and k != l and l != j and R[k][j] < 4.0 and R[k][l] < 4.0 and R[k][i] < 4.0 and sin(phi[j][k][l]) != 0:
					e_jkl_x = (ey[k][j]*ez[k][l])-(ez[k][j]*ey[k][l])
					e_jkl_y = (ez[k][j]*ex[k][l])-(ex[k][j]*ez[k][l])
					e_jkl_z = (ex[k][j]*ey[k][l])-(ey[k][j]*ex[k][l])
#for i in range(7):
#    for j in range(7):
#        for k in range(7):
#	    for l in range(7):
#		if i != j and i != k and i != l and j != k and k != l and R[k][j] < 4.0 and R[k][l] < 4.0 and R[k][i] < 4.0 and sin(phi[j][k][l]) != 0:
					e_xx = e_jkl_x*ex[k][i]
					e_yy = e_jkl_y*ey[k][i]
					e_zz = e_jkl_z*ez[k][i]
					theta = (e_xx+e_yy+e_zz)/sin(phi[j][k][l])
					if theta < -1.0:
						theta = degrees(asin(-1.0))
					elif theta > 1.0:
						theta = degrees(asin(1.0))
					else:
						theta = degrees(asin(theta))
					print (i,j,k,l,theta)
print ()
# ************* Dihedral Angles ****************

for i in range(7):
    for j in range(7):
        for k in range(7):
            for l in range(7):
	        if i != j and i != k and  i != l and j != k  and j != l and k != l and R[j][i] < 4.0 and R[j][k] < 4.0 and R[k][j] < 4.0 and  R[k][l] < 4.0 and sin(phi[i][j][k]) != 0 and sin(phi[j][k][l]) != 0:
		   e1_x = (ey[j][i]*ez[j][k])-(ez[j][i]*ey[j][k]) 
		   e1_y = (ez[j][i]*ex[j][k])-(ex[j][i]*ez[j][k])
                   e1_z = (ex[j][i]*ey[j][k])-(ey[j][i]*ex[j][k])                  
		   e2_x = (ey[k][j]*ez[k][l])-(ez[k][j]*ey[k][l])
                   e2_y = (ez[k][j]*ex[k][l])-(ex[k][j]*ez[k][l])
                   e2_z = (ex[k][j]*ey[k][l])-(ey[k][j]*ex[k][l])
		   e3 = (e1_x)*(e2_x)+(e1_y)*(e2_y)+(e1_z)*(e2_z)
		   tau = e3/(sin(phi[i][j][k]))*(sin(phi[j][k][l]))
                   if tau < -1.0:
                      tau = degrees(acos(-1.0))
                   elif tau > 1.0:
                      tau = degrees(acos(1.0))
                   else:
		      tau = degrees(acos(tau))
                   print i,j,k,l,tau

# ************* Center of mass *****************

m = [12.000000000,12.000000000,15.994914619,1.007825032,1.007825032,1.007825032,1.007825032]
i,m_sum,X_cm,Y_cm,Z_cm = 0,0,0,0,0
while i < 7:
      m_sum = m_sum + m[i]
      X_cm =  X_cm + float(x[i])*m[i]
      Y_cm =  Y_cm + float(y[i])*m[i]
      Z_cm =  Z_cm + float(z[i])*m[i]
      i += 1
X_cm = X_cm/m_sum
Y_cm = Y_cm/m_sum
Z_cm = Z_cm/m_sum
outdata.write('\n'+'Molecular center of mass: '+'%.8f'%X_cm+' '*3+'%.8f'%Y_cm+' '*3+'%.8f'%Z_cm)

#************* Moment of inertia ***************

Iner = np.zeros((3,3))
i = 0
while i < 7:
    Iner[0][0] = Iner[0][0] + m[i]*(((float(y[i])-Y_cm)**2)+((float(z[i])-Z_cm)**2))
    Iner[1][1] = Iner[1][1] + m[i]*((float(x[i])-X_cm)**2+(float(z[i])-Z_cm)**2)
    Iner[2][2] = Iner[2][2] + m[i]*((float(x[i])-X_cm)**2+(float(y[i])-Y_cm)**2)
    Iner[0][1] = Iner[1][0] = Iner[0][1] + m[i]*(float(x[i])-X_cm)*(float(y[i])-Y_cm)
    Iner[0][2] = Iner[2][0] = Iner[0][2] + m[i]*(float(x[i])-X_cm)*(float(z[i])-Z_cm)
    Iner[1][2] = Iner[2][1] = Iner[1][2] + m[i]*(float(y[i])-Y_cm)*(float(z[i])-Z_cm)
    i += 1
outdata.write(2*'\n'+'Moment of inertia tensor (amu bohr^2):')
outdata.write('\n'+'%.10f'%Iner[0][0]+3*' '+'%.10f'%Iner[0][1]+3*' '+'%.10f'%Iner[0][2])
outdata.write('\n'+'%.10f'%Iner[0][1]+' '*3+'%.10f'%Iner[1][1]+' '*3+'%.10f'%Iner[1][2])
outdata.write('\n'+'%.10f'%Iner[0][2]+' '*3+'%.10f'%Iner[1][2]+' '*3+'%.10f'%Iner[2][2])

# Diagonalization of matrix

dia_Iner = np.linalg.eigvals(Iner)
I = sorted(dia_Iner)
outdata.write(2*'\n'+'Principal moments of inertia (amu * bohr^2):\n')
outdata.write(9*' '+'%.6f'%I[0]+3*' '+'%.6f'%I[1]+3*' '+'%.6f'%I[2])
outdata.write(2*'\n'+'Principal moments of inertia (amu * AA^2):\n')
conv = 0.529177249 * 0.529177249
I1 = [i*conv for i in I]
outdata.write(9*' '+'%.6f'%I1[0]+3*' '+'%.6f'%I1[1]+3*' '+'%.6f'%I1[2])
outdata.write(2*'\n'+'Principal moments of inertia (g * cm^2):\n')
conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8
I2 = [i*conv for i in I]
outdata.write(9*' '+'%10.6E'%I2[0]+3*' '+'%10.6E'%I2[1]+3*' '+'%10.6E'%I2[2]+2*'\n')

# type of rotor

if I[0]==I[1] and I[1]==I[2]:
   outdata.write('Molecule is a spherical top.\n')
elif I[0]==0 and I[1]==I[2]:
   outdata.write('Molecule is a linear rotor.\n')
elif I[0]==I[1] and I[1]<I[2]:
   outdata.write('Molecule is an oblate symmetric top.\n')
elif I[1]==I[2] and I[1]>I[0]:
   outdata.write('Molecule is a prolate symmetric top.\n')
else:
   outdata.write('Molecule is an asymmetric top.\n')

# *************** Rotational Constants ***************

outdata.write('\nRotational constants (cm-1):\n')
A = ((6.6260755E-34)*(1E7))/(8*((pi)**2)*(2.99792458E10)*(I2[0]))
B = ((6.6260755E-34)*(1E7))/(8*((pi)**2)*(2.99792458E10)*(I2[1]))
C = ((6.6260755E-34)*(1E7))/(8*((pi)**2)*(2.99792458E10)*(I2[2]))
outdata.write(9*' '+'A = '+'%.4f'%A+3*' '+'B = '+'%.4f'%B+3*' '+'C = '+'%.4f'%C+'\n'*2)

outdata.write('\nRotational constants (MHz):\n')
A = A*(2.99792458E10)*(1E-6)
B = B*(2.99792458E10)*(1E-6)
C = C*(2.99792458E10)*(1E-6)
outdata.write(9*' '+'A = '+'%.3f'%A+3*' '+'B = '+'%.3f'%B+3*' '+'C = '+'%.3f'%C)
# ******************************************************************
