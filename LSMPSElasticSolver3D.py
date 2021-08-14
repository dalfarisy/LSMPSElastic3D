#LSMPSElasticSolver3D v.1.1
#Created by Devvskiii

import numpy as np
import array as arr
import scipy as sc
import scipy.sparse.linalg as SClinalg
import math
from neighbourfind3d import neighbourfind
from neighbourfind3d import neighbourfindlimited
from LSMPSgeneral import LSMPS
from LSMPSgeneral import LSMPSconst
from LSMPSgeneral import LSMPSHrs
import datetime


start = datetime.datetime.now()

#%% Opening the geometry data file

file = open("geom.txt","r")
readdata = file.readlines()
file.close()
E = float(readdata[0])
v = float(readdata[1])
h = float(readdata[2])

x = []
y = []
z = []
nx = []
ny = []
nz = []
dispX = []
dispY = []
dispZ = []
forceX = []
forceY = []
forceZ = []
bx = []
by = []
bz = []

for i in range(3,len(readdata)):
    temp = readdata[i].split("\t")
    x.append(float(temp[0]))
    y.append(float(temp[1]))
    z.append(float(temp[2]))
    if temp[3] != 'nan':
        nx.append(float(temp[3]))
    elif temp[3] == 'nan':
        nx.append('nan')
    if temp[4] != 'nan':
        ny.append(float(temp[4]))
    elif temp[4] == 'nan':
        ny.append('nan')
    if temp[5] != 'nan':
        nz.append(float(temp[5]))
    elif temp[5] == 'nan':
        nz.append('nan')    
    if temp[6] != 'nan':
        dispX.append(float(temp[6]))
    elif temp[6] == 'nan':
        dispX.append('nan')
    if temp[7] != 'nan':
        dispY.append(float(temp[7]))
    elif temp[7] == 'nan':
        dispY.append('nan')
    if temp[8] != 'nan':
        dispZ.append(float(temp[8]))
    elif temp[8] == 'nan':
        dispZ.append('nan')    
    if temp[9] != 'nan':
        forceX.append(float(temp[9]))
    elif temp[9] == 'nan':
        forceX.append('nan')
    if temp[10] != 'nan':
        forceY.append(float(temp[10]))
    elif temp[10] == 'nan':
        forceY.append('nan')
    if temp[11] != 'nan':
        forceZ.append(float(temp[11]))
    elif temp[11] == 'nan':
        forceZ.append('nan')
    if temp[12] != 'nan':
        bx.append(float(temp[12]))
    elif temp[12] == 'nan':
        bx.append('nan')
    if temp[13] != 'nan':
        by.append(float(temp[13]))
    elif temp[13] == 'nan':
        by.append('nan')
    if temp[14] != 'nan' + '\n':
        bz.append(float(temp[14]))
    elif temp[14] == 'nan' + '\n':
        bz.append('nan')

maxnumx = (max(x) - min(x))/h
maxnumy = (max(y) - min(y))/h
maxnumz = (max(z) - min(z))/h
maxnum = max([maxnumx,maxnumy,maxnumz])

print(f"Number of Particles: {len(x)}")
    
#%% Solver Control

scattersize = 250**(2)/(maxnum)**(2)
# The sizing of each node in the scatter plot

neighboursearchmethod = 2
# 1 for constant cutoff radius
# 2 for constant number of neighbours

scaledcutoffradius = 6.5001
# The cutoff radius for neighbour search

targetneighbournumber = 50
# Number of neighbour (for neighboursearchmethod = 2 only)

boundaryplotmethod = 1
# 1 for plotting all the nodes
# 2 for plotting without dirichlet(displacement) boundary nodes
# 3 for plotting without both dirichlet and neumann boundary nodes

dispscale = 200
# The displacement scale for plotting

datatype = 'd'
# 'f' for float
# 'd' for double

convertunits = 2
# 1 to keep current units
# 2 to convert unit from mm to m

writeresult = 2
# 1 to not write the result file
# 2 to write the result file

#%% Material data calculation

lamda = v*E/((1+v)*(1-2*v))
Mu = E/(2*(1+v))

#%% Boundary Search

neighbourstarttime = datetime.datetime.now()
if neighboursearchmethod == 1:
    neighbour = neighbourfind(x,y,z,scaledcutoffradius*h)
if neighboursearchmethod == 2:
    neighbour = neighbourfindlimited(x,y,z,scaledcutoffradius*h,targetneighbournumber)
else:
    print("Neighbour search method selection INVALID")

for i in range(len(x)):
    neighbour[i].append(i)
  
neighbourendtime = datetime.datetime.now()

#%% Solving for LSMPS coefficients

calcstarttime = datetime.datetime.now()
e = h
# P = [1,x,y,z,xy,xz,yz,x^2,y^2,z^2]
LSMPSconstant = LSMPSconst(x,y,z,neighbour,e)
LSMPSHRS = LSMPSHrs(e)
 
#%% Rearranging The Coefficients in Matrix Form

print('Rearranging the Matrix ...')

b = np.zeros(len(x)*3)
row = arr.array('I',[])
col = arr.array('I',[])
data = arr.array(datatype,[])

for i in range(len(x)):
    
    # Dirichlet Boundary
    if dispX[i] != 'nan':
        
        row.append(3*i)
        col.append(3*i)
        data.append(1.0)
        
        b[3*i] = dispX[i]
        
    if dispY[i] != 'nan':
        
        row.append(3*i+1)
        col.append(3*i+1)
        data.append(1.0)
        
        b[3*i+1] = dispY[i]
        
    if dispZ[i] != 'nan':
        
        row.append(3*i+2)
        col.append(3*i+2)
        data.append(1.0)
        
        b[3*i+2] = dispZ[i]
        
    # Neumann Boundary
    if forceX[i] != 'nan':
        
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        
        for j in range(len(neighbour[i])):
            
            temp1 = temp1 + (lamda + 2*Mu)*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
            temp2 = temp2 + lamda*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1]
            
            temp3 = temp3 + lamda*nx[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nz[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1]
            
            row.append(3*i)
            col.append(3*neighbour[i][j])
            data.append((lamda + 2*Mu)*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3])
        
            row.append(3*i)
            col.append(3*neighbour[i][j]+1)
            data.append(lamda*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1])
            
            row.append(3*i)
            col.append(3*neighbour[i][j]+2)
            data.append(lamda*nx[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nz[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1])
            
            b[3*i] = forceX[i]
            
            if j == len(neighbour[i])-1:
                
                row.append(3*i)
                col.append(3*i)
                data.append(temp1)
                
                row.append(3*i)
                col.append(3*i+1)
                data.append(temp2)
                
                row.append(3*i)
                col.append(3*i+2)
                data.append(temp3)
            
    if forceY[i] != 'nan':
        
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        
        for j in range(len(neighbour[i])):
            
            temp1 = temp1 + lamda*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2]
            
            temp2 = temp2 + (lamda + 2*Mu)*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
            temp3 = temp3 + lamda*ny[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nz[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2]
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j])
            data.append(lamda*ny[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2])
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j]+1)
            data.append((lamda + 2*Mu)*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3])
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j]+2)
            data.append(lamda*ny[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nz[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2])
            
            b[3*i+1] = forceY[i]
            
            if j == len(neighbour[i])-1:
                
                row.append(3*i+1)
                col.append(3*i)
                data.append(temp1)
                
                row.append(3*i+1)
                col.append(3*i+1)
                data.append(temp2)
                
                row.append(3*i+1)
                col.append(3*i+2)
                data.append(temp3)
                
    if forceZ[i] != 'nan':
        
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        
        for j in range(len(neighbour[i])):
            
            temp1 = temp1 + lamda*nz[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
            temp2 = temp2 + lamda*nz[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3]
            
            temp3 = temp3 + (lamda + 2*Mu)*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2]
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j])
            data.append(lamda*nz[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*nx[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3])
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j]+1)
            data.append(lamda*nz[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2] + Mu*ny[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3])
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j]+2)
            data.append((lamda + 2*Mu)*nz[i]*LSMPSconstant[i][3][j]*LSMPSHRS[3] + Mu*nx[i]*LSMPSconstant[i][1][j]*LSMPSHRS[1] + Mu*ny[i]*LSMPSconstant[i][2][j]*LSMPSHRS[2])
            
            b[3*i+2] = forceZ[i]
            
            if j == len(neighbour[i])-1:
                
                row.append(3*i+2)
                col.append(3*i)
                data.append(temp1)
                
                row.append(3*i+2)
                col.append(3*i+1)
                data.append(temp2)
                
                row.append(3*i+2)
                col.append(3*i+2)
                data.append(temp3)
            
    if forceX[i] == 'nan' and forceY[i] == 'nan' and forceZ[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan' and dispZ[i] == 'nan':
        
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        temp4 = 0.0
        temp5 = 0.0
        temp6 = 0.0
        temp7 = 0.0
        temp8 = 0.0
        temp9 = 0.0
        
        for j in range(len(neighbour[i])):
            
            
            # X direction
            temp1 = temp1 + (lamda + 2*Mu)*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][8][j]*LSMPSHRS[8] + Mu*LSMPSconstant[i][9][j]*LSMPSHRS[9]
            
            temp2 = temp2 + (lamda + Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4]
            
            temp3 = temp3 + (lamda + Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5]

            row.append(3*i)
            col.append(3*neighbour[i][j])
            data.append((lamda + 2*Mu)*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][8][j]*LSMPSHRS[8] + Mu*LSMPSconstant[i][9][j]*LSMPSHRS[9])
            
            row.append(3*i)
            col.append(3*neighbour[i][j]+1)
            data.append((lamda + Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4])
            
            row.append(3*i)
            col.append(3*neighbour[i][j]+2)
            data.append((lamda + Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5])
            
            b[3*i] = bx[i]
            
            
            # Y direction
            temp4 = temp4 + (lamda + Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4]
            
            temp5 = temp5 + (lamda + 2*Mu)*LSMPSconstant[i][8][j]*LSMPSHRS[8] + Mu*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][9][j]*LSMPSHRS[9]
            
            temp6 = temp6 + (lamda + Mu)*LSMPSconstant[i][6][j]*LSMPSHRS[6]
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j])
            data.append((lamda + Mu)*LSMPSconstant[i][4][j]*LSMPSHRS[4])
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j]+1)
            data.append((lamda + 2*Mu)*LSMPSconstant[i][8][j]*LSMPSHRS[8] + Mu*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][9][j]*LSMPSHRS[9])
            
            row.append(3*i+1)
            col.append(3*neighbour[i][j]+2)
            data.append((lamda + Mu)*LSMPSconstant[i][6][j]*LSMPSHRS[6])
            
            b[3*i+1] = by[i]
            
            
            # Z direction
            temp7 = temp7 + (lamda + Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5]
            
            temp8 = temp8 + (lamda + Mu)*LSMPSconstant[i][6][j]*LSMPSHRS[6]
            
            temp9 = temp9 + (lamda + 2*Mu)*LSMPSconstant[i][9][j]*LSMPSHRS[9] + Mu*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][8][j]*LSMPSHRS[8]
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j])
            data.append((lamda + Mu)*LSMPSconstant[i][5][j]*LSMPSHRS[5])
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j]+1)
            data.append((lamda + Mu)*LSMPSconstant[i][6][j]*LSMPSHRS[6])
            
            row.append(3*i+2)
            col.append(3*neighbour[i][j]+2)
            data.append((lamda + 2*Mu)*LSMPSconstant[i][9][j]*LSMPSHRS[9] + Mu*LSMPSconstant[i][7][j]*LSMPSHRS[7] + Mu*LSMPSconstant[i][8][j]*LSMPSHRS[8])
            
            b[3*i+2] = bz[i]
            
            
            if j == len(neighbour[i])-1:
                
                row.append(3*i)
                col.append(3*i)
                data.append(temp1)
                
                row.append(3*i)
                col.append(3*i+1)
                data.append(temp2)
                
                row.append(3*i)
                col.append(3*i+2)
                data.append(temp3)
                
                row.append(3*i+1)
                col.append(3*i)
                data.append(temp4)
                
                row.append(3*i+1)
                col.append(3*i+1)
                data.append(temp5)
                
                row.append(3*i+1)
                col.append(3*i+2)
                data.append(temp6)
                
                row.append(3*i+2)
                col.append(3*i)
                data.append(temp7)
                
                row.append(3*i+2)
                col.append(3*i+1)
                data.append(temp8)
                
                row.append(3*i+2)
                col.append(3*i+2)
                data.append(temp9)
            
#%% Calculating Displacement Vector and Stresses
print('Matrix Rearrangement Complete!')
A = sc.sparse.csr_matrix((data, (row, col)), shape=(len(x)*3, len(x)*3))
print('Solving...')
startsolve = datetime.datetime.now()
U = SClinalg.spsolve(A,b,use_umfpack=True)
endsolve = datetime.datetime.now()

print('Solving Complete!')
Ux = arr.array(datatype,[])
Uy = arr.array(datatype,[])
Uz = arr.array(datatype,[])
for i in range(len(U)):
    if i%3 == 0:
        Ux.append(U[i])
    if i%3 == 1:
        Uy.append(U[i])
    if i%3 == 2:
        Uz.append(U[i])
    
dUx = np.zeros((len(x),10))
#dUxdY = np.zeros(len(x))
dUy = np.zeros((len(x),10))
#dUydY = np.zeros(len(x))
dUz = np.zeros((len(x),10))

epsxx = arr.array(datatype,[])
epsyy = arr.array(datatype,[])
epszz = arr.array(datatype,[])
epsxy = arr.array(datatype,[])
epsxz = arr.array(datatype,[])
epsyz = arr.array(datatype,[])
sigmax = arr.array(datatype,[])
sigmay = arr.array(datatype,[])
sigmaz = arr.array(datatype,[])
sigmaxy = arr.array(datatype,[])
sigmaxz = arr.array(datatype,[])
sigmayz = arr.array(datatype,[])
vonmises = arr.array(datatype,[])

dUx = LSMPS(x,y,z,Ux,neighbour,e)
#dUxdY = LSMPS(x,y,Ux,neighbour,e)
dUy = LSMPS(x,y,z,Uy,neighbour,e)
#dUydY = LSMPS(x,y,Uy,neighbour,e)
dUz = LSMPS(x,y,z,Uz,neighbour,e)

# =============================================================================
# print(f"Max dUxdX:{max(dUxdX)}")
# print(f"Min dUxdX:{min(dUxdX)}")
# print(f"Max dUxdY:{max(dUxdY)}")
# print(f"Min dUxdY:{min(dUxdY)}")
# print(f"Max dUydX:{max(dUydX)}")
# print(f"Min dUydX:{min(dUydX)}")
# print(f"Max dUydY:{max(dUydY)}")
# print(f"Min dUydY:{min(dUydY)}")
# =============================================================================

j = 0

if boundaryplotmethod == 1:
    for i in range(len(x)):
        epsxx.append(dUx[i][1])
        epsyy.append(dUy[i][2])
        epszz.append(dUz[i][3])
        epsxy.append(0.5*(dUx[i][2] + dUy[i][1]))
        epsxz.append(0.5*(dUx[i][3] + dUz[i][1]))
        epsyz.append(0.5*(dUy[i][3] + dUz[i][2]))
        
        sigmax.append(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsxx[j]))
        sigmay.append(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epsyy[j]))
        sigmaz.append(lamda*(epsxx[j] + epsyy[j] + epszz[j]) + 2*Mu*(epszz[j]))
        sigmaxy.append(2*Mu*epsxy[j])
        sigmaxz.append(2*Mu*epsxz[j])
        sigmayz.append(2*Mu*epsyz[j])
        vonmises.append(np.sqrt(sigmax[j]**(2) + sigmay[j]**(2) + sigmaz[j]**(2) - sigmax[j]*sigmay[j] - sigmax[j]*sigmaz[j] - sigmay[j]*sigmaz[j] + 3*(sigmaxy[j]**(2) + sigmaxz[j]**(2) + sigmayz[j]**(2))))
        j = j + 1  

else:
    print("Boundary plot method selection INVALID")
    
#%% Writing the result file

if writeresult == 2:
    
    file = open("result.txt","w")
    file.write(f"{h}")
    file.write("\n")
    j = 0
    
    if boundaryplotmethod == 1:
        for i in range(len(x)):
            file.write(f"{x[i]}" + "\t" + f"{y[i]}"+ "\t" + f"{z[i]}" + "\t" + f"{vonmises[j]}" + "\t" + f"{sigmax[j]}" + "\t" + f"{sigmay[j]}" + "\t" + f"{sigmaz[j]}" +"\t" + f"{sigmaxy[j]}" +"\t" + f"{sigmaxz[j]}" +"\t" + f"{sigmayz[j]}" + "\t" + f"{Ux[i]}" + "\t" + f"{Uy[i]}" + "\t" + f"{Uz[i]}")
            file.write("\n")
            j = j + 1
                
    
    file.close()
    
#%% Calculating the final position
    
Xnew = arr.array(datatype,[])
Ynew = arr.array(datatype,[])
Znew = arr.array(datatype,[])
        
if boundaryplotmethod == 1:
    for i in range(len(x)):
        Xnew.append(x[i] + Ux[i]*dispscale)
        Ynew.append(y[i] + Uy[i]*dispscale)
        Znew.append(z[i] + Uz[i]*dispscale)
elif boundaryplotmethod == 2:
    for i in range(len(x)):
        if dispX[i] == 'nan' and dispY[i] == 'nan':
            Xnew.append(x[i] + Ux[i]*dispscale)
            Ynew.append(y[i] + Uy[i]*dispscale)
            Znew.append(z[i] + Uz[i]*dispscale)
elif boundaryplotmethod == 3:
    for i in range(len(x)):
        if nx[i] == 'nan' and ny[i] == 'nan' and dispX[i] == 'nan' and dispY[i] == 'nan':
            Xnew.append(x[i] + Ux[i]*dispscale)
            Ynew.append(y[i] + Uy[i]*dispscale)
            Znew.append(z[i] + Uz[i]*dispscale)

#%% Converting to SI units (m and Pa)
            
if convertunits == 2:
    for i in range(len(Xnew)):
        Xnew[i] = Xnew[i]/1000
        Ynew[i] = Ynew[i]/1000
        Znew[i] = Ynew[i]/1000
        sigmax[i] = sigmax[i]*1000
        sigmay[i] = sigmay[i]*1000
        sigmaz[i] = sigmaz[i]*1000
        sigmaxy[i] = sigmaxy[i]*1000
        sigmaxz[i] = sigmaxz[i]*1000
        sigmayz[i] = sigmayz[i]*1000
        vonmises[i] = vonmises[i]*1000
    for i in range(len(x)):
        x[i] = x[i]/1000
        y[i] = y[i]/1000
        z[i] = z[i]/1000
        Ux[i] = Ux[i]/1000
        Uy[i] = Uy[i]/1000
        Uz[i] = Uz[i]/1000

calcendtime = datetime.datetime.now()

#%% Plotting The Contour

neighbourtime = neighbourendtime - neighbourstarttime
calculationtime = calcendtime - calcstarttime
solvetime = endsolve - startsolve
totaltime = calcendtime - start
print(f"Neighbour Search Elapsed Time: {int(neighbourtime.total_seconds()*1000)} ms")
print(f"Calculation Elapsed Time: {int(calculationtime.total_seconds()*1000)} ms")
print(f"Solving Elapsed Time: {int(solvetime.total_seconds()*1000)} ms")
print(f"Total Elapsed Time: {int(totaltime.total_seconds()*1000)} ms")
print(f"Max Positive X displacement: {max(Ux)}")
print(f"Max Negative X displacement: {min(Ux)}")
print(f"Max Positive Y displacement: {max(Uy)}")
print(f"Max Negative Y displacement: {min(Uy)}")
print(f"Max Positive Z displacement: {max(Uz)}")
print(f"Max Negative Z displacement: {min(Uz)}")






