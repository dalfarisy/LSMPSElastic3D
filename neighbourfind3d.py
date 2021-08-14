"""
LSMPSNeighbourfind3D v1.0
By Dew
"""

import numpy as np
import array as arr

#%% Domain Divider

def domaindivide():
    print("Dividing Domain...")
    global NcellX,NcellY,NcellZ,DivDomain
    xmin = min(x)-0.00000001
    xmax = max(x)+0.00000001
    ymin = min(y)-0.00000001
    ymax = max(y)+0.00000001
    zmin = min(z)-0.00000001
    zmax = max(z)+0.00000001    
    
    dY = rc
    dX = rc
    dZ = rc
    
    NcellY = int(np.ceil((ymax-ymin)/dY))
    NcellX = int(np.ceil((xmax-xmin)/dX))
    NcellZ = int(np.ceil((zmax-zmin)/dZ))    
    
    DivDomain = [0]*NcellX
    
    for i in range(len(DivDomain)):
        DivDomain[i] = [0]*NcellY
        for j in range(len(DivDomain[i])):
            DivDomain[i][j]=[0]*NcellZ
            for k in range(len(DivDomain[i][j])):
                DivDomain[i][j][k] = [[],[],[]]
    
    # Finding the domain's neighbour
    for i in range(NcellX):
        for j in range(NcellY):
            for k in range(NcellZ):
                if i != 0:
                    DivDomain[i][j][k][0].append([i-1,j,k])
                if j != 0:
                    DivDomain[i][j][k][0].append([i,j-1,k])
                if k != 0:
                    DivDomain[i][j][k][0].append([i,j,k-1])
                if i != 0 and j != 0:
                    DivDomain[i][j][k][0].append([i-1,j-1,k])
                if i != 0 and k != 0:
                    DivDomain[i][j][k][0].append([i-1,j,k-1])
                if j != 0 and k != 0:
                    DivDomain[i][j][k][0].append([i,j-1,k-1])
                if i!=0 and j != 0 and k != 0:
                    DivDomain[i][j][k][0].append([i-1,j-1,k-1])
                if i != NcellX-1:
                    DivDomain[i][j][k][0].append([i+1,j,k])
                if j != NcellY-1:
                    DivDomain[i][j][k][0].append([i,j+1,k])
                if k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i,j,k+1])
                if i != NcellX-1 and j != 0:
                    DivDomain[i][j][k][0].append([i+1,j-1,k])
                if i != 0 and j != NcellY-1:
                    DivDomain[i][j][k][0].append([i-1,j+1,k])
                if i != NcellX-1 and k != 0:
                    DivDomain[i][j][k][0].append([i+1,j,k-1])
                if i != 0 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i-1,j,k+1])
                if j != NcellY-1 and k != 0:
                    DivDomain[i][j][k][0].append([i,j+1,k-1])
                if j != 0 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i,j-1,k+1])
                if i != NcellX-1 and j != NcellY-1:
                    DivDomain[i][j][k][0].append([i+1,j+1,k])
                if i != NcellX-1 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i+1,j,k+1])
                if j != NcellY-1 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i,j+1,k+1])
                if i != NcellX-1 and j != NcellY-1 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i+1,j+1,k+1])
                if i != 0 and j != NcellY-1 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i-1,j+1,k+1])  
                if i != NcellX-1 and j != 0 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i+1,j-1,k+1])
                if i != NcellX-1 and j != NcellY-1 and k != 0:
                    DivDomain[i][j][k][0].append([i+1,j+1,k-1])
                if i != 0 and j != 0 and k != NcellZ-1:
                    DivDomain[i][j][k][0].append([i-1,j-1,k+1]) 
                if i != 0 and j != NcellY-1 and k != 0:
                    DivDomain[i][j][k][0].append([i-1,j+1,k-1])
                if i != NcellX-1 and j != 0 and k != 0:
                    DivDomain[i][j][k][0].append([i+1,j-1,k-1])   
    
    # Finding the minimum and maximum coordinates for each domain
    for i in range(NcellX):
        for j in range(NcellY):
            for k in range(NcellZ):
                xstart = xmin + i*dX
                xend = xmin + (i+1)*dX
                ystart = ymin + j*dY
                yend = ymin + (j+1)*dY
                zstart = zmin + k*dZ
                zend = zmin + (k+1)*dZ
                DivDomain[i][j][k][1] = [[xstart,xend],[ystart,yend],[zstart,zend]]
            
    # Filling each domain with particles
    for k in range(len(x)):
        for i in range(NcellX):
            for j in range(NcellY):
                for l in range(NcellZ):
                    if x[k] > DivDomain[i][j][l][1][0][0] and x[k] < DivDomain[i][j][l][1][0][1] and y[k] > DivDomain[i][j][l][1][1][0] and y[k] < DivDomain[i][j][l][1][1][1] and z[k] > DivDomain[i][j][l][1][2][0] and z[k] < DivDomain[i][j][l][1][2][1]:
                        DivDomain[i][j][l][2].append(k)
                    
#%% Sorting Algorithm
                    
def sorting(sorter,sortdata):
    c = [x for _,x in sorted(zip(sorter,sortdata))]
    return c

#%% New Neighbour Find Algorithm With Domain Division

def neighbourfindalgorithm():
    print("Finding Neighbour...")
    for i in range(len(x)):
        neighbour.append([])
    for i in range(NcellX):
        for j in range(NcellY):
            totaldomain = []
            for l in range(len(DivDomain[i][j][2])):
                totaldomain.append(DivDomain[i][j][2][l])
            for k in range(len(DivDomain[i][j][0])):
                m = DivDomain[i][j][0][k][0]
                n = DivDomain[i][j][0][k][1]
                for l in range(len(DivDomain[m][n][2])):
                    totaldomain.append(DivDomain[m][n][2][l])
            for k in range(len(DivDomain[i][j][2])):
                for l in range(len(totaldomain)):
                    if DivDomain[i][j][2][k] != totaldomain[l]:
                        dx = x[DivDomain[i][j][2][k]] - x[totaldomain[l]]
                        dy = y[DivDomain[i][j][2][k]] - y[totaldomain[l]]
                        if np.sqrt(dx**(2) + dy**(2)) < rc:
                            neighbour[DivDomain[i][j][2][k]].append(totaldomain[l])

#%% New Neighbour Find Algorithm With Domain Division and maximum number of neighbour limit

def neighbourfindalgorithmmaxlimit(targetneighbour):
    print("Finding Neighbour...")
    for i in range(len(x)):
        neighbour.append([])
    for i in range(NcellX):
        for j in range(NcellY):
            for h in range(NcellZ):
                totaldomain = []
                for l in range(len(DivDomain[i][j][h][2])):
                    totaldomain.append(DivDomain[i][j][h][2][l])
                for k in range(len(DivDomain[i][j][h][0])):
                    m = DivDomain[i][j][h][0][k][0]
                    n = DivDomain[i][j][h][0][k][1]
                    o = DivDomain[i][j][h][0][k][2]
                    for l in range(len(DivDomain[m][n][o][2])):
                        totaldomain.append(DivDomain[m][n][o][2][l])
                for k in range(len(DivDomain[i][j][h][2])):
                    temp = []
                    temp1 = []
                    sortedtemp = []
                    for l in range(len(totaldomain)):
                        if DivDomain[i][j][h][2][k] != totaldomain[l]:
                            dx = x[DivDomain[i][j][h][2][k]] - x[totaldomain[l]]
                            dy = y[DivDomain[i][j][h][2][k]] - y[totaldomain[l]]
                            dz = z[DivDomain[i][j][h][2][k]] - z[totaldomain[l]]
                            radius = np.sqrt(dx**(2) + dy**(2) + dz**(2))
                            if radius < rc:
                                temp.append(totaldomain[l])
                                temp1.append(radius)
                            
                    sortedtemp = sorting(temp1,temp)
                    for p in range(targetneighbour):
                        neighbour[DivDomain[i][j][h][2][k]].append(sortedtemp[p])
                    
        
                            
                    
#%% Main Function
                            
def neighbourfind(x1,y1,z1,rc1):
    global x,y,z,neighbour,rc
    neighbour = []
    x = x1
    y = y1
    z = z1
    rc = rc1
    domaindivide()
    neighbourfindalgorithm()
    neighbournum = []
    for i in range(len(neighbour)):
        neighbournum.append(len(neighbour[i]))
    print(f"Maximum Neighbour: {max(neighbournum)}")
    print(f"Minimum Neighbour: {min(neighbournum)}")
    return neighbour

def neighbourfindlimited(x1,y1,z1,rc1,targetnumber):
    global x,y,z,neighbour,rc
    neighbour = []
    x = x1
    y = y1
    z = z1
    rc = rc1
    domaindivide()
    neighbourfindalgorithmmaxlimit(int(targetnumber))
    neighbournum = []
    for i in range(len(neighbour)):
        neighbournum.append(len(neighbour[i]))
    print(f"Maximum Neighbour: {max(neighbournum)}")
    print(f"Minimum Neighbour: {min(neighbournum)}")
    return neighbour
