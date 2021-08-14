##LSMPS General v1.1 by Devvskiii


import numpy as np
import math
import datetime
import array as arr

#%% Standard general LSMPS function
def LSMPS(x,y,z,w,neighbour,Rs):
    
    print("LSMPS Operator Calculating ...")
    
    start = datetime.datetime.now()
    
    result = []
    
    Hrs = np.zeros((10,10))
    # P = [1,x,y,z,xy,xz,yz,x^2,y^2,z^2]
    m = 0; n = 0; o = 0
    Hrs[0][0] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 0; o = 0
    Hrs[1][1] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 1; o = 0
    Hrs[2][2] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 0; o = 1
    Hrs[3][3] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 1; o = 0
    Hrs[4][4] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 0; o = 1
    Hrs[5][5] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 1; o = 1
    Hrs[6][6] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 2; n = 0; o = 0
    Hrs[7][7] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 2; o = 0
    Hrs[8][8] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 0; o = 2
    Hrs[9][9] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    
    for i in range(len(x)):
        neighbourlen = len(neighbour[i])
        
        Rdata = arr.array('d',[])
        for j in range(neighbourlen):
            Rdata.append(np.sqrt((x[i]-x[neighbour[i][j]])**(2) + (y[i]-y[neighbour[i][j]])**(2) + (z[i]-z[neighbour[i][j]])**(2)))
        
        Rmax = max(Rdata)
        
        weight = np.zeros(neighbourlen)
        for j in range(neighbourlen):
            xn = x[neighbour[i][j]]
            yn = y[neighbour[i][j]]
            zn = z[neighbour[i][j]]
            
            weight[j] = weighting(x[i],xn,y[i],yn,z[i],zn,Rmax)
        
        M = np.zeros((10,10))
        P = np.zeros((10,1))
        b = np.zeros((10,1))
        
        for j in range(neighbourlen):
            # P = [1,x,y,z,xy,xz,yz,x^2,y^2,z^2]
            
            xn = x[neighbour[i][j]]
            yn = y[neighbour[i][j]]
            zn = z[neighbour[i][j]]
            
            dx1 = xn - x[i]
            dy1 = yn - y[i]
            dz1 = zn - z[i]
            
            dx2 = dx1/Rs
            dy2 = dy1/Rs
            dz2 = dz1/Rs
            
            P[0][0] = 1
            P[1][0] = dx2
            P[2][0] = dy2
            P[3][0] = dz2
            P[4][0] = dx2*dy2
            P[5][0] = dx2*dz2
            P[6][0] = dy2*dz2
            P[7][0] = dx2**(2)
            P[8][0] = dy2**(2)
            P[9][0] = dz2**(2)
            
            M = M + weight[j]*np.matmul(P,P.T)
            
            b = b + weight[j]*P*(w[neighbour[i][j]]-w[i])
            
        
        Minv = np.linalg.inv(M) 
        
        tempresult = np.matmul(Hrs,np.matmul(Minv,b))
        
        result.append(tempresult)
        
    end = datetime.datetime.now()
    
    elapsedtime = end - start
    
    print(f"LSMPS Operator Time: {int(elapsedtime.total_seconds()*1000)} ms")
        
    return result

#%% Weighting Function
def weighting(x1,x2,y1,y2,z1,z2,Re):
    R = np.sqrt((x2-x1)**(2) + (y2-y1)**(2) + (z2-z1)**(2))
    if R < Re:
        #weight = Re/R - 1
        weight = (1 - R/Re)**(2)
    elif R >= Re:
        weight = 0
    return weight

#%% Cubic Spline Kernel
def cubespline(x1,x2,y1,y2,z1,z2,Re):
    Rnorm = np.sqrt((x2-x1)**(2) + (y2-y1)**(2) + (z2-z1)**(2))/Re
    if Rnorm > 1:
        weight = 0
    elif Rnorm <= 1/2:
        weight = 2/3 - 4*Rnorm**(2) + 4*Rnorm**(3)
    else:
        weight = 4/3 - 4*Rnorm + 4*Rnorm**(2) - 4/3*Rnorm**(3)
    return weight

#%% LSMPS Implicit Form
def LSMPSconst(x,y,z,neighbour,Rs):
    
    print("LSMPS Constants Calculation ...")
    
    start = datetime.datetime.now()
    
    result = []
    
    for i in range(len(x)):
        neighbourlen = len(neighbour[i])
        
        Rdata = arr.array('d',[])
        for j in range(neighbourlen):
            Rdata.append(np.sqrt((x[i]-x[neighbour[i][j]])**(2) + (y[i]-y[neighbour[i][j]])**(2) + (z[i]-z[neighbour[i][j]])**(2)))
        
        Rmax = max(Rdata)
        
        weight = np.zeros(neighbourlen)
        
        for j in range(neighbourlen):
            xn = x[neighbour[i][j]]
            yn = y[neighbour[i][j]]
            zn = z[neighbour[i][j]]
            
            weight[j] = weighting(x[i],xn,y[i],yn,z[i],zn,Rmax)
        
        M = np.zeros((10,10))
        P = np.zeros((10,1))
        B = np.zeros((10,neighbourlen))
        
        B1 = np.zeros((10,neighbourlen))
            
        const = np.zeros((10,neighbourlen))
        
        for j in range(neighbourlen):
            # P = [1,x,y,xy,x^2,y^2,x^2y,xy^2,x^2y^2]
            
            xn = x[neighbour[i][j]]
            yn = y[neighbour[i][j]]
            zn = z[neighbour[i][j]]
            
            dx1 = xn - x[i]
            dy1 = yn - y[i]
            dz1 = zn - z[i]
            
            dx2 = dx1/Rs
            dy2 = dy1/Rs
            dz2 = dz1/Rs
            
            P[0][0] = 1
            P[1][0] = dx2
            P[2][0] = dy2
            P[3][0] = dz2
            P[4][0] = dx2*dy2
            P[5][0] = dx2*dz2
            P[6][0] = dy2*dz2
            P[7][0] = dx2**(2)
            P[8][0] = dy2**(2)
            P[9][0] = dz2**(2)
            
            M = M + weight[j]*np.matmul(P,P.T)
            
        Minv = np.linalg.inv(M)    
                    
        for j in range(neighbourlen):
            
            xn = x[neighbour[i][j]]
            yn = y[neighbour[i][j]]
            zn = z[neighbour[i][j]]
            
            dx1 = xn - x[i]
            dy1 = yn - y[i]
            dz1 = zn - z[i]
            
            dx2 = dx1/Rs
            dy2 = dy1/Rs
            dz2 = dz1/Rs
            
            B[0][j] = 1*weight[j]
            B[1][j] = dx2*weight[j]
            B[2][j] = dy2*weight[j]
            B[3][j] = dz2*weight[j]
            B[4][j] = dx2*dy2*weight[j]
            B[5][j] = dx2*dz2*weight[j]
            B[6][j] = dy2*dz2*weight[j]
            B[7][j] = dx2**(2)*weight[j]
            B[8][j] = dy2**(2)*weight[j]
            B[9][j] = dz2**(2)*weight[j]
            
            for k in range(10):
                for l in range(10):
                    B1[k][j] = B1[k][j] + B[l][j]*Minv[k][l]
                                        
        result.append(B1)
        
    end = datetime.datetime.now()
    
    elapsedtime = end - start
    
    print(f"LSMPS Operator Time: {int(elapsedtime.total_seconds()*1000)} ms")
        
    return result

def LSMPSHrs(Rs):
    
    print("LSMPS Hrs Calculation ...")
    
    start = datetime.datetime.now()
    
    Hrs = np.zeros(10)
    # P = [1,x,y,z,xy,xz,yz,x^2,y^2,z^2]
    m = 0; n = 0; o = 0
    Hrs[0] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 0; o = 0
    Hrs[1] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 1; o = 0
    Hrs[2] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 0; o = 1
    Hrs[3] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 1; o = 0
    Hrs[4] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 1; n = 0; o = 1
    Hrs[5] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 1; o = 1
    Hrs[6] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 2; n = 0; o = 0
    Hrs[7] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 2; o = 0
    Hrs[8] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    m = 0; n = 0; o = 2
    Hrs[9] = Rs**(-(m+n+o)) *(math.factorial(m)*math.factorial(n)*math.factorial(o))
    
    end = datetime.datetime.now()
    
    elapsedtime = end - start
    
    print(f"LSMPS Hrs Time: {int(elapsedtime.total_seconds()*1000)} ms")
        
    return Hrs