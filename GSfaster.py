#!/usr/bin/env python


#   ==========================================================================================================
#   Description:  A program to calculate the Grad-Shafranov equation to find equilibria models
#                 Faster using numpy arrays                  
#
#   Version: 4.0
#   Date: June 2020
#   

#   Equation: GS A = -r**2*sin(theta)*(electron denisty) + Beta*d(Beta)/dA
#                B = strength * (A - A(R,theta))**(power)
#               
#   The code requires you to put an input strength and power
#          
#   Author: Ankan Sur
#   Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
#   ==========================================================================================================


import numpy as np
from math import *
import math
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import sys, os
from optparse import OptionParser
from scipy import special
from matplotlib.ticker import AutoMinorLocator
import timeit

#define grid
u = np.linspace(-1,1,101)
r = np.linspace(0.0,2.0,101)

#find stellar radius
rid = np.where(r>1)[0][0]-1
uid = np.where(u>0)[0][0]-1
        
Nr = len(r)
Nu= len(u)

dr = r[2]-r[1]
du = u[2]-u[1]

#input parameters
strength = 20
power = 1.1

tNum = 10000000

def init(A):
    for i in range(Nr):
        for j in range(Nu):
            A[i,j] = (1-u[j]**2)
    return A


def source(S):
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if r[i]<1.0:
                S[i,j] = (r[i])**2.0*(1-u[j]**2)*(1-r[i]**2)*1.574
            else:
                S[i,j] = 0.0
    return S
    
def source2(Ac,rid,uid):
    SS = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if Ac[i,j]>Ac[rid,uid]:
                SS[i,j] = strength**2*power*(Ac[i,j]-Ac[rid,uid])**(2.0*power-1.0)  
            else:
                SS[i,j] = 0.0    
    return SS
    
    
def derivS(AA,rid,uid):
    dS = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if AA[i,j]>AA[rid,uid]:
                dS[i,j] = strength**2*power*(2*power-1)*(AA[i,j]-AA[rid,uid])**(2.0*power-2.0)  
                #dS[i,j] = strength**2*power*(2*power-1)
            else:
                dS[i,j] = 0.0    
    return dS
    
    
def den(d):
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            d[i,j] = 2.0/dr**2 + 2*(1-u[j]**2)/r[i]**2/du**2
    return d
    
    
def num(nn):
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            nn[i,j] = (1-u[j]**2)/r[i]**2/du**2
    return nn
             
    
def solve_boundary(A):
    for i in range(Nr):
        A[i,0] = 0
        A[i,Nu-1] = 0
    
    for j in range(Nu):
        A[0,j] = 0.0
        A[Nr-1,j] = 0.0#(1-u[j]**2)/r[Nr-1]/2.0#(5*r[0]**3/r[Nr-1] - 3*r[0]**5/r[Nr-1] + 3*r[Nr-1]**4 - 5*r[Nr-1]**2)*(1-u[j]**2)/30.0
           
    return A
    
       
def new_den(d, dS):
    dnew = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if dS[i,j]<0:
                dnew[i,j] = d[i,j]-dS[i,j]
            else:
                dnew[i,j] = d[i,j]  
                
    return dnew 
    

       
def new_source(Ac, S2, dS):
    Snew = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if dS[i,j]>0:
                Snew[i,j] = S2[i,j]
            else:
                Snew[i,j] = S2[i,j]-dS[i,j]*Ac[i,j]  
                
    return Snew 
                     
       
A = np.zeros((Nr,Nu))
A = init(A)
A = solve_boundary(A)

S1 = np.zeros((Nr,Nu))
S1 = source(S1)

d = np.zeros((Nr,Nu))
d = den(d)


nn = np.zeros((Nr,Nu))
nn = num(nn)

Aaux = np.zeros((Nr,Nu))
def L2_error(p, pn):
    return np.sqrt(np.sum((p - pn)**2)/np.sum(pn**2))

start = timeit.default_timer()
w = 1

for n in range(tNum):

    Ac = A.copy()
           
    S2 = source2(Ac,rid,uid)
    dS = derivS(Ac,rid,uid)
    dnew = new_den(d, dS)
    Snew = new_source(Ac, S2, dS)
        
    A[1:-1,1:-1] = (1-w)*Ac[1:-1,1:-1] + w*((A[2:,1:-1]+A[0:-2,1:-1])/dr**2 + nn[1:-1,1:-1]*(A[1:-1,2:]+A[1:-1,0:-2]) + (S1[1:-1,1:-1]+Snew[1:-1,1:-1]))/dnew[1:-1,1:-1]
        
    A = solve_boundary(A)

    error = L2_error(A, Ac)
    
    
    print(error)
    if error<1e-8:
        break 
    if math.isnan(error):
        print("Exiting loop. Nan found")
        break
        

stop = timeit.default_timer()

print('Time: ', stop - start) 

out_folder = 'Output'
if not os.path.exists(out_folder):
  os.makedirs(out_folder)

np.savetxt(out_folder+'A_nr_%s_%s.txt'%(str(strength),str(power)), A)


