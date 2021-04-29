#!/usr/bin/env python


#   ==========================================================================================================
#   Description:  A program to calculate the Grad-Shafranov equation to find equilibria models
#
#   Version: 5.0
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
strength = 10
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
                S[i,j] = (r[i])**2.0*(1-u[j]**2)
            else:
                S[i,j] = 0.0
    return S
    

  
def solve_boundary(A):
    for i in range(Nr):
        A[i,0] = 0
        A[i,Nu-1] = 0
    
    for j in range(Nu):
        A[0,j] = 0.0
        A[Nr-1,j] = 0.0#(1-u[j]**2)/r[Nr-1]/2.0#(5*r[0]**3/r[Nr-1] - 3*r[0]**5/r[Nr-1] + 3*r[Nr-1]**4 - 5*r[Nr-1]**2)*(1-u[j]**2)/30.0
           
    return A
    
    
def cal_B_mag():
    
    
    
A = np.zeros((Nr,Nu))
A = init(A)
A = solve_boundary(A)

S = np.zeros((Nr,Nu))
S = source(S)


Aaux = np.zeros((Nr,Nu))
def L2_error(p, pn):
    return np.sqrt(np.sum((p - pn)**2)/np.sum(pn**2))


start = timeit.default_timer()

w = 0.1

for n in range(tNum):

    Ac = A.copy()
               
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if Ac[i][j]>Ac[rid,uid]:
                S2 = strength**2*power*(Ac[i,j]-Ac[rid,uid])**(2.0*power-1.0)  
                dS2 = strength**2*power*(2*power-1)*(Ac[i][j] - Ac[rid,uid])**(2.0*power-2.0)
            else:
                S2 = 0.0
                dS2 = 0.0
            Q = (S2-dS2*Ac[i,j]) + max(0,dS2)*Ac[i,j]
            den = (2.0/dr**2 + 2.0*(1-u[j]**2)/r[i]**2/du**2) - min(0,dS2)
            A[i][j] = (1-w)*Ac[i,j] + w*((A[i+1][j]+A[i-1][j])/dr**2 + (1-u[j]**2)*(A[i][j+1]+A[i][j-1])/r[i]**2/du**2 + Q + S[i][j])/den
        
    #A = solve_boundary(A)

    error = L2_error(A, Ac)
    
    print(error)
    if error<1e-6:
        break 
    if math.isnan(error):
        print("Exiting loop. Nan found")
        break
        

stop = timeit.default_timer()

print('Time: ', stop - start) 

out_folder = '1Nov/'
if not os.path.exists(out_folder):
  os.makedirs(out_folder)

np.savetxt(out_folder+'A_n1_%s_%s.txt'%(str(strength),str(power)), A)