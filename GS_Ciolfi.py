#!/usr/bin/env python


#   ==========================================================================================================
#   Description:  A program to calculate the Grad-Shafranov equation to find equilibria models
#
#   Version: 4.0
#   Date: June 2020
#   

#   Equation: Ciolfi Toroidal
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
u = np.linspace(-1,1,51)
r = np.linspace(0.0,2.0,51)

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
k = 0.1

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
    
    
def source2(Ac,rid,uid):
    SS = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if np.abs(Ac[i,j]/Ac[rid,uid])>1:
                SS[i,j] = strength*Ac[i,j]*(np.abs(Ac[i,j]/Ac[rid,uid])-1)
            else:
                SS[i,j] = 0.0
    return SS
    
    
def Falpha(Ac, rid, uid):
    ff = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if np.abs(Ac[i,j]/Ac[rid,uid])>k: 
                ff[i,j] = r[i]**2.0*(1-u[j]**2)*((1.0-np.abs(Ac[i,j]/Ac[rid,uid]))**4-k)
            else:
                ff[i,j] = 0.0
    return ff
   
    
def den1(d1):
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            d1[i,j] = r[i]**2*du**2
    return d1
    
    
    
def den2(d2):
    for i in range(Nr):
        for j in range(Nu):
            d2[i,j] = dr**2*(1-u[j]**2)    
    return d2
      
    
def solve_boundary(A):
    for i in range(Nr):
        A[i,0] = 0
        A[i,Nu-1] = 0
    
    for j in range(Nu):
        A[0,j] = 0.0
        A[Nr-1,j] = 0.0#(1-u[j]**2)/r[Nr-1]/2.0#(5*r[0]**3/r[Nr-1] - 3*r[0]**5/r[Nr-1] + 3*r[Nr-1]**4 - 5*r[Nr-1]**2)*(1-u[j]**2)/30.0
           
    return A
    
A = np.zeros((Nr,Nu))
A = init(A)
A = solve_boundary(A)

S = np.zeros((Nr,Nu))
S = source(S)

d2 = np.zeros((Nr,Nu))
d2 = den2(d2)
d1 = np.zeros((Nr,Nu))
d1 = den1(d1)
d3 = np.zeros((Nr,Nu))

Aaux = np.zeros((Nr,Nu))

def L2_error(p, pn):
    return np.sqrt(np.sum((p - pn)**2)/np.sum(pn**2))

d3 = (d1+d2)*2 

start = timeit.default_timer()

w = 0.1

for n in range(tNum):

    Ac = A.copy()
    F = Falpha(Ac,rid,uid)
               
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if Ac[i][j]>Ac[rid,uid]:
                S2 = strength*Ac[i,j]*(np.abs(Ac[i,j]/Ac[rid,uid])-1) 
                dS2 = strength*((2*Ac[i,j]/Ac[rid,uid])-1)
                #dS = 1.1*strength**2*1.2
                dS1 = 0
                
            else:
                dS1 = -4*strength*(1-Ac[i,j]/Ac[rid,uid])**3/Ac[rid,uid]
                S2 = 0.0
                dS2 = 0.0
            dS = dS1+dS2
            Q = (S2+S[i][j]-dS*Ac[i,j]) + max(0,dS)*Ac[i,j]
            den = (2.0/dr**2 + 2.0*(1-u[j]**2)/r[i]**2/du**2) - min(0,dS)
            A[i][j] = (1-w)*Ac[i,j] + w*((A[i+1][j]+A[i-1][j])/dr**2 + (1-u[j]**2)*(A[i][j+1]+A[i][j-1])/r[i]**2/du**2 + Q + (F[i][j]+S2*dS2)*S[i][j])/den
        
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

out_folder = 'Output/'
if not os.path.exists(out_folder):
  os.makedirs(out_folder)

np.savetxt(out_folder+'A_C_%s_%s.txt'%(str(strength),str(power)), A)
#np.savetxt('A_nr_%s_%s.txt'%(str(strength),str(power)), A)
#np.savetxt('A_n1.txt', A)