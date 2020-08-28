#!/usr/bin/env python


#   ==========================================================================================================
#   Description:  A program to calculate the Grad-Shafranov equation to find equilibria models
#
#   Version: 2.0
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
strength = 0
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
    
def source2(Ac,rid,uid):
    SS = np.zeros((Nr,Nu))
    for i in range(1,Nr-1):
        for j in range(1,Nu-1):
            if Ac[i,j]>Ac[rid,uid]:
                SS[i,j] = strength**2*power*(Ac[i,j]-Ac[rid,uid])**(2.0*power-1.0)  
            else:
                SS[i,j] = 0.0    
    return SS
    
    
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
errors = [1e-4,1e-5,1e-6,1e-7,1e-8, 1e-9]
l=0
itr = []
takes =[]

start = timeit.default_timer()

for n in range(tNum):

    #error = np.zeros((Nr,Nu))
    Ac = A.copy()
    
    ## d1 = r**2*du**2, d2 = dr**2*(1-u**2),  d3 = 2(r**2*du**2 + dr**2*(1-u**2))
       
    S2 = source2(Ac,rid,uid)
        
    A[1:-1,1:-1] = (d1[1:-1,1:-1]*(Ac[2:,1:-1]+Ac[0:-2,1:-1])+d2[1:-1,1:-1]*(Ac[1:-1,2:]+Ac[1:-1,0:-2])+dr**2*d1[1:-1,1:-1]*(S[1:-1,1:-1]+S2[1:-1,1:-1]))/d3[1:-1,1:-1]
        
    A = solve_boundary(A)

    error = L2_error(A, Ac)
    
    if error<errors[l]:
        l = l+1
        itr.append(n)
    
    print(error)
    if error<1e-8:
        break 
    if math.isnan(error):
        print("Exiting loop. Nan found")
        break
        

stop = timeit.default_timer()

print('Time: ', stop - start)

print itr 

def plot2D(x, y, p, case):

    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    if case==1:
        pyplot.savefig('SimA.pdf')
    else:
        pyplot.savefig('Analy.pdf')
    pyplot.show()

np.savetxt('A_n1_%s_%s.txt'%(str(strength),str(power)), A)
#np.savetxt('A_n1_%s_%s.txt'%(str(strength),str(power)), A)
#np.savetxt('A_n1.txt', A)