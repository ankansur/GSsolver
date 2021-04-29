import numpy as np
from matplotlib import pyplot as plt
import sys, os
from scipy.interpolate import griddata
from scipy import special
from scipy.special import legendre
import glob
import matplotlib.cm as cm
import matplotlib as mpl

strength = 25
strengths = [0,5,10,15,20]
power = 1.1
rmin = 0.0

def spherical_to_cartesian2D(r, theta):
    
    #Given 1D arrays for r and theta, the function makes a spherical (r,theta)
    #grid and then transforms it to cartesian coordinates. It outputs two 2D
    #arrays, one for x and one for z.
    
    theta_matrix, radius_matrix = np.meshgrid(theta, r)
    x = radius_matrix*np.sin(theta_matrix)
    z = radius_matrix*np.cos(theta_matrix)
    return x, z

def spherical_to_cartesian_bfield_2D(N_r, N_theta, theta, B_r, B_theta):
    """
    Given a vector field (B_r,B_theta) in spherical coordinates, the function
    outputs a field in cartesian coordinates.
    """

    B_x = np.zeros((N_r, N_theta))
    B_z = np.zeros((N_r, N_theta))

    for i in range(0, N_theta):
        B_x[:,i] = B_r[:,i]*np.sin(theta[i]) + B_theta[:,i]*np.cos(theta[i])
        B_z[:,i] = B_r[:,i]*np.cos(theta[i]) - B_theta[:,i]*np.sin(theta[i])
    return B_x, B_z

def stream_function(r, theta, B_r, B_theta, order=1):
    """
    Given a spherical (r,theta) grid which can be nonuniform, and a
    divergenceless axisymmetric vector field B(r, theta), this function returns
    a scalar function psi(r, theta)called the Stokes stream function, such that
    B = rot(A), where A is vector potential given by
    A = psi(r,theta)/r*sin(theta) * e_phi. The contours of psi are the field
    lines of the original vector field B.

    Routine originally written for ASH data by ???
    Modified for PLUTO by S. Matt (2 Dec., 2011)
    Translated to Python by F. Bartolic(August, 2016)
    """

    N_r = len(r)
    N_theta = len(theta)

    #calculate derivatives of the stream function
    dpsi_dr = np.zeros((N_r,N_theta))
    dpsi_dtheta = np.zeros((N_r,N_theta))

    for i in range(0, N_theta):
        dpsi_dr[:, i] = -r*np.sin(theta[i])*B_theta[:, i]
        dpsi_dtheta[:, i] = r*r*np.sin(theta[i])*B_r[:,i]

    # start at integrating at pole at inner radius, do all along pole, then do
    # for each theta
    psi = np.zeros((N_r,N_theta))
    psi_2 = np.zeros((N_r, N_theta))
    dtheta = np.zeros(N_theta)
    dr = np.zeros(N_r)
    if order >= 0:
        dr[1:] = r[1:] - r[:-1]
        dtheta[1:] = theta[1:] - theta[:-1]
        psi[1:, 0] = psi[:-1, 0] + dpsi_dr[1:, 0]*dr[1:]

        for i in range(1, N_theta):
             psi[:, i] = psi[:, i-1] + dpsi_dtheta[:, i]*dtheta[i]

    if order <= 0:
        dr[:-1] = r[:-1] - r[1:]
        dtheta[:-1] = theta[:-1] - theta[1:]

        for i in range(N_r-2, -1, -1):
            psi_2[i, N_theta - 1] = psi_2[i + 1, N_theta - 1] +\
                                dpsi_dr[i, N_theta - 1]*dr[i]
        for i in range(N_theta-2, -1, -1):
            psi_2[:, i] = psi_2[:, i + 1] + dpsi_dtheta[:, i]*dtheta[i]

        if order < 0:
            return psi_2
        else:
            psi = 0.5*(psi + psi_2)  # Avg of the two
    return psi


def fieldlines(r, theta, B_r, B_th, magB):
    x, z = spherical_to_cartesian2D(r, theta)
    theta0 = np.linspace(0.0,np.pi,Nth)
    x_r = np.cos(theta0)
    y_r = np.sin(theta0)
    w = max(r)
    Z_l, X_l = np.mgrid[-w:w:100j, 0:w:100j]
    px_l = x.flatten()
    pz_l = z.flatten()
    px_r = -x.flatten()
    pz_r = z.flatten()
    f = plt.figure(figsize=(10,12))
    ax = f.add_subplot(111)
    B_x_l, B_z_l = spherical_to_cartesian_bfield_2D(Nr, Nth, theta, B_r, B_th)
    pu_l = B_x_l.flatten()
    pv_l = B_z_l.flatten()
    gu_l = griddata(zip(px_l,pz_l), pu_l, (X_l,Z_l))
    gv_l = griddata(zip(px_l,pz_l), pv_l, (X_l,Z_l))
    ax.streamplot(X_l,Z_l,gu_l,gv_l, color='blue', linewidth=1.5)
    Bmag = ax.contourf(x, z, magB, 50, cmap='Purples', extend='both')
    #cbar = plt.colorbar(Bmag,fraction=0.046, pad=0.02)
    #cbar.set_label(r"B$_{\rm tor}$ / B$_{\rm 0}$", size=25)
    #ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'x / R',size=25)
    ax.set_ylabel(r'z / R',size=25)
    #ticklabs = cbar.ax.get_yticklabels()
    #cbar.ax.set_yticklabels(ticklabs, fontsize=15)
    #cbar.ax.tick_params(labelsize=15)
    ax.set_xlim(0.0,1.5)
    ax.set_ylim(-1.5,1.5) 
    #ax.set_title(r't = %d'%(t),size=25)
    #plt.savefig('fig.png', bbox_inches='tight')
    plt.show()
    plt.close(f)
    
def Beta(A,strength,power,rid,uid):
    b = np.zeros((Nr,Nth))
    for i in range(Nr):
        for j in range(Nth):
            if A[i,j]>A[rid,uid]:
                b[i,j] = strength*(A[i,j]-A[rid,uid])**power  
            else:
                b[i,j] = 0.0 
                
                
    return b
    
def Beta2(A,strength,power):
    b = np.zeros((Nr,Nth))
    for i in range(Nr):
        for j in range(Nth):
            b[i][j] = strength**2*A[i][j] 
                
                
    return b
    
def plot_field(A,B,strength,power,lvs3):
    lvs1 = np.linspace(np.amin(A),np.amax(A),12)
    #print lvs1
    #lvs2 = [0,0.0048,0.01,0.016,0.018,0.02,0.0244]
    #Btor = np.sqrt(B)
    theta0 = np.linspace(0.0,2*np.pi,100)
    x_r = np.cos(theta0)
    y_r = np.sin(theta0)
    x, z = spherical_to_cartesian2D(r, theta)
    Psi = stream_function(r, theta, B_r, B_th)
    f = plt.figure(figsize=(10,12))
    ax = f.add_subplot(111)
    Psimin = np.amin(Psi)
    Bmag = ax.contourf(x,z, B, 200, cmap='GnBu', alpha=0.7, vmin=0.0)
    cbar = plt.colorbar(Bmag,fraction=0.15, pad=0.02)
    #Bfieldl = ax.contour(x, z, A, levels=lvs1, linewidths=1, origin='lower', colors='black',linestyles='solid') 
    Bfieldl = ax.contour(x, z, A, levels=lvs1, linewidths=1, origin='lower', colors='black',linestyles='solid') 
    Bfieldl = ax.contour(x, z, A, levels=lvs3, linewidths=1, origin='lower', colors='k',linestyles='solid') 
    cbar.set_label(r'$\beta$', size=25)
    cbar.ax.tick_params(labelsize='x-large')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.plot(y_r,x_r, c='red',linewidth=4,linestyle='--')
    #ax.plot(0.8*y_r,0.8*x_r, c='g',linewidth=0.5,linestyle='--')
    ax.set_xlabel(r'x / R',size=25)
    ax.set_ylabel(r'z / R',size=25)
    #ax.text(1.0,0.8,r'$s=%d$'%strength,fontsize=20)
    #ax.text(1.0,0.9,r'$p=%s$'%str(power),fontsize=20)
    ax.set_xlim(0.0,1.1)
    ax.set_ylim(-1,1)
    #plt.savefig('An1_field_strength%d.pdf'%strength, bbox_inches='tight')
    plt.savefig('A_n1_25_1.1.pdf', bbox_inches='tight')
    plt.show()
    #plt.savefig('Pol2_nr.pdf',bbox_inches='tight')
    plt.close(f)
    

def Tor_energy(br, bth, bphi, r, theta):
    idx=-1
    tor_e = 0
    pol_e = 0
    dr = r[2]-r[1]
    dth = theta[1]-theta[2]
    #br = br[0:idx,:]
    #bth = bth[0:idx,:]
    #bphi = bphi[0:idx,:]
    nr = len(br)
    nth = len(br[0])
    for i in range(nr):
        for j in range(nth):
            tor_e += bphi[i,j]**2*r[i]**2*np.sin(theta[j])*dr*dth
            #tor_e += bphi[i,j]**2*dr*du
            pol_e += (br[i,j]**2 + bth[i,j]**2)*r[i]**2*np.sin(theta[j])*dr*dth
            #pol_e += (br[i,j]**2 + bth[i,j]**2)*dr*du
    return pol_e,tor_e
        
        
def Bpoleavg(br, bth, bphi, r, theta):
    idx=-1
    pol_B = 0
    avg_B = 0
    dr = r[2]-r[1]
    dth = theta[1]-theta[2]
    #br = br[0:idx,:]
    #bth = bth[0:idx,:]
    #bphi = bphi[0:idx,:]
    nr = len(br)
    nth = len(br[0])
    vol = 0
    for i in range(nr):
        for j in range(nth):
            avg_B += (br[i,j]**2 + bth[i,j]**2 + bphi[i,j]**2)*r[i]**2*np.sin(theta[j])*dr*dth
            vol+=r[i]**2*np.sin(theta[j])*dr*dth
            
    for i in range(rid-6,rid+6,1):
        for j in range(nth):
            pol_B += (br[i,j]**2 + bth[i,j]**2 + bphi[i,j]**2)*r[i]**2*np.sin(theta[j])*dr*dth
    return 2*pol_B/avg_B

A = np.loadtxt('A_n1_%s_%s.txt'%(str(strength),str(power)))
u = np.linspace(-1,1,len(A[0]))
r = np.linspace(rmin,2.0,len(A))
theta = np.arccos(u)
Nr = len(r)
Nth= len(theta)
rid = np.where(r>1)[0][0]-1
uid = np.where(u>0)[0][0]-1       
dr = r[2]-r[1]
dth = theta[1]-theta[2]
du = u[1]-u[0]
B_r = np.zeros((Nr,Nth))
B_th = np.zeros((Nr,Nth))
B_phi = np.zeros((Nr,Nth))

#A = np.loadtxt('A_n1_%s_%s.txt'%(strength,power))
#A = np.loadtxt('A_n1_mixed10_1.1.txt')
#A = np.loadtxt('test/A_nr2_0_1.1.txt')
B = Beta(A,strength,power,rid,uid) 
#np.savetxt('B_nr_%s_%s.txt'%(strength,power),B)
lvs3 = np.amax(A)  
plot_field(A,B,strength,power,lvs3)

for i in range(1,Nr-1):
    for j in range(1,Nth-1):
        #Solve each component
        #r component
        B_r[i][j] = 1.0/r[i]/r[i]/np.sin(theta[j])*(A[i][j+1]-A[i][j-1])/2/dth
        #th component
        B_th[i][j] = -1.0/r[i]/np.sin(theta[j])*(A[i+1][j]-A[i-1][j])/2/dr
        #phi component
        B_phi[i][j] = 1.0/r[i]/np.sin(theta[j])*(B[i][j])


poloidal_energy,toroidal_energy = Tor_energy(B_r, B_th, B_phi, r, theta)
#print toroidal_energy/(poloidal_energy+toroidal_energy)*100
#print np.sqrt(poloidal_energy+toroidal_energy)
print Bpoleavg(B_r, B_th, B_phi, r, theta)
#print toroidal_energy/(poloidal_energy)*100