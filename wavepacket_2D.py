import numpy as N
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import sparse
import scipy.linalg as SLA
from scipy.sparse.linalg import isolve
import time
import os

os.system('mkdir frames')

# main program PARAMS
Nx   = 200           
Ny   = Nx           # only tried for Nx=Ny
tmax = 400          # number of frames
xmax = 1.;ymax = 1.
dx   = 1./(Nx-1)
a    = 1j/2/dx**2   # time term: Adds to diagonal of FD Hamiltonian
b    = 1/4/dx**2.   # off-diagonal elements of FD Hamiltonian
TOL  = 1e-6         # tolerance in iterative linear solver

# initial Gaussian and its momentum
sig  = .07                # width of Gaussian
cwp  = [Nx/3,Ny/4]        # start position
pwp  = [-4*Nx/10,6*Nx/10] # momentum (2pi/lambda)
    
# potential
V = N.zeros((Nx,Ny))
# define whatever potential you wish here

# Gaussian wave packet
psi = [N.exp(1j*(pwp[0]*(ii-cwp[0])+pwp[1]*(jj-cwp[1]))*dx-\
             ((ii-cwp[0])**2+(jj-cwp[1])**2)*dx**2/(2*sig**2))\
       for ii in range(Nx) for jj in range(Ny)]
psi = N.array(psi).reshape(Nx,Ny)
psi = psi/SLA.norm(psi)

# create FD Hamiltonian for one slice of the domain
lhs = N.zeros((3,Nx))
lhs[0,:] = b*N.ones(Nx) # upper off-diagonal
lhs[2,:] = lhs[0,::-1]  # lower off-diagonal

# returns FD Ham in sparse csr format to be used in an iterative solver below
# only the diagonal is updated in LHS/LHS2

def LHS(idx):
    mtx = sparse.spdiags([lhs[0],2*a-2*b-V[:,idx]/2,lhs[2]],[1,0,-1],Nx,Nx)
    return sparse.csr_matrix(mtx) 
def LHS2(idx):
    mtx = sparse.spdiags([lhs[0],2*a-2*b-V[idx,:]/2,lhs[2]],[1,0,-1],Nx,Nx)
    return sparse.csr_matrix(mtx)
def RHS(idx):
    return  N.array((2*a+2*b+V[:,idx]/2)*psi[:,idx]-\
                    b*N.array([psi[ii,min(Nx-1,idx+1)]+psi[ii,max(0,idx-1)] for ii in range(Nx)]))
def RHS2(idx):
    return N.array((2*a+2*b+V[idx,:]/2)*psi[:,idx]-\
                   b*N.array([psi[jj,min(Nx-1,idx+1)]+psi[jj,max(0,idx-1)] for jj in range(Nx)]))

def plotandsave(tt,psi):
    plt.imshow(N.abs(psi)**2)
    # add possible line plots, etc. here
    plt.savefig('frames/fig'+f"{tt:04d}.png") 

t1 = time.process_time()
for tt in range(tmax):
    plotandsave(tt,psi)
    psi = N.array([isolve.gmres(LHS(jj),RHS(jj),tol=TOL)[0] for jj in range(Nx)])
    psi = N.array([isolve.gmres(LHS2(ii),RHS2(ii),tol=TOL)[0] for ii in range(Nx)])
    psi = psi/SLA.norm(psi)
    if N.mod(tt,round(tmax/10)) == 0:
        print(round(100*tt/tmax),'% of the frames saved')
tm = time.process_time()-t1
print('\nAll frames saved in ',round(tm,1),'s')
