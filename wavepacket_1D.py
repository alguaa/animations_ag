import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as N   
import scipy.linalg as SLA
import os

os.system('mkdir frames')

mpl.rcParams['figure.figsize'] = 16/9*6,6

ymin = -2
xmin = ymin*16/9
xmax = abs(xmin);ymax = abs(ymin)

NN   = 2000        # number of x coordinates
pwp  = 40          # momentum
sig  = .15         # width of wave packet
Vmax = pwp**2/2*1. # height of tunneling barrier
noframes = 200     # number of frames

dx = (2*xmax-2*xmin)/NN
xvec = N.arange(xmin*2,xmax*2,dx)
b  = -1/2/dx**2

# finite difference Laplacian 
diag   = -N.diagflat(2*b*N.ones(NN))
offdiag = N.diagflat(b*N.ones(NN-1),1)
D2 = diag+offdiag+offdiag.transpose()

# add potential to diagonal
pot = N.zeros(NN)
for ii in range(NN):
    if (ii>int(.495*NN) and ii<int(.505*NN)):
        pot[ii] = Vmax

# Hamiltonian; H = D2 + V
H = D2+N.diagflat(pot)

ev,evec = SLA.eigh(H)
idx = ev.argsort()
ev,evec = ev[idx],evec[idx]

initwf = N.exp(-(xvec+2)**2/sig-1j*pwp*(xvec+2)) #initial Gaussian
init = ev[0]
cm = N.array([N.sum(initwf*evec[:,mm]) for mm in range(len(ev))])
cm = cm/SLA.norm(cm)
psi = [cm[jj]*evec[:,jj] for jj in range(len(ev))]

fig = plt.figure()
for tt in range(noframes):
    tt2 = tt/noframes
    Psi = psi[0]*N.exp(1j*ev[0]*tt2)

    for jj in range(1,len(psi)):
        Psi+=psi[jj]*N.exp(1j*ev[jj]*tt2)

    plt.xlim([xmin,xmax])
    plt.ylim([ymin*.3,ymax*.3])
    plt.gca().add_patch(Rectangle((-4,-10),10,10,linewidth=1,edgecolor=None,facecolor='gray'))

    if max(pot)>0:
        plt.fill(xvec,pot/abs(b)*.5,linewidth=1,color='gray')

    else:
        plt.fill(xvec,pot/abs(b)*.5,linewidth=1,color='black')
        
    plt.plot(xvec,5*N.abs(Psi)**2,linewidth=1,color='red')
    plt.savefig('frames/fig'+f"{tt:04d}.png",
                bbox_inches='tight',pad_inches=0)
    plt.clf()
    plt.close("fig")
