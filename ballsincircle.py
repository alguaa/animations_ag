import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as N
from   scipy import interpolate
from   sympy import N as symN
from   sympy.solvers import solve
from   sympy import sqrt,atan,sin,cos,Symbol,re,im,diff
import scipy.linalg as SLA
import os

os.system('mkdir frames')

def mm(* args):
    tmp=N.dot(args[0],args[1])
    for ii in range(len(args)-2):
        tmp=N.dot(tmp,args[ii+2])
    return tmp

def myangle(xp1,yp1):
    if xp1<0 and yp1<0:
        ang = N.pi+N.arctan(yp1/xp1)
    elif xp1>0 and yp1<0:
        ang = N.arctan(yp1/xp1)
    return ang

NN=10000
mycircle = N.array([[1.01*N.cos(2*N.pi*ii/NN),\
                     1.01*N.sin(2*N.pi*ii/NN)] for ii in range(NN)]).T

#Parameters
g             = 9.81 #acceleration
startposition = -.5  #along x axis. 0 is center of circle
noballs       = 10   #number of balls
jumps         = 5    #number of bounces
dt            = .005 #time increment, should be small
eps           = 1e-2 #difference between adjacent balls


def findpaths(xp1,yp1,v0):
    for kk in range(jumps):
        angle = myangle(xp1,yp1)
        dvec = [-N.sin(angle),N.cos(angle)]
        nvec = [-dvec[1],dvec[0]]
        nvec = nvec/N.sqrt(mm(nvec,nvec))
        if kk == 0:
            ivec = [0,-1]
            ivec = N.array([float(ivec[0]),float(ivec[1])])
            ivec = ivec/SLA.norm(ivec)
        else:
            ivec = [1,slope]
            ivec = N.array([float(ivec[0]),float(ivec[1])])
            ivec = ivec/SLA.norm(ivec)
            if slope>0 and xp1<0:
                ivec[0] = -abs(ivec[0])
                ivec[1] = -abs(ivec[1])
            if slope<0 and xp1<0:
                ivec[0] = -abs(ivec[0])
                ivec[1] = abs(ivec[1])
            if slope>0 and xp1>0:
                ivec[0] = -abs(ivec[0])
                ivec[1] = -abs(ivec[1])
            if slope<0 and xp1>xp0:
                ivec[0] = abs(ivec[0])
                ivec[1] = -abs(ivec[1])
            if slope>0 and xp1>xp0:    
                ivec[0] = abs(ivec[0])
                ivec[1] = abs(ivec[1])
        rvec = (ivec-mm(2*mm(ivec,nvec),nvec))
        rvec = N.array([float(rvec[0]),float(rvec[1])])
        rvec = rvec/SLA.norm(rvec)
        v0 = v0*rvec
        vx = v0[0]
        vy0 = v0[1]
        x = Symbol('x')
        sol = solve(x**2 + (vy0/vx*(x-xp1)-g/2/vx**2*(x-xp1)**2+yp1)**2 - 1,x)
        for ii in range(len(sol)):
            if abs(re(sol[ii])-xp1)>1e-8 and abs(im(sol[ii]))<1e-4: 
                soln = re(sol[ii])
        xi,xe = xp1,soln
        xp0 = xp1
        yp0 = yp1
        pathx = N.arange(xi,xe,dt*vx)
        pathy = vy0/vx*(pathx-xi)-g/2/vx**2*(pathx-xi)**2+yp1
        xs = Symbol('xs')
        slope = diff(vy0/vx*(xs-xi)-g/2/vx**2*(xs-xi)**2+yp1).evalf(subs={xs:xe})
        yp1 = vy0/vx*(xe-xi)-g/2/vx**2*(xe-xi)**2+yp1
        yp1 = float(yp1)
        xp1 = float(xe)
        slope = float(slope)
        Etot = g*(1+y0)
        Ep = g*(1+yp1)
        v0 = (2*(Etot-Ep))**.5
        ballpaths.append([pathx,pathy])
    return ballpaths

ball = []
for ll in range(noballs):
    print('ball',ll)
    x0 = startposition+ll*eps 
    y0 = .0
    xp1 = x0
    yp1 = -N.sqrt(1-xp1**2)
    ds = N.sqrt((y0-yp1)**2)
    t0 = N.sqrt(2*ds/g)
    v0 = g*t0
    xtot0,ytot0 = [],[]
    tt = 0
    while dt*tt<t0:
        ytot0.append(y0-g*(dt*tt)**2/2)
        xtot0.append(x0)
        tt+=1
    xtot,ytot = N.array(xtot0),N.array(ytot0)
    ballpaths = []
    p = findpaths(xp1,yp1,v0)
    xtot = N.concatenate([xtot,p[0][0]])
    ytot = N.concatenate([ytot,p[0][1]])
    #ytot = p[0][1]
    for ii in range(1,jumps):
        xtot = N.concatenate([xtot,p[ii][0]])
        ytot = N.concatenate([ytot,p[ii][1]])
    ball.append([xtot,ytot])

print(len(ball))
for ii in range(len(xtot)):
    print(ii,'of',len(xtot))
    plt.plot(mycircle[0],mycircle[1])
    for bb in range(noballs):
        plt.plot(ball[bb][0][ii],ball[bb][1][ii],marker='.')
    plt.savefig('frames/fig'+f"{ii:05d}.png",bbox_inches='tight',pad_inches=0)
    plt.clf()
