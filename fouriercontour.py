import numpy as N
import matplotlib.pyplot as plt
import os,time

os.system('mkdir ./frames') #Creates a folder in which the frames are saved

timesteps = 61
dt = 1/timesteps    #incremental time between subsequent frames
factor = 3;         #number of periods
tsRange = timesteps*factor

nn = 1000000000    #number of points defining a function when made here
mm = 7       #number of basis functions in the F-expansion
lw = 1       #line width of plots

#Define function in complex plane
function2=[]
for ii in range(nn):
    if ii<int(nn/4)+1:
        function2.append(1+1j*4*ii/nn)
    if ii>int(nn/4) and ii<int(2*nn/4)+1:
        function2.append(1+4*(ii-int(nn/4))/nn+1j)
    if ii >int(2*nn/4) and ii<int(3*nn/4)+1:
        function2.append(2+1j-1j*4*(ii-int(2*nn/4))/nn)
    if ii>int(3*nn/4):
        function2.append(2-4*(ii-int(3*nn/4))/nn)
function2 = N.array(function2)+.1*1j

#Returns positions for all rods for each time stamp
def myfourier(function,nn,mm):
    cf = [1./nn*N.sum([function[ii]*N.exp(-2*N.pi*1j*ii*jj/nn)\
                       for ii in range(0,nn)]) for jj in range(-mm,mm+1)]
    fourier = [N.sum([cf[ii]*N.exp(2*N.pi*1j*(ii-mm)*tt*dt)\
                      for ii in range(0,2*mm+1)]) for tt in range(0,timesteps+1)]
    #sort frequencies to match coefficients
    mmAsUsed = [ii for ii in range(-mm,mm+1)];
    mmSorted = sorted(mmAsUsed,key=abs)
    idx = [(N.abs(mmSorted[ii]-N.array(mmAsUsed))).argmin() for ii in range(2*mm+1)]
    apprY = [];Ytvec = [];barvec = []
    for ttmp in range(tsRange):
        bars = []
        bars.append(cf[idx[0]])
        for ii in range(1,2*mm+1):
            bars.append(bars[ii-1]+cf[idx[ii]]*N.exp(2*N.pi*1j*(idx[ii]-mm)*ttmp*dt))
        barvec.append(bars)
        apprY.append(N.imag(bars[ii]))
        Ytvec.append(N.real(bars[ii]))    
    return barvec,apprY,Ytvec
t1 = time.process_time()
fig2 = myfourier(N.array(function2),nn,mm)
tm = time.process_time()-t1
print('\nFourier expansion done in:',round(tm,2),'s')


for ttmp in range(tsRange):
    plt.plot(N.array(fig2[2][:ttmp+1]),N.array(fig2[1][:ttmp+1]),color=[0,1,0],linewidth=lw)
    for ii in range(1,2*mm+1):
        plt.plot(N.real(fig2[0][ttmp][ii-1:ii+1]),N.imag(fig2[0][ttmp][ii-1:ii+1]),color='magenta',linewidth=lw)
    plt.xlim([.7,2.3])
    plt.ylim([-.1,1.3])
    plt.savefig('./frames/fig'+f"{ttmp:05d}.png")
    plt.close()
    if N.mod(ttmp,round(tsRange/10)) == 0:
        print(round(100*ttmp/tsRange),'% of the frames saved')
print('All',tsRange,' frames saved')
