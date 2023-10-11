import numpy as N
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

os.system('mkdir frames')

# customized colormap
mycolor = LinearSegmentedColormap.from_list('mycmap',['gray','black'])

# set image resolution for frames to be saved
yr = 10
xr = int(yr*16/9)
mpl.rcParams['figure.figsize'] = xr,yr

# function to determine whether coordinate belongs to MB set or not
def mandelbrot_set(xmin, xmax, ymin, ymax, width, height, max_iterations):
    x = N.linspace(xmin, xmax, width).reshape((1, width))
    y = N.linspace(ymin, ymax, height).reshape((height, 1))
    c = x + y * 1j
    z = N.zeros_like(c)
    for i in range(max_iterations):
        z = z**2 + c
    mandelbrot = (abs(z) < 2).astype(int)
    return mandelbrot

height = int(2160/2)            # resolution of image (does not need to match the resolution of the savec frames)
width  = int(16/9*height)       # aspect ratio = 16/9

Nmax = 25       # number of frames between the initial and final frame (i.e., sets the zoom speed)
zfac = 20/Nmax
itermax = 600   # this number depends on how deep the zoom is. Deeper means larger itermax
myiters = N.linspace(25,600,Nmax)

for ii in range(0,Nmax):
    
    print(ii,'of ',Nmax,' frames saved')

    # present frame zoom (XMIN, XMAX, etc.)
    XMIN = -1*16/9*2**(-ii*zfac) - .74877   # .74877 is the x value of the point we approach
    XMAX =  1*16/9*2**(-ii*zfac) - .74877
    YMIN = -1*2**(-ii*zfac)      - .065176  # .065176 is the y value of the point we approach
    YMAX =  1*2**(-ii*zfac) - .065176
    ext = [XMIN,XMAX,YMIN,YMAX] 
    
    # function that aims to provide a sufficient number of iterations at a specific zoom level
    iters = int(15+zfac*30/abs(XMIN-XMAX)**(.25/zfac)-0*4*N.log(abs(XMIN-XMAX)))

    iters = int(myiters[ii])
    mandelbrot = mandelbrot_set(ext[0],ext[1],ext[2],ext[3],width,height,iters)
    plt.imshow(mandelbrot, cmap=mycolor, extent=ext)
    plt.axis('off')
    plt.savefig('frames/fig'+f"{ii:06d}.png",bbox_inches='tight',pad_inches=0)
    plt.close("fig")
    plt.clf()
    gc.collect()
