# ag_animations

Short python scripts that generate frames to be used to animate specific concepts in physics and mathematics.

All scripts run on python3 and need numpy, matplotlib, and possibly other libraries evident from what a specific script tries to import.

The scripts:


    'fourier' generates time dependent plots on epicycles (connected rotating rods) where the apex of the outermost rod draws an approximative path of a predefined contour in the complex plane. Here the contour is exemplified by a square. The contour may have any shape.


    'wavepacket_1D' generates time dependent plots of a wave packet moving in one dimension. The script exemplifies a Gaussian wave packet scattering off a rectangular potential, qualitatively demonstrating quantum tunneling.


    'wavepacket_2D' generates time dependent plots of a wave packet moving in two dimensions. The script exemplifies a Gaussian wave packet moving in a rectangular domain having imprenetable walls.


    'ballsincircle' generates time dependent plots of balls dropped on a circle's interior under gravity. The SymPy library is here used to compute exact reflection angles.


FFmpeg may be used to generate a video of the frames by a simple command line in the Terminal by e.g. executing the command while being in a frames folder: 
       ffmpeg -r 60 -f image2 -pattern_type glob -i "*.png"  -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" tmp.mov