# animations_ag

- [Introduction](#introduction)
- [Usage](#usage)
- [Description](#description)

## Introduction
Short easy-to-use python scripts that generate frames for animating  specific concepts in physics and mathematics.

## Usage
To run a script, open Terminal, stand in the directory of the script and type:

       python (or python3) <script>.py

All scripts run on python3 and need numpy, matplotlib, and typically other libraries evident from what a specific script tries to import.

Execution of a script generates a 'frames' folder located at the same path as the script. If not (e.g., if executed through Visual Studio), check your home directory.

FFmpeg may be used to generate a video of the frames by a simple command line in the Terminal by executing the command below while being in the folder containg the frames:

       ffmpeg -pattern_type glob -i "*.png"  -pix_fmt yuv420p movie.mov


## Description

### fouriercontour.py  
Generates time dependent plots on epicycles (connected rotating rods) where the apex of the outermost rod draws an approximative path of a predefined contour in the complex plane. Here the contour is exemplified by a square. The contour may have any shape.

### wavepacket_1D.py 
Generates time dependent plots of a wave packet moving in one dimension. The script exemplifies a Gaussian wave packet scattering off a rectangular potential, qualitatively demonstrating quantum tunneling.

### wavepacket_2D.py  
Generates time dependent plots of a wave packet moving in two dimensions. The script exemplifies a Gaussian wave packet moving in a rectangular domain having impenetrable walls.

### ballsincircle.py  
Generates time dependent plots of balls dropped on a circle's interior under gravity. The SymPy library is here used to compute exact reflection angles. This script is particularly badly coded, but at least it works.

### mandelbrotzoom.py  
Generates plots of gradually increasing zoom around a particular coordinate of the Mandelbrot set. The script exemplifies zooming into the so-called seahorse valley.

### stringart.py
Generates a single stringart image based on an input .jpg image. Just type

	  python (or python3) stringart.py -img </path/to/image>.jpg

to obtain a first draft of your generated stringart. The default values seem working reasonably well on portraits; however further elaboration with the flags is likely needed. The script auto crops the input image based on its smallest side; hence, to properly frame your image, a substantially quadratic input image is recommended. Needs the PIL library (pip install Pillow). For further usage information, type

   	  python (or python3) stringart.py --h

### softbody.py
Generates frames of simple rectangularly arranged grid points connected by springs falling under gravity and bouncing on a floor. Aimed to be a starting point sandbox. Just execute by

	  python (or python3) softbody.py

and elaborate with the parameters in the script.