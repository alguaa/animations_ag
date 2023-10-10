# ag_animations
<style>
.gray-text {
    color: gray;
}
</style>

Short python scripts that generate frames to be used to animate specific concepts in physics and mathematics.

To run a script, open Terminal and stand in the directory of the script and type:

       python [or python3] <script>.py

All scripts run on python3 and need numpy, matplotlib, and typically other libraries evident from what a specific script tries to import.

Execution of a script generates a 'frames' folder located at the same path as the script. If not (e.g., if executed through Visual Studio), check your home directory.

FFmpeg may be used to generate a video of the frames by a simple command line in the Terminal by executing the command below while being in the folder containg the frames:

       ffmpeg -pattern_type glob -i "*.png"  -pix_fmt yuv420p movie.mov


*** Short description of the scripts ***

**fouriercontour.py** generates time dependent plots on epicycles (connected rotating rods) where the apex of the outermost rod draws an approximative path of a predefined contour in the complex plane. Here the contour is exemplified by a square. The contour may have any shape.

`<span class="gray-text">wavepacket_1D.py</span>` generates time dependent plots of a wave packet moving in one dimension. The script exemplifies a Gaussian wave packet scattering off a rectangular potential, qualitatively demonstrating quantum tunneling.

wavepacket_2D.py generates time dependent plots of a wave packet moving in two dimensions. The script exemplifies a Gaussian wave packet moving in a rectangular domain having imprenetable walls.

ballsincircle.py generates time dependent plots of balls dropped on a circle's interior under gravity. The SymPy library is here used to compute exact reflection angles. This script is particularly badly coded, but at least it works.