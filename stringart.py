#!/usr/bin/env python3

import numpy as N
import matplotlib.pyplot as plt
import matplotlib as mpl
import random, time, sys, os
import argparse
from PIL import Image, ImageOps, ImageFilter, ImageEnhance

T, F = True, False

class stringart:
	def __init__(self,scalefactor,weight,lines,linewidth,\
		invertimage,contour,autosaturation,resolution,image_path):
		
		self.scalefactor = scalefactor
		self.weight = weight
		self.lines = lines
		self.linewidth = linewidth
		self.invertimage = invertimage
		self.contour = contour
		self.autosaturation = autosaturation
		self.resolution = resolution
		self.image_path = image_path
		
	def prepareimage(self):
		image = Image.open(self.image_path)
		image = ImageOps.grayscale(image)
		image = N.array(image)
		image = image/N.max(image)
		inputimage = image
		if self.invertimage == 'True':
			image = 1 - image
						
		dimension = N.min([N.shape(image)[0],N.shape(image)[1]])		
		importedimage = image[:dimension,:dimension]

		if self.scalefactor == None:
			sf = int(dimension/90)
		else:
			sf = self.scalefactor
		#sf = self.scalefactor
		newdim = int(dimension/sf)
		sampledimage = N.zeros((newdim,newdim))
		for ii in range(newdim):
			for jj in range(newdim):
				sampledimage[ii,jj] = N.mean(importedimage[sf*ii:sf*(ii+1),sf*jj:sf*(jj+1)])

		self.sf = sf
		image = sampledimage
		dimension = newdim		

		return image, dimension, inputimage

	def circumferencecoordinates(self,image,dim):
		circumference = N.zeros((dim,dim))
		circumference_coords = []
		# circle
		if self.contour == 1:
			for ii in range(dim):
				for jj in range(dim):
					circ1 = (ii-int(dim/2))**2+(jj-int(dim/2))**2
					circ2 = (ii-int(dim/2))**2+(jj-int(dim/2))**2 + dim

					if circ1 >= int(dim/2)**2:
						image[ii,jj] = 0
						
					if circ1 < int(dim/2)**2 and circ2 >= int(dim/2)**2:
						circumference[ii,jj] = 1
						circumference_coords.append([ii,jj])
		# rectangle
		if self.contour == 2:
			circumference_coords = []
			for ii in range(dim):
				for jj in range(dim):
					if ii == 0 or jj == 0 or ii == dim-1 or jj == dim-1:
						circumference_coords.append([ii,jj])
		
		return circumference_coords

	# apply one-pixel padding (for computational purposes)
	def pad_matrix(self,matrix):
		original_height, original_width = matrix.shape
		padded_matrix = N.zeros((original_height + 2, original_width + 2), dtype=float)
		padded_matrix[1:original_height + 1, 1:original_width + 1] = matrix

		return padded_matrix

	def getlinecoordinates(self, p1, p2, dim):
		highlight_color = 1
		mat = N.zeros((dim,dim))
		x1, y1 = p1
		x2, y2 = p2
		dx = abs(x2 - x1)
		dy = -abs(y2 - y1)
		sx = 1 if x1 < x2 else -1
		sy = 1 if y1 < y2 else -1
		error = dx + dy

		while True:
			mat[x1, y1] = highlight_color  # Color the element
			if x1 == x2 and y1 == y2:
			    break
			e2 = 2 * error
			if e2 >= dy:
				error += dy
				x1 += sx
			if e2 <= dx:
				error += dx
				y1 += sy
		mat = mat.T

		return mat


	def choosepair(self,idx, jjpair, pairs_list, jjlinemat):
		kk = 0
		while jjpair[idx[kk]] in pairs_list:
			kk += 1
			pass

		return kk

	def calculatestrings(self):
		image,dim = self.prepareimage()[:2]
		circumference_coords = self.circumferencecoordinates(image,dim)
		image = self.pad_matrix(image)
		if self.scalefactor == None:
			print(f'Calculated scalefactor:		{self.sf}')	
			print(f'Dimensions, prepared image: 	{dim}x{dim} pixels (set -SF flag manually for other dimensions)')
		else:
			print(f'Dimensions, prepared image: 	{dim}x{dim} pixels')
		dim = dim + 2
		nocoords = len(circumference_coords) # number of coords on circumference
		rint = random.randint(0,nocoords-1)  # random starting coordinate
		print(f'Circumference elements:		{nocoords} coordinates (recommended 250-500)\n')
		
		linemat = N.zeros((dim,dim))
		linecoords = []
		linecoordsmat = []
		diff = []
		pairs = []
		t1 = time.process_time()

		for ii in range(self.lines):
			magnitudevec = [] # how much the presently lined image differs from image
			jjlinemat = []     
			jjidx = []
			jjpair = []	
			pi = circumference_coords[rint]

			for jj in range(nocoords):
				pf = circumference_coords[jj]
				if N.sqrt((pi[0]-pf[0])**2 + (pi[1]-pf[1])**2) > 1:
					tmplinemat = linemat + self.weight*self.getlinecoordinates(pi,pf,dim)
					magnitudevec.append(N.sqrt(N.mean((image-tmplinemat)**2)))
					jjlinemat.append(tmplinemat)
					jjpair.append([pi,pf])
					jjidx.append(jj)

			magnitudevec = N.array(magnitudevec)
			idx = N.argsort(magnitudevec)
			magnitudevec = magnitudevec[idx]
			
			if ii == 0:
				pair = jjpair[idx[0]]
				pairs.append(pair)
			else:
				kk = self.choosepair(idx,jjpair,pairs,jjlinemat)
				pairs.append(jjpair[idx[kk]])
				linemat = jjlinemat[idx[kk]]
				linecoords.append(jjpair[idx[kk]])
				rint = jjidx[idx[kk]]
				
			diff.append(N.sqrt(N.mean((image-linemat)**2)))
			
			if self.autosaturation == 'True':
				if len(diff) > 200:
					if abs(diff[100] - diff[0]) > 100*abs(diff[-1] - diff[-100]):
						print('autosaturation chosen (-AS True); generation stopped.')
						break

			if ii % int(self.lines/25) == 0:
				print(f'{round(ii/(self.lines)*100)}% done ({ii} of {self.lines} lines), diff: {round(magnitudevec[idx[0]],3)}')

		return linecoords, diff, image



import argparse

flags = [
    {'name': 'scalefactor', 'short': '-sf', 'long': '--scalefactor', 'type': int, 'metavar': '', 'help': 'Downsampling factor of the original image. Recommended sampling should be such that the image is just above 100x100 pixels. If not specified, calculation will proceed with an estimated value.', 'default': None},
    {'name': 'weight', 'short': '-wt', 'long': '--weight', 'type': float, 'metavar': '', 'help': 'Parameter used upon minimization between lines and image. Elaboration with values between 0.02 and 0.5 is recommended. Larger weight = faster saturation (with possibly worse result).', 'default': 0.02},
    {'name': 'lines', 'short': '-ls', 'long': '--lines', 'type': int, 'metavar': '', 'help': 'Maximum number of lines used to generate the string art.', 'default': 3000},
    {'name': 'linewidth', 'short': '-lw', 'long': '--linewidth', 'type': float, 'metavar': '', 'help': 'String linewidth.', 'default': 0.1},
    {'name': 'invertimage', 'short': '-ii', 'long': '--invertimage', 'metavar': '', 'help': 'Invert image; True or False', 'default': 'True'},
    {'name': 'contour', 'short': '-co', 'long': '--contour', 'type': int, 'metavar': '', 'help': 'Choose between circle (1) or square (2) as contour.', 'default': 1},
    {'name': 'autosaturation', 'short': '-as', 'long': '--autosaturation', 'metavar': '', 'help': 'Auto finish string art upon reaching calculated saturation.', 'default': 'True'},
    {'name': 'resolution', 'short': '-res', 'long': '--resolution', 'type': str, 'metavar': '', 'help': 'Approximate resolution of the final image. Choose between 4k (~2160x2160), HD (~1080x1080), and VGA (~540x540)', 'default': 'HD'},
    {'name': 'image_path', 'short': '-img', 'long': '--image_path', 'type': str, 'metavar': '', 'help': 'Path to the image file. If the script and image are in the same folder, just type the filename + extension.', 'default': ''}
]

parser = argparse.ArgumentParser(description="The script aims to draw a string art representation of an imported image.\
    Elaboration with the flag parameters may likely be needed depending on your input image. The default flag settings \
    appear performing decently on portraits.")


for flag in flags:    
    if 'type' in flag and flag['type'] != bool:
        parser.add_argument(flag['short'], flag['long'], dest=flag['name'], metavar=flag['metavar'], type=flag['type'], help=f"{flag['help']} (default = {flag['default']})", default=flag['default'])
    else:
        parser.add_argument(flag['short'], flag['long'], dest=flag['name'], metavar=flag['metavar'], help=f"{flag['help']} (default = {flag['default']})", default=flag['default'])
args = parser.parse_args()
print("Chosen options:")
for flag in flags:
	if getattr(args, flag['name']) == None:
	    print(f"{flag['name'].capitalize()} ({flag['short']}): Not specified")
	else:
	    print(f"{flag['name'].capitalize()} ({flag['short']}): {getattr(args, flag['name'])}")
if args.image_path:
    if os.path.exists(args.image_path) and os.path.isfile(args.image_path):
        print(f"\nTreating image: {args.image_path}")
    else:
        print(f"Invalid image path: {args.image_path}\n")
        sys.exit(0)
else:
    print("\nYou must specify image path.\n")
    sys.exit(0)

sa = stringart(	args.scalefactor, 
				args.weight,
				args.lines,
				args.linewidth,
				args.invertimage, 
				args.contour,
				args.autosaturation,
				args.resolution,
				args.image_path)

image, dim, inputimage = sa.prepareimage()
inputdimension = N.min([N.shape(inputimage)[0],N.shape(inputimage)[1]])
print(f'Input image dimensions 		{N.shape(inputimage)[0]}x{N.shape(inputimage)[1]}')
print(f'Input image trimmed dim. 	{inputdimension}x{inputdimension}')

circ_coords = N.array(sa.circumferencecoordinates(image,dim)).T
t1 = time.process_time()

linecoords, diff, image = sa.calculatestrings()

t2 = time.process_time()
linecoords = N.array(linecoords)

print(f'\nGeneration time: {round(t2-t1,1)} s.')
print(f'Strings used: {len(linecoords)}')

original_image_path = args.image_path
filename, extension = os.path.splitext(os.path.basename(original_image_path))

singlelinex = [linecoords[ii].T[0] for ii in range(len(linecoords))]
singleliney = [linecoords[ii].T[1] for ii in range(len(linecoords))]

resolutions = {'4k':1,'HD':2,'VGA':4}
chosenres = [resolutions.get(sa.resolution,1) if sa.resolution in resolutions else 2][0]
if sa.resolution not in resolutions:
	print('Invalid resolution chosen: Default HD res applied')
yr = 28.2/chosenres
xr = yr
mpl.rcParams['figure.figsize'] = xr,yr

fig, ax = plt.subplots()
sa.invertimage = sa.invertimage.lower() == 'true'
fig.patch.set_facecolor('white' if sa.invertimage else 'black')
line_color = 'black' if sa.invertimage else 'white'
plt.plot(singlelinex, singleliney, '-', color=line_color, linewidth=sa.linewidth)

plt.gca().invert_yaxis()
plt.gca().set_aspect('equal')
ax.axis("off")
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig(f'stringart_{filename}.jpg', bbox_inches='tight', pad_inches=0)
print(f'\nProgram finished. Image stringart_{filename}.png saved.\n')

