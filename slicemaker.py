#!/usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt

#if len(sys.argv)==1:

filename = sys.argv[1]
fileroot = filename.replace('.fits','')
fileroot = fileroot.replace('.gz','')
print(fileroot)
cube = fits.getdata(filename)
vmn = 0.0
vmx = np.percentile(cube[np.isfinite(cube)],99.9)
planes = np.arange(10)+300
#for plane in planes:
for plane in np.arange(cube.shape[0]):
    plt.imshow(cube[plane,:,:],vmin=vmn,vmax=vmx,
               cmap='Greys',origin='lower',
               interpolation='nearest')
    plt.axis('off')
    plt.savefig(fileroot+'.{0}.png'.format(str(plane).zfill(3)),bbox_inches='tight')
    plt.clf()

