import astropy.wcs as wcs
import astropy.io.fits as fits
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian2DKernel
from skimage.morphology import dilation,disk
import glob
from astropy.table import Table
import pdb 
cmcdir = '/srv/astro/erickoch/gould_belt/exclude_from_filpaperanalysis/'
cmcfile = 'aur_cntr_scanam_spire350_0.fits'

planckfile = 'HFI_SkyMap_857_2048_R1.10_nominal.fits'
planckdata = hp.read_map(planckfile)

cmcdir = '/srv/astro/erickoch/gould_belt/'
flist = glob.glob(cmcdir+'*350.fits')

#    = 'california_cntr-350.fits'
#cmcfile = 'chamaeleonI-350.fits'

t = Table(names=('region','filename','offset','slope'),
          dtype = ['a60','a60','f8','f8'])


#flist = ['/srv/astro/erickoch/gould_belt/aquilaM2-350.fits']

for cmcfile in flist:
    t.add_row()
    t['filename'][-1] = cmcfile
    root = (cmcfile.split('/'))[-1]
    root = (root.split('-'))[0]    
    print('Now processing {0}'.format(root))
    orig_img = fits.getdata(cmcfile)
    hdr = fits.getheader(cmcfile)
    if hdr['BUNIT']=='Jy/beam':
        jypb2Mjysr = 1e-6/(2*np.pi*(25.16*u.arcsec)**2/(8*np.log(2))).to(u.sr)
    else:
        jypb2Mjysr = 1.0
    orig_img = orig_img*jypb2Mjysr

    width = 5/60./hdr['CDELT2']/np.sqrt(8*np.log(2))

    kernel = Gaussian2DKernel(width)
    img = convolve_fft(orig_img,kernel,min_wt=4.0,
                       interpolate_nan=False,normalize_kernel=True)
    badmask = np.isnan(orig_img)
    struct = disk(width*4)
    badmask = dilation(badmask,struct)
    img[np.where(badmask)]=np.nan

    w = wcs.WCS(hdr)
    x,y = np.meshgrid(np.arange(img.shape[1]),
                      np.arange(img.shape[0]),
                      indexing='xy')

    c = wcs.utils.pixel_to_skycoord(x,y,w,0)
    l = c.galactic.l.radian
    b = c.galactic.b.radian

    #idx = hp.ang2pix(2048,l,b)
    planckmap = (hp.pixelfunc.get_interp_val(planckdata,np.pi/2-b,l,nest=False))
    #planckidx = hp.pixelfunc.ang2pix(2048,b,l)

    good = np.where(np.isfinite(planckmap)*np.isfinite(img))
    xmin = np.percentile(planckmap[good],1)
    xmax = np.percentile(planckmap[good],99)
    ymin = np.percentile(img[good],1)
    ymax = np.percentile(img[good],99)
    m, b = np.polyfit(planckmap[good],img[good],1)
    plt.hexbin(planckmap.ravel(),img.ravel()-b,cmap='gray_r',
                bins='log',gridsize=[100,100])
    print('Offset is {0} MJy/sr'.format(b))
    t['offset'][-1] = b
    t['slope'][-1] = m
    t['region'][-1] = root
    plt.plot(np.linspace(0,1000,20),np.linspace(0,1000,20))
    plt.xlabel(r'Planck 350 $\mu$m')
    plt.ylabel(r'Herschel 350 $\mu$m')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin-b,ymax-b])
    
    plt.savefig('planck_vs_{0}.png'.format(root))
    plt.close()
    plt.clf()

    fig = plt.figure(figsize=(10,5))
    plt.subplot(121)
    plt.imshow(img,vmax=ymax,vmin=ymin)
    
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(planckmap,vmax=xmax,vmin=xmin)
    plt.colorbar()
    plt.savefig('maps_{0}.png'.format(root))
    plt.close()
    plt.clf()
    
#    pdb.set_trace()

t.write('offset_values.csv',format='ascii.csv')
