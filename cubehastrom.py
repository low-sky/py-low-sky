from astropy.io import fits
from skimage import transform as tf
import astropy.wcs as wcs
import numpy as np

def cubehastrom(WarpHDU, TemplateHDU, anchors = 5, order = 2):

    normal = np.linspace(0,1,anchors)
    naxis = WarpHDU.header.get('NAXIS', 0)
    if naxis == 1:
        xx = normal
    elif naxis == 2:
        xx,yy = np.meshgrid(normal,normal,indexing='xy')
        bundle = (xx,yy)
    elif naxis == 3: 
        xx,yy,zz = np.meshgrid(normal,normal,normal,indexing='xy')
        bundle = (xx,yy,zz)
    if naxis > 0:
        size = ()
        for idx in range(naxis):
            size +=(WarpHDU.header['NAXIS' + str(idx + 1)],)

    naxisTemplate = TemplateHDU.header.get('NAXIS',0)
    if naxisTemplate > 0:
        sizeTemplate = ()
        for idx in range(naxis):
            sizeTemplate = sizeTemplate + (TemplateHDU.header['NAXIS' + str(idx + 1)],)

    # Build Output Image
    OutputSize = ((np.c_[size,sizeTemplate]).max(axis=1))
    # Switch first two elements to make array sizes pythonic
    OutputSize[0],OutputSize[1],OutputSize[2] = OutputSize[2],OutputSize[1],OutputSize[0]
    
    OutputImg = np.zeros(OutputSize)

    # Calculate XY positions for anchor points in image
    XY = ()
    for elt in zip(size,bundle):
        XY += (elt[0]*elt[1],)
    XY += (0,) # This appends the final argument to WCS functions

    # Put WCS coordinates into place in the new image
    WarpWCS = wcs.WCS(WarpHDU.header)
    TemplateWCS = wcs.WCS(TemplateHDU.header)
    WarpWorld = WarpWCS.wcs_pix2world(*XY)
    WarpWorld += (0,)
    import pdb; pdb.set_trace()
    TemplateXY = TemplateWCS.wcs_world2pix(*WarpWorld)
    # Assemble mapping coordinates
    src = np.c_[XY[0].flatten(),XY[1].flatten(),XY[2].flatten()]
    dst = np.c_[TemplateXY[0].flatten(),TemplateXY[1].flatten(),TemplateXY[2].flatten()]

    # Exectue Transform
    xform = tf.estimate_transform('polynomial', dst, src, order = order)
    import pdb; pdb.set_trace()
    OutputImg[0:size[2],0:size[1],0:size[0]] = WarpHDU.filled_data[:]
    OutputImg = tf.warp(OutputImg, inverse_map=xform) 

    # Crop image to template size
    OutputImg = OutputImg[0:sizeTemplate[1],0:sizeTemplate[0]]

    TemplateHdr = TemplateWCS.to_header()
    OutputHdr = WarpHDU.header.copy()
    for keyword in TemplateHdr.keys():
        OutputHdr[keyword] = TemplateHdr[keyword]
    OutputHDU = fits.PrimaryHDU(data = OutputImg,header = OutputHdr)
    return(OutputHDU)
