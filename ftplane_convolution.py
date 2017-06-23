from spectral_cube import SpectralCube, VaryingResolutionSpectralCube
from radio_beam import Beam
import numpy as np
import astropy.units as u
import astropy.utils.console as console
import copy
import warnings

def ftconvolve(ImageIn, major = 1.0, minor = 1.0,
               angle = 0.0):
    NanMaskFlag = False
    nanmask = np.isnan(ImageIn)
    image = np.copy(ImageIn)
    if np.any(nanmask):
        wtimg = np.ones_like(ImageIn)
        NanMaskFlag = True
        image[nanmask] = 0.0
        wtimg[nanmask] = 0.0
        ftwtimg = np.fft.fftn(wtimg)    

    ftimg = np.fft.fftn(image)

    if major == 0.0:
        sigmau = np.inf
    else:
        sigmau = 0.5/(np.pi * major)

    if minor == 0.0:
        sigmav = np.inf
    else:
        sigmav = 0.5/(np.pi * minor)

    FTPA = angle + np.pi

    a = 0.5 * (np.cos(FTPA)**2 / sigmau**2 +
               np.sin(FTPA)**2 / sigmav**2)
    c = 0.5 * (np.sin(FTPA)**2 / sigmau**2 +
               np.cos(FTPA)**2 / sigmav**2)
    b = 0.25 * np.sin(2 * FTPA) * (1.0 / sigmav**2 - 1.0 / sigmau**2)

    vv, uu = np.meshgrid(np.fft.fftfreq(ftimg.shape[0]),
                         np.fft.fftfreq(ftimg.shape[1]),
                         indexing='ij')

    FTkernel = np.exp(-a*(uu)**2 -c*(vv)**2 +2*b*(uu*vv))
    ConvolvedImage = (np.fft.ifftn(ftimg * FTkernel)).real

    if NanMaskFlag:
        ConvolvedMask = (np.fft.ifftn(ftwtimg * FTkernel)).real
        ConvolvedImage /= ConvolvedMask
        ConvolvedImage[nanmask] = np.nan
    return(ConvolvedImage)

def MakeRoundBeam(incube, outfile=None):

    '''

    This takes a FITS file or a SpectralCube and outputs a cube with a
    single round beam of maximum size.  No provision for robust
    statistics is currently made so a single plane with a very large
    beam will end up being the target.

    Parameters
    ----------
    filename : `string` or `SpectralCube`
       Input spectral cube
    
    Returns
    -------
    cube : `SpectralCube`

    '''
    if isinstance(incube,str):
        cube = SpectralCube.read(incube)

    if isinstance(incube, VaryingResolutionSpectralCube):
        cube = incube

    if not isinstance(cube, VaryingResolutionSpectralCube):
        warnings.warn("No information about multiple beams")
        return(None)
    
    beams = cube.beams
    major = np.array([bm.major.to(u.deg).value for bm in beams])
    minor = np.array([bm.minor.to(u.deg).value for bm in beams])
    PAs = np.array([bm.pa.to(u.rad).value for bm in beams])
    target = np.array(major.max())
    # Let's assume square pixels
    pixsize = cube.wcs.pixel_scale_matrix[1,1]
    output = np.zeros(cube.shape)
    fwhm2sigma = np.sqrt(8*np.log(2))

    with console.ProgressBar(cube.shape[0]) as bar:
        for ii, plane in enumerate(cube.filled_data[:]):
            majpix = (target**2 - major[ii]**2)**0.5 /\
                     pixsize / fwhm2sigma
            minpix = (target**2 - minor[ii]**2)**0.5 /\
                     pixsize / fwhm2sigma
            output[ii,:,:] = ftconvolve(plane,
                                        major = majpix,
                                        minor = minpix,
                                        angle = PAs[ii])
            bar.update()

    hdr = copy.copy(cube.header)
    hdr['CASAMBM'] = False
    hdr['BMAJ'] = float(target)
    hdr['BMIN'] = float(target)
    hdr['BPA'] = 0.0
    outcube = SpectralCube(output, cube.wcs, header=hdr)
    if outfile:
        outcube.write(outfile)
        return None
    return(outcube)


