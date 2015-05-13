import astropy
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
#from astroy.wcs import WCS

class Galaxy:
    def __init__(self,Name=None):
        if not Name:
            self.Name = 'M33'
            self.Distance = 8.4e5*u.pc
            self.CenterPosition = SkyCoord(23.461667,30.660194,unit=(u.deg,u.deg),frame='fk5')
            self.PositionAngle = 202*u.deg
            self.Inclination = 56*u.deg
            self.Vsys = -179*u.km/u.s
    
def DiskCoord(ra = None, dec = None, fitsheader = None, galaxy = Galaxy(),
              Xgal = None, Ygal = None, returnXY = False):
#    if fitsheader:
#        w = WCS(fitsheader)
#        whichaxes = ['celestial' in axis['coordinate_type'] for axis in w.get_axis_types()]
#        xaxis = next(i+1 for i in range(w.naxis) if whichaxes[i])        
    Offsets = SkyCoord(ra,dec,unit=(u.deg,u.deg))
    PAs = galaxy.CenterPosition.position_angle(Offsets)
    GalPA = PAs - galaxy.PositionAngle
    GCDist = Offsets.separation(galaxy.CenterPosition)
    # Transform into galaxy plane
    Rplane = galaxy.Distance*np.tan(GCDist)
    Xplane = Rplane * np.cos(GalPA)
    Yplane = Rplane * np.sin(GalPA)
    Xgal = Xplane
    Ygal = Yplane / np.cos(galaxy.Inclination)
    Rgal = (Xgal**2+Ygal**2)**0.5
    if returnXY:
        return (Xgal,Ygal)
    else:
        return Rgal
    
