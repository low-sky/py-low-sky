#!/usr/bin/env python
import sys
import astropy.io.fits as fits
from astropy.coordinates import FK5,SkyCoord
import astropy.units as u
from math import cos


def bary2lsr(hdr):
    Obs_Direction = SkyCoord(hdr['CRVAL1'],hdr['CRVAL2'],frame='fk5',
                             unit="deg")
    LSR_Direction = SkyCoord("18h03m50.29s +30d00m16.8s",frame='fk5')
    sep = Obs_Direction.separation(LSR_Direction)
    mag = 20*u.km/u.s*cos(sep.value)
    return(mag)

if __name__ == "__main__":
    if len(sys.argv)>1:
        inputfile = sys.argv[1]
        hdr = fits.getheader(inputfile)
        magnitude = bary2lsr(hdr)
        print("To convert BARYCENT to LSRK, add {0} to the velocity.".format(magnitude))
