import numpy as np

def median(theta,nan=False):
    """
    Calculates the median of a set of angles on the circle and returns
    a value on the interval from (-pi, pi].  Angles expected in radians.

    Setting the nan kwarg to True uses NaN-robust numpy functions.
    """
    if nan:
        meanfunc = np.nanmedian
    else:
        meanfunc = np.median
    medangle = np.arctan2(meanfunc(np.sin(theta)),meanfunc(np.cos(theta)))
    return(medangle)

def mean(theta,nan=False):
    """
    Calculates the mean of a set of angles on the circle and returns
    a value on the interval from (-pi, pi].  Angles expected in radians.

    Setting the nan kwarg to True uses NaN-robust numpy functions.
    """

    if nan:
        meanfunc = np.nanmean
    else:
        meanfunc = np.mean
    meanangle = np.arctan2(meanfunc(np.sin(theta)),meanfunc(np.cos(theta)))
    return(meanangle)

def std(theta,nan=False):
    """
    Calculates the standard deviation of a set of angles on the circle
    and returns a value on the interval from (-pi, pi].  Angles
    expected in radians.

    Setting the nan kwarg to True uses NaN-robust numpy functions.
    """

    if nan:
        meanfunc = np.nanmean
    else:
        meanfunc = np.mean
    meanangle = np.arctan2(meanfunc(np.sin(theta)),meanfunc(np.cos(theta)))
    diffcos = np.cos(theta)*np.cos(meanangle)+np.sin(theta)*np.sin(meanangle)
    diffsin = np.sin(theta)*np.cos(meanangle)-np.cos(theta)*np.sin(meanangle)
    diffangle = np.arctan2(diffsin,diffcos)
    stdangle = np.sqrt(meanfunc(diffangle**2))
    return(stdangle)

def mad(theta,nan=False):
    """
    Calculates the median absolute deviation of a set of angles on the
    circle and returns a value on the interval from (-pi, pi].  Angles
    expected in radians.

    Setting the nan kwarg to True uses NaN-robust numpy functions.
    """
    if nan:
        meanfunc = np.nanmedian
    else:
        meanfunc = np.median
    medangle = np.arctan2(meanfunc(np.sin(theta)),meanfunc(np.cos(theta)))
    diffcos = np.cos(theta)*np.cos(medangle)+np.sin(theta)*np.sin(medangle)
    diffsin = np.sin(theta)*np.cos(medangle)-np.cos(theta)*np.sin(medangle)
    madangle = meanfunc(np.abs(np.arctan2(diffsin,diffcos)/0.6745))
    return(madangle)

def percentile(theta,pct,nan=False):
    """
    Calculates the percentile rank of an angle in a distribution from
    set of angles on the circle and returns a value on the interval
    from (-pi, pi].  Angles expected in radians.

    Setting the nan kwarg to True uses NaN-robust numpy functions.
    """
    if nan:
        meanfunc = np.nanmedian
    else:
        meanfunc = np.median
    medangle = np.arctan2(meanfunc(np.sin(theta)),meanfunc(np.cos(theta)))
    diffcos = np.cos(theta)*np.cos(medangle)+np.sin(theta)*np.sin(medangle)
    diffsin = np.sin(theta)*np.cos(medangle)-np.cos(theta)*np.sin(medangle)
    diffangle = np.arctan2(diffsin,diffcos)
    pctangle = np.percentile(diffangle,pct)+medangle
    pctangle = np.arctan2(np.sin(pctangle),np.cos(pctangle))
    return(pctangle)

def difference(theta,phi):
    """
    Calculates the difference between two broadcastable arrays of
    angles and returns a difference in radians.
    """
    diffcos = np.cos(theta)*np.cos(phi)+np.sin(theta)*np.sin(phi)
    diffsin = np.sin(theta)*np.cos(phi)-np.cos(theta)*np.sin(phi)
    diffangle = np.arctan2(diffsin,diffcos)
    return(diffangle)
