import numpy as np



def PAmoment(weight, x,    y, ReturnSizes = False):
    """
    Calculates the position angle (measured counter-clockwise from the
    x-axis) for a set of x,y data with weights.  This is the
    moment-of-inertia approach.

    Parameters
    ----------
    weight : float
        Weight for each point.  Should be non-negative.
    x, y : float
        coordinates of each point
    ReturnSizes : bool
        If True, return major and minor axis sizes along with position angle.

    Returns
    -------
    pa : float
        Position angle measured in radians from the x-axis
    major : float
        Major axis size of the points
    minor : float
        Minor axis size of the points
    """


    sumwts = np.nansum(weight)
    xcen = np.nansum(x*weight)/np.nansum(weight)
    ycen = np.nansum(y*weight)/np.nansum(weight)
    matrix = 1/sumwts*np.array([[np.nansum(weight*(x-xcen)**2),
                                 np.nansum(weight*(x-xcen)*(y-ycen))],
                                [np.nansum(weight*(x-xcen)*(y-ycen)),
                                 np.nansum(weight*(y-ycen)**2)]])
    determ = np.linalg.det(matrix)
    if ~np.isfinite(determ) or determ == 0:
        return(np.nan)
    evals, evecs = np.linalg.eigh(matrix)
    bigvec = evecs[-1,:]
    pa = np.arctan2(bigvec[1],bigvec[0])
    if ReturnSizes:
        major = evals[-1]
        minor = evals[0]
        return (pa,major,minor)
    return(pa)
