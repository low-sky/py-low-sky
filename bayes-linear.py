#!/usr/bin/env python
"""
bayes-linear.py -- Attempt to implement the Gibbs sampler in pytho
"""

import numpy as np
import scipy.stats as s
import matplotlib.pyplot as p
import emcee 

# Generate some data to fit
x = s.uniform.rvs(size=((30)),scale=10)
y = 3*x+2.2+s.norm.rvs(size=((30)),scale=0.1)


def lnprob(p,x=x,y=y):
#    x = np.array([ 8.09379819,  8.96357024,  7.66253652,  1.44041006,  9.72833801,
#                   9.95047447,  9.70450228,  2.32545068,  6.32867266,  9.31727735])
#    y= np.array([ 35.75998072,  35.48902694,  32.28364158,   7.62223271,
#                  38.38233809,  44.9036255 ,  42.09450263,  11.16334563,
#                  27.45946439,  38.08215968])

    lp = -np.sum((y-(p[1]*x+p[0]))**2)*p[2]*0.5+s.invgamma.logpdf(p[2],[1,20])
    return lp

ndim, nwalkers = 3,100
#ivar = s.invgamma.rvs([1,2],size=((ndim)),scale=2.0,loc=1)
p0 = [np.random.rand(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)#, args=[x,y])
sampler.run_mcmc(p0, 1000)

#p.plot(x,y,'r+')
#p.show()
# Begin a bayes loop



