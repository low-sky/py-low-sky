#!/usr/bin/env python
"""

samples = BayesRatio(X, Y, Xerror, Yerror, nWalkers=100,nBurn=100,nSample=1000)

X,Y -- Data vectors (same length, can be length 1)
Xerror, Yerror -- 1 sigma errors for X,Y
nWakers = 10 (default), number of walkers in the sampler (>2 required)
nBurn = 100 (default), number of steps to burn in chain for
nSample = 1000 (default), number of steps to sample chain with

Returns:
samples = a numpy array of samples of the ratio of the data
"""
import numpy as np
import emcee
import matplotlib.pyplot as p

def logprob(p,x,y,x_err,y_err):
    theta = p[0]
    if np.abs(theta-np.pi/4)>np.pi/4:
        return -np.inf
    Delta = (np.cos(theta)*y - np.sin(theta)*x)**2
    Sigma = (np.sin(theta))**2*x_err**2+(np.cos(theta))**2*y_err**2
    lp = -0.5*np.nansum(Delta/Sigma)-0.5*np.nansum(np.log(Sigma))

    return lp

def BayesRatio(x,y,x_err,y_err,nWalkers=10,nBurn=100,nSample=1000):
    ndim = 1
    p0 = np.zeros((nWalkers,ndim))
    p0[:,0] = np.pi/4+np.random.randn(nWalkers)*0.1

    sampler = emcee.EnsembleSampler(nWalkers, ndim, logprob, 
                                    args=[x,y,x_err,y_err])
    pos, prob, state = sampler.run_mcmc(p0, nBurn)
    sampler.reset()
    sampler.run_mcmc(pos, nSample)
    return(np.tan(sampler.flatchain))
    
