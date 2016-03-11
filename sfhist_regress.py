#!/usr/bin/env python -tt

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plot

newparams = np.zeros(10000)

highztable = np.genfromtxt('error_fit_highz')

sfruv = highztable[:,2]
sfrir = highztable[:,3]
sfrha = highztable[:,0]
sfrhaerr = highztable[:,1]

uvesc = highztable[:,2]/(highztable[:,2]+highztable[:,3])
haesc = highztable[:,0]/(highztable[:,2]+highztable[:,3])
haescerr = highztable[:,1]/(highztable[:,2]+highztable[:,3])

converttable = np.genfromtxt('KFUV_KIR.dat')

time = 10.**(converttable[:,0])

kuvtime = 10.**(converttable[:,1])

kuvfinal = 1.4e-28*3.826e33/(3.e18/1600.)

uvcorrtime = kuvtime/kuvfinal

uvzerofind = interp1d(time, uvcorrtime-1.)

uvduration = brentq(uvzerofind, time[0], time[-1])

kirtime = 10.**converttable[:, 2:]

def escregress(x, q):
  return np.power(x, q)

param0 = .55
  
#param, paramcov = curve_fit(escregress, uvesc, haesc, sigma = haescerr, \
#absolute_sigma = True, p0 = param0)

mark_uvesc = 10.**(-.4*np.linspace(start = .55, stop = 5.5, num = 10))

closest_uvesc = np.abs(np.subtract.outer(mark_uvesc, uvesc)).argmin(0)

for j in xrange(10000):

  uvnewdurations = np.random.uniform(low = uvduration, \
  high = uvduration + 9.e7, size = len(uvesc))
  
  uvcorr_newdurations = uvzerofind(uvnewdurations)+1.
  
  newsfruv = sfruv*uvcorr_newdurations

  newsfrir = np.zeros(sfrir.shape)
  
  for i in xrange(len(mark_uvesc)):
    irinterp = interp1d(time, kirtime[:, i])
  #  irduration = brentq(irzerofind, time[0], time[-1])
  #  num_gals = (closest_uvesc == i).sum()
  #  irnewdurations = np.random.uniform(low = irduration - 7.e7, \
  #  high = uvduration + 7.e7, size = num_gals)
    ind_gals = (closest_uvesc == i)
    irnewdurations = uvnewdurations[ind_gals]
    kirtime_gals = irinterp(irnewdurations)
    kirfinal_gals = irinterp(uvduration)
  #  kirfinal_gals = 4.5e-44*3.826e33
    newsfrir[ind_gals] = sfrir[ind_gals]*kirtime_gals/kirfinal_gals
  
  newuvesc = newsfruv/(newsfruv+newsfrir)
  newhaesc = sfrha/(newsfruv+newsfrir)
  newhaescerr = sfrhaerr/(newsfruv+newsfrir)
  
  newparam, newparamcov = curve_fit(escregress, newuvesc, newhaesc, \
  sigma = newhaescerr, absolute_sigma = True, p0 = param0)
  
  newparams[j] = newparam
  
plot.hist(newparams, bins = 20)
plot.show()