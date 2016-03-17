#!/usr/bin/env python -tt

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import matplotlib.pyplot as plot
from matplotlib.ticker import FormatStrFormatter

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

kuvfinal = 1.313e-28*3.826e33/(3.e18/1500.)

uvcorrtime = kuvtime/kuvfinal

uvzerofind = interp1d(time, uvcorrtime-1.)
uvtimeinterp = interp1d(time, kuvtime)

uvduration = brentq(uvzerofind, time[0], time[-1])

kirtime = 10.**converttable[:, 2:]

def escregress(x, q):
  return np.power(x, q)

param0 = .55
  
#param, paramcov = curve_fit(escregress, uvesc, haesc, sigma = haescerr, \
#absolute_sigma = True, p0 = param0)

#print param

mark_uvesc = 10.**(-.4*np.linspace(start = .55, stop = 5.5, num = 10))

closest_uvesc = np.abs(np.subtract.outer(mark_uvesc, uvesc)).argmin(0)

for j in xrange(10000):

  uvnewdurations = np.random.uniform(low = 10.**7.3, \
  high = 10.**9.3, size = len(uvesc))
  
  #uvnewdurations = np.zeros(len(uvesc))+uvduration
  
  uvcorr_newdurations = uvtimeinterp(uvnewdurations)/kuvfinal
  
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
  #  kirfinal_gals = irinterp(uvduration)
    kirfinal_gals = 4.5e-44*3.826e33
    newsfrir[ind_gals] = sfrir[ind_gals]*kirtime_gals/kirfinal_gals
  
  newuvesc = newsfruv/(newsfruv+newsfrir)
  newhaesc = sfrha/(newsfruv+newsfrir)
  newhaescerr = sfrhaerr/(newsfruv+newsfrir)
  
  newparam, newparamcov = curve_fit(escregress, newuvesc, newhaesc, \
  sigma = newhaescerr, absolute_sigma = True, p0 = param0)
  
  newparams[j] = newparam

plot.rcParams['figure.figsize'] = [7.32, 7.32]
majorFormatter = FormatStrFormatter('%0.3f')

fig, ax = plot.subplots()
plot.hist(newparams, bins = 16, histtype = 'step', color = 'black', \
range = (0.52, 0.56))
ax.set_xlabel(r'best-fit q for f$_{esc}$(H$\alpha$) = f$_{esc}$(0.16$\mu$m)$^{q}$')
ax.set_ylabel(r'Number')
ax.xaxis.set_major_formatter(majorFormatter)
ax.set_xticks(np.arange(.52, .56, .01))
ax.set_xticks(np.arange(.52, .56, .0025), minor = True)

plot.savefig('dist_hauv_age.eps')