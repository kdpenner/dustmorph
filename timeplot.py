#!/usr/bin/env python -tt

import numpy
import matplotlib.pyplot as plot
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.ticker import MultipleLocator
from scipy.optimize import brentq
from scipy.interpolate import interp1d


table = numpy.genfromtxt('KFUV_KIR.dat')

uvfescs = 10.**(-.4*numpy.linspace(start = .55, stop = 5.5, num = 10))

sfruv = uvfescs

#sfruv = numpy.linspace(start = .0001, stop = .0002, num = 10)

sfrir = 1-sfruv

#sfrir = sfruv/uvfescs - sfruv

kuvtime = 10.**(table[:,1])

kuvfinal = 1.313e-28*3.826e33/(3.e18/1500.)

#kuvfinal = kuvtime[-1]

#print kuvfinal

uvcorrtime = kuvtime/kuvfinal

newsfruv = numpy.dot(numpy.matrix(uvcorrtime).transpose(), numpy.matrix(sfruv))

time = 10.**table[:,0]

uvzerofind = interp1d(time, uvcorrtime-1.)

uvduration = brentq(uvzerofind, time[0], time[-1])

kirtime = 10.**table[:, 2:]

#kirfinal = 4.5e-44*3.826e33

#kirfinal = kirtime[-1]

#print kirfinal

kirfinal = numpy.zeros(uvfescs.shape)

for i in xrange(len(uvfescs)):
  irinterp = interp1d(time, kirtime[:, i])
#  kirfinal[i] = irinterp(uvduration)
  kirfinal[i] = 4.5e-44*3.826e33

ircorrtime = kirtime/kirfinal

sfrirm = numpy.identity(sfrir.size)
repind = numpy.diag_indices(sfrir.size)
sfrirm[repind] = sfrir

newsfrir = numpy.dot(numpy.matrix(ircorrtime), numpy.matrix(sfrirm))

newsfrtot = newsfruv+newsfrir

goodtimes = numpy.where((time > 10.**7.3) & (time < 10.**9.3))[0]

newsfrgood = newsfrtot[goodtimes, :]
newsfruvgood = newsfruv[goodtimes, :]

rangeuvfescs = (newsfruvgood/newsfrgood).transpose()

hafescs = uvfescs**0.537

sfrha = (sfruv+sfrir)*hafescs

rangehafescs = (sfrha/newsfrgood).transpose()

plot.rcParams['figure.figsize'] = [7.32, 7.32]
color1 = '#d95f02'
color2 = '#1b9e77'

#fig, ax1 = plot.subplots()

#ax1.plot(uvfescs, rangehafescs, marker = 'o', color = 'k')

#ax1.plot(uvfescs, rangehafescs[:, -1], marker = 'o', color = color2)

#ax1.plot([.001, 1], numpy.power([.001, 1], 0.55), linestyle = '-', \
#color = color1, linewidth = 4, marker = None, \
#label = r'f$_{esc}$(H$\alpha$) = f$_{esc}$(0.16$\mu$m)$^{0.55}$')

#ax1.set_ylim([.001, 1])
#ax1.set_xlim([.001, 1])
#ax1.set_aspect('equal')
#ax1.set_xlabel(r'f$_{esc}$(0.16$\mu$m) (L$_{emergent}$/L$_{intrinsic}$)')
#ax1.set_ylabel(r'f$_{esc}$(H$\alpha$) (L$_{emergent}$/L$_{intrinsic}$)')
#ax1.set_xscale('log')
#ax1.set_yscale('log')

#ax1.legend(loc = 'lower right', prop = {'size': 18})
  
#ax1.set_yticks([.001, .01, .1, 1])
#ax1.set_yticklabels([r'$10^{-3}$', r'10$^{-2}$', r'10$^{-1}$', r'  1  '])
#ax1.set_xticks([.001, .01, .1, 1])
#ax1.set_xticklabels([r'$10^{-3}$', r'10$^{-2}$', r'10$^{-1}$', r'1'])

#ax2 = ax1.twinx()
#ax2.set_ylim([-2.5*numpy.log10(.001),0])
#ax2.yaxis.set_major_locator(MultipleLocator(2))
#ax2.yaxis.set_minor_locator(MultipleLocator(.5))
#ax2.set_ylabel(r'A$_{H\alpha}$ (mag)')

#plot.savefig('spread_ha_age1.eps')

#plot.close()

fig, ax1 = plot.subplots()

ax1.plot([.001, 1], numpy.power([.001, 1], 0.537), linestyle = '-', \
color = color1, linewidth = 4, marker = None, \
label = r'f$_{esc}$(H$\alpha$) = f$_{esc}$(0.16$\mu$m)$^{0.537}$')

ax1.plot(rangeuvfescs, rangehafescs, marker = 'o', color = 'k', \
linestyle = 'None')

input = mlines.Line2D(uvfescs, hafescs, marker = 'o', color = color2, \
linestyle = 'None', label = r'given f$_{esc}$', markeredgecolor = color2)

ax1.add_line(input)

ax1.set_ylim([.001, 1])
ax1.set_xlim([.001, 1])
ax1.set_aspect('equal')
ax1.set_xlabel(r'f$_{esc}$(0.16$\mu$m) (L$_{emergent}$/L$_{intrinsic}$)')
ax1.set_ylabel(r'f$_{esc}$(H$\alpha$) (L$_{emergent}$/L$_{intrinsic}$)')
ax1.set_xscale('log')
ax1.set_yscale('log')

handles, labels = ax1.get_legend_handles_labels()

#print handles, labels

extra = mlines.Line2D([], [], color = 'k', marker = 'o', linestyle = 'None')

handles.append(extra)
labels.append('true f$_{esc}$ for a different duration\nof a star formation episode')

handdict = {extra: HandlerLine2D(numpoints = 1), \
input: HandlerLine2D(numpoints = 1)}

ax1.legend(handles, labels, loc = 'lower right', prop = {'size': 18}, \
frameon = False, handler_map=handdict)

ax1.set_yticks((.001, .01, .1, 1))
ax1.set_yticklabels((r'$10^{-3}$', r'10$^{-2}$', r'10$^{-1}$', r'  1  '))
ax1.set_xticks((.001, .01, .1, 1))
ax1.set_xticklabels((r'$10^{-3}$', r'10$^{-2}$', r'10$^{-1}$', r'1'))

ax2 = ax1.twinx()
ax2.set_ylim([-2.5*numpy.log10(.001),0])
ax2.yaxis.set_major_locator(MultipleLocator(2))
ax2.yaxis.set_minor_locator(MultipleLocator(.5))
ax2.set_ylabel(r'A$_{H\alpha}$ (mag)')

plot.savefig('spread_hauv_age.eps')
