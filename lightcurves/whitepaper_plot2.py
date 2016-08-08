#!/usr/bin/env python
from __future__ import print_function
+__author__ = 'stefano valenti'
import numpy as np
import pylab as pl
import glob
import os
from numpy import interp as ninterp
import socket

#pl.ion()
#FBB if running on an interactive shell activate pl.ion() within the shel before execing the code
pl.rcParams['figure.figsize'] = 17, 6
#fig = pl.figure(frameon=True)
ax = pl.axes([.06, .1, .27, .8])
ax2 = pl.axes([.39, .1, .27, .8])
ax3 = pl.axes([.72, .1, .27, .8])

#ax = fig.add_subplot(2,1,1,frame_on=True)
#ax2 = fig.add_subplot(2,1,2,frame_on=True)

def colorrange(N):
    N = N / 3.
    color1 = {}
    for i in range(0, int(N) + 1, 1):
        color1[int(i)] = (0.0, i / N, 1)
        color1[int(i + N)] = (i / N, 1, 1 - i / N)
        color1[int(i + N * 2)] = (1, 1 - i / N, 0.)

    color2 = {}
    for j in color1.keys():
        print ("color counter", j)
        color2[j] = color1[int(len(color1) - j - 1)]
    return color2

########################################

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def readmodels(model='shen'):
    models = {}
    if model == 'shen':
        _dir = 'models_shen2010/'
        files = glob.glob(_dir + '*2.dat')
        print ("input files: ", files)
        for f in files:
            data = np.genfromtxt(f, 'float')
            t, M_bol, M_U, M_B, M_V, M_R, M_I, M_J, M_H, M_K = zip(*data)
            models[f.split('/')[-1]] = {}
            models[f.split('/')[-1]]['days'] = np.array(t)
            models[f.split('/')[-1]]['mag'] = (np.array(M_V) * (-1)) + np.max(M_V)
    elif model == 'temp':
        _dir = 'SN_template/'
        for asci in ['snIa_K211b.dat', 'snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat', 'PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat', \
                     'snIc_2002ap_V.dat', 'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat', \
                     'SNIa_CSM_r.dat']:
            data = np.genfromtxt(_dir + asci, 'float')
            print ("this input is:", asci)
            if asci in ['snIa_K211b.dat', 'snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat', 'PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat', \
                          'snIc_2002ap_V.dat', 'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat', \
            'SNIa_CSM_r.dat']:
                ph, V = zip(*data)
            else:
                ph, u, B, V, r, i = zip(*data)
            ph = np.array(ph)
            V = np.array(V)
            # I'm not sure about the explosion of these objects
            if asci in ['snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat', 'snIc_2002ap_V.dat']:
                ph = ph[:] + 1
            elif asci in ['snIa_K211b.dat']:
                ph = ph[1:]
                V = V[1:]
            elif asci in ['cv_sscyg_V.dat']:
                ph = ph[0:] - 12
                V = V[0:]

            models[asci] = {}
            models[asci]['days'] = ph
            models[asci]['mag'] = V - np.min(V)
    elif model == 'xxx':
        # in the case we want to use othe models
        pass
    return models
#########################################################

models1 = readmodels(model='shen')
models2 = readmodels(model='temp')
models3 = {}
models = {}
for key in models1:
    models[key] = models1[key]
for key in models2:
    models[key] = models2[key]
for key in models3:
    models[key] = models3[key]

_color = colorrange(len(models))
ii = 0

lab = {'snIa_K211b.dat': 'SN Ia',
     'snIc_2002ap_V.dat': 'SN Ic',
     'snIb_2009jf_V.dat': 'SN Ib',
     'snIb_ptf13bvn_V.dat': 'SN Ib',
     'PTF12cod_R.dat': 'SN II',
     'PTF11htj_R.dat': 'SN II',
     'PTF12bro_R.dat': 'SN II',
     'cv_tpyx_V.dat': 'CV',
     'cv_v1324sco_V.dat': 'CV',
     'cv_sscyg_V.dat': 'CV',
     'shen_0_6_0_2.dat': '.Ia',
     'shen_1_2_0_02.dat': '.Ia',
     'SNIa_CSM_r.dat': 'Ia+shock',
     'SNIa_RG1Msunprogenitor_maxinteraction_r.dat': 'Ia+CSM'}

for key in ['snIa_K211b.dat', 'snIc_2002ap_V.dat', 'snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat', 'PTF12cod_R.dat', 'PTF11htj_R.dat', \
            'PTF12bro_R.dat', 'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat', \
            'shen_0_6_0_2.dat', 'shen_1_2_0_02.dat', \
            'SNIa_CSM_r.dat']:
    ph = models[key]['days']
    mag = models[key]['mag']

    hours = .5
    sigma = 1
    ph1 = ph + hours / 24.
    mag1 = ninterp(ph1, ph, mag)
    diff1 = mag1 - mag

    hours = 2
    sigma = 1
    ph2 = ph + hours / 24.
    mag2 = ninterp(ph2, ph, mag)
    diff2 = mag2 - mag

    hours = 24
    sigma = 1
    ph3 = ph + hours / 24.
    mag3 = ninterp(ph3, ph, mag)
    diff3 = mag3 - mag


    if key in ['cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat']:
        tt = '-.'
        _lw = 4
    elif key in ['PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat']:
        tt = '--'
        _lw = 3
    elif key in models1:
        tt = ':'
        _lw = 2
    else:
        tt = '-'
        _lw = 1
    #ax2.plot(ph,mag,tt,color=_color[ii],label = key,lw=_lw)
    #ax2.plot(ph,mag,'o',color=_color[ii],label='')

    ax.plot(ph, diff1, tt, color=_color[ii], label=lab[key], lw=_lw)
    ax2.plot(ph, diff2, tt, color=_color[ii], label=lab[key], lw=_lw)
    ax3.plot(ph, diff3, tt, color=_color[ii], label=lab[key], lw=_lw)
    #    t = np.linspace(np.min(ph), np.max(ph), (int(np.max(ph))-int(np.min(ph))+1)*100)
    #    dy3 = np.array(mag)-np.array(mag) +.1
    #    flux=[]
    #    mvec=[]
    #    cvec=[]
    #    for mu in t:
    #        gaus1 = gaussian(ph,mu,sigma)
    #        gaus2 = gaus1/(dy3**2)
    #        m, c = np.polynomial.polynomial.polyfit(ph, mag, 1, w = [ty for ty in gaus2], full=False)
    #        flux.append(mu*c+m)
    #        mvec.append(m)
    #        cvec.append(c)

    ax.legend(ncol=1, loc=(.7, .1), numpoints=1, fontsize=10)
    ii = ii + 1

ax.set_xlim(0, 10)
ax.set_ylim(-.15, .05)
ax2.set_xlim(0, 10)
ax2.set_ylim(-.15, .05)
ax3.set_xlim(0, 10)
ax3.set_ylim(-1, .2)
ax3.set_ylabel('$\Delta$mag')
ax2.set_ylabel('$\Delta$mag')
ax.set_ylabel('$\Delta$mag')
#ax.set_ylabel('slope [mag/day]')
#ax2.set_ylabel('magnitude')
ax.set_title('30 minutes')
ax2.set_title('2 hours')
ax3.set_title('24 hours')
#ax2.set_title(' Light Curve')
#ax.set_title(' Slope')
ax.set_xlabel('Days after explosion')
ax2.set_xlabel('Days after explosion')
ax3.set_xlabel('Days after explosion')
ax.fill([0, 4, 4, 0], [-.02, -.02, -.01, -.01], 'g', alpha=0.1)
ax2.fill([0, 4, 4, 0], [-.08, -.08, -.03, -.03], 'g', alpha=0.1)
ax3.fill([0, 4, 4, 0], [-.9, -.9, -.3, -.3], 'g', alpha=0.1)
pl.show()
