#!/usr/bin/env python

import numpy as np
import pylab as pl
import glob
import os
import string
from numpy import interp as ninterp

pl.ion()
pl.rcParams['figure.figsize'] = 12, 5
ax = pl.axes([.08,.1,.41,.8])
ax2 = pl.axes([.55,.1,.41,.8])

def colorrange(N):
    N=N/3.
    color1={}
    for i in range(0,int(N)+1,1):
        color1[int(i)]=(0.0,i/N,1)
        color1[int(i+N)]=(i/N,1,1-i/N)
        color1[int(i+N*2)]=(1,1-i/N,0.)

    color2={}
    for j in color1.keys():
        print j
        color2[j]=color1[int(len(color1)-j-1)]
    return color2

########################################

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def readmodels(model='shen'):
    models={}
    if model == 'shen':
        _dir= 'models_shen2010/'
        files = glob.glob(_dir+'*2.dat')
        print files
        for file in files:
            data = np.genfromtxt(file,'float')
            t, M_bol,  M_U,  M_B, M_V, M_R,  M_I,  M_J,  M_H,  M_K = zip(*data)
            models[string.split(file,'/')[-1]] = {}
            models[string.split(file,'/')[-1]]['days'] = np.array(t)
            models[string.split(file,'/')[-1]]['mag'] = (np.array(M_V) * (-1)) + np.max(M_V)
    elif model == 'temp':
        _dir='SN_template/'
        for asci in ['snIa_K211b.dat','snIb_2009jf_V.dat','snIb_ptf13bvn_V.dat','PTF12cod_R.dat','PTF11htj_R.dat','PTF12bro_R.dat',\
                     'snIc_2002ap_V.dat','cv_tpyx_V.dat','cv_v1324sco_V.dat','cv_sscyg_V.dat']:
            data = np.genfromtxt(_dir + asci,'float')
            print asci
            if asci in ['snIa_K211b.dat','snIb_2009jf_V.dat','snIb_ptf13bvn_V.dat','PTF12cod_R.dat','PTF11htj_R.dat','PTF12bro_R.dat',\
                          'snIc_2002ap_V.dat','cv_tpyx_V.dat','cv_v1324sco_V.dat','cv_sscyg_V.dat']:
                ph,V = zip(*data)
            else:
                ph,u,B,V,r,i = zip(*data)
            ph = np.array(ph)
            V = np.array(V)
            # I'm not sure about the explosion of these objects
            if asci in ['snIb_2009jf_V.dat','snIb_ptf13bvn_V.dat','snIc_2002ap_V.dat']: 
                ph=ph[:]+1
            elif asci in ['snIa_K211b.dat']: 
                ph = ph[1:]
                V = V[1:]
            elif asci in ['cv_sscyg_V.dat']: 
                     ph = ph[0:]-12
                     V = V[0:]

            models[asci]={}
            models[asci]['days'] = ph
            models[asci]['mag'] = V - np.min(V)
    elif model == 'xxx':
        # in the case we want to use othe models
        pass
    return models

#########################################################

models1 = readmodels(model='shen')
models2 = readmodels(model='temp')
models3={}
models = {}
for key in models1:
    models[key]=models1[key]
for key in models2:
    models[key]=models2[key]
for key in models3:
    models[key]=models3[key]

_color = colorrange(len(models))
ii=0

lab={'snIa_K211b.dat':'SN Ia',
     'snIc_2002ap_V.dat':'SN Ic',
     'snIb_2009jf_V.dat':'SN Ib',
     'snIb_ptf13bvn_V.dat':'SN Ib',
     'PTF12cod_R.dat':'SN II',
     'PTF11htj_R.dat':'SN II',
     'PTF12bro_R.dat':'SN II',
     'cv_tpyx_V.dat':'CV',
     'cv_v1324sco_V.dat':'CV',
     'cv_sscyg_V.dat':'CV',
     'shen_0_6_0_2.dat':'.Ia', 
     'shen_1_2_0_02.dat':'.Ia'}

for key in ['snIa_K211b.dat','snIc_2002ap_V.dat','snIb_2009jf_V.dat','snIb_ptf13bvn_V.dat',\
            'PTF12cod_R.dat','PTF11htj_R.dat','PTF12bro_R.dat',\
            'cv_tpyx_V.dat','cv_v1324sco_V.dat','cv_sscyg_V.dat',\
            'shen_0_6_0_2.dat', 'shen_1_2_0_02.dat']:
    ph = models[key]['days']
    mag = models[key]['mag']
    
    hours = .5
    sigma = 1
    ph1 = ph + hours/24.
    mag1 = ninterp(ph1,ph,mag)
    diff1 = mag1-mag

    hours = 2
    sigma = 1
    ph2 = ph + hours/24.
    mag2 = ninterp(ph2,ph,mag)
    diff2 = mag2-mag

    hours = 24
    sigma = 1
    ph3 = ph + hours/24.
    mag3 = ninterp(ph3,ph,mag)
    diff3 = mag3-mag

    if key in ['cv_tpyx_V.dat','cv_v1324sco_V.dat','cv_sscyg_V.dat']:
        tt = '-.'
        _lw = 4
    elif key in ['PTF12cod_R.dat','PTF11htj_R.dat','PTF12bro_R.dat']: 
        tt = '--'
        _lw = 3
    elif key in models1:
         tt = ':'
         _lw = 2
    else:
        tt = '-'
        _lw = 1

    ax.plot(ph,mag,tt,color=_color[ii],label = '',lw=_lw)

    t = np.linspace(np.min(ph), np.max(ph), (int(np.max(ph))-int(np.min(ph))+1)*10)
    dy3 = np.array(mag)-np.array(mag) +.1
    flux=[]
    mvec=[]
    cvec=[]
    for mu in t:
        gaus1 = gaussian(ph,mu,sigma)
        gaus2 = gaus1/(dy3**2)
        m, c = np.polynomial.polynomial.polyfit(ph, mag, 1, w = [ty for ty in gaus2], full=False)
        flux.append(mu*c+m)
        mvec.append(m)
        cvec.append(c)

    ax2.plot(t,cvec,tt,color=_color[ii],label= lab[key],lw=_lw)
    ax2.legend(ncol=1,loc=(.8,.1),numpoints=1,fontsize = 10)
    ii=ii+1

ax2.set_ylim(-5,1)
ax.set_xlim(0,15)
ax2.set_xlim(0,15)
ax.set_ylim(8,-1)
ax2.set_ylabel('slope [mag/day]')
ax.set_ylabel('magnitude')
ax.set_title('Light curves')
ax2.set_title('Slope')
ax2.set_xlabel('Days after explosion')
ax.set_xlabel('Days after explosion')
