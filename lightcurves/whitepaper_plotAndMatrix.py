
#!/usr/bin/env python
from __future__ import print_function
__author__ = 'federica bianco'
### many chunks are taken from whitepaper_plot2.py written by stefano vlenti
import numpy as np
import pylab as pl
import glob
import os
import itertools
from numpy import interp as ninterp
import socket

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
alllcvs = ['snIa_K211b.dat', 'snIc_2002ap_V.dat',
           'snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat',
           'PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat',
           'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat', \
           'shen_0_6_0_2.dat', 'shen_1_2_0_02.dat', \
           'SNIa_CSM_r.dat']
    
#pl.ion()
#FBB if running on an interactive shell activate pl.ion() within the shel before execing the code
pl.rcParams['figure.figsize'] = 17, 6
#fig = pl.figure(frameon=True)
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
            models[f.split('/')[-1]]['mag'] = (np.array(M_V) * (-1)) +\
                                              np.max(M_V)
    elif model == 'temp':
        _dir = 'SN_template/'
        for asci in ['snIa_K211b.dat', 'snIb_2009jf_V.dat',
                     'snIb_ptf13bvn_V.dat',
                     'PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat', 
                     'snIc_2002ap_V.dat',
                     'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat', 
                     'SNIa_CSM_r.dat']:
            data = np.genfromtxt(_dir + asci, 'float')
            print ("this input is:", asci)
            if asci in ['snIa_K211b.dat', 'snIb_2009jf_V.dat',
                        'snIb_ptf13bvn_V.dat',
                        'PTF12cod_R.dat', 'PTF11htj_R.dat', 'PTF12bro_R.dat', 
                        'snIc_2002ap_V.dat',
                        'cv_tpyx_V.dat', 'cv_v1324sco_V.dat', 'cv_sscyg_V.dat',
                        'SNIa_CSM_r.dat']:
                
                ph, V = zip(*data)
            else:
                ph, u, B, V, r, i = zip(*data)
            ph = np.array(ph)
            V = np.array(V)
            # I'm not sure about the explosion of these objects
            if asci in ['snIb_2009jf_V.dat', 'snIb_ptf13bvn_V.dat',
                        'snIc_2002ap_V.dat']:
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


def calcDmag(time_lcv):
    hours, ph, mag = time_lcv
    # calculates the mag difference between 2 points of a lightcurve
    # defined by phase <ph> and magnitude <mag>
    # fiven a time interval <hours>
    sigma = 1
    mag1 = []
    hours = np.array(hours)
    
    ph1 = ph + hours / 24.
    mag1 = ninterp(ph1, ph, mag)
    return mag1 - mag


def setUp():
    models1 = readmodels(model='shen')
    models2 = readmodels(model='temp')
    models3 = {}
    models = {}
    for key in models1:
        print (key)
        models[key] = models1[key]
    for key in models2:
        models[key] = models2[key]
    for key in models3:
        models[key] = models3[key]
    return models


#########################################################


def wpplot2(gaps):

    print (gaps)
    models = setUp()
    _color = colorrange(len(models))
    fig = pl.figure()
    ax = [fig.add_subplot(1,len(gaps),i+1) for i in range(len(gaps))]
    for ii, key in enumerate(alllcvs):
        ph = models[key]['days']
        mag = models[key]['mag']
        
        diffs = map (calcDmag,
                     itertools.izip(gaps,
                                    itertools.repeat(ph),
                                    itertools.repeat(mag)))
        
        


        if 'cv_' in key:
            tt = '-.'
            _lw = 4
        elif 'PTF' in key:
            tt = '--'
            _lw = 3
        elif 'shen' in key:
            tt = ':'
            _lw = 2
        else:
            tt = '-'
            _lw = 1
        #ax2.plot(ph,mag,tt,color=_color[ii],label = key,lw=_lw)
        #ax2.plot(ph,mag,'o',color=_color[ii],label='')

        for i in range(len(gaps)):
            ax[i].plot(ph, diffs[i], tt, color=_color[ii], label=lab[key], lw=_lw)
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
        
        ax[0].legend(ncol=1, loc=4, numpoints=1, fontsize=10)

    for i,g in enumerate(gaps):
        ax[i].set_xlim(0, 10)
        if g <= 5:
            ax[i].set_ylim(-.15, .05)
        else:
            ax[i].set_ylim(-1, .2)
       
        ax[i].set_title('%.1f hours'%g)
        ax[i].set_xlabel('Days after explosion')
        if g == 0.5:
            ax[i].fill([0, 4, 4, 0], [-.02, -.02, -.01, -.01], 'g', alpha=0.1)
        elif g == 2:
            ax[i].fill([0, 4, 4, 0], [-.08, -.08, -.03, -.03], 'g', alpha=0.1)
        elif g == 24:
            ax[i].fill([0, 4, 4, 0], [-.9, -.9, -.3, -.3], 'g', alpha=0.1)
    ax[0].set_ylabel('$\Delta$mag')
    ax[-1].yaxis.set_label_position("right")
    ax[-1].set_ylabel('$\Delta$mag', rotation=270, labelpad=20)
     #pl.show()
    

def diffmatrix(gap, ax = None):
    # calculates confusion matrix for a given gap between observations
    # for events starting any 12 hours within 4 days

    #done cheaply in C style for loops. should be pythonized, when time allows
    
    models = setUp()
    Nm = len(models)
    _color = colorrange(Nm)
    Nobs = 4 # possible observations every 12 hours for 4 days
    Dobs = 12 #hours between events
    nEvents = np.arange(Nobs * 24.0/Dobs)
    nE = len(nEvents)
    print (nE)
    diffMatrix = np.zeros((Nm * nE, Nm * nE))
    xticks = []
    yticks = []    
    for ii, key0 in enumerate(alllcvs):
        [xticks.append(lab[key0]) if i==nE/2 else xticks.append('')\
         for i in range(nE)]
        print (xticks)
        yticks = []  
        for jj, key1 in enumerate(alllcvs):        
            [yticks.append(lab[key1]) if i==nE/2 else yticks.append('')\
             for i in range(nE)]
        
            ph0 = models[key0]['days']
            mag0 = models[key0]['mag']

            ph1 = models[key1]['days']
            mag1 = models[key1]['mag']
        
            diff0 = calcDmag((gap, ph0, mag0))
            diff1 = calcDmag((gap, ph1, mag1))
            for kk,t1 in enumerate(nEvents * Dobs):
                for ll,t2 in enumerate(nEvents * Dobs):
                   diffMatrix[ii*nE+kk][jj*nE+ll] =\
                                    np.abs(ninterp(t1,
                                                ph0, diff0) -\
                                        ninterp(t2,
                                                ph1, diff1))
    return diffMatrix, (xticks, yticks)


def plotDiffMatrix(diffMatrix, allticks, norm, gap, ax, log=True):
    # plots the confusion matrix as a similarity matrix in log scale
    # normalizing it as directed by the argument norm
    
    #print (diffMatrix)
    xticks, yticks = allticks

    # similarity matrix: normalized (max  - diff):
    # if diff=0 sim=1; if diff=max, sim=0
    diffMatrix = (norm[1] - diffMatrix) / (norm[1] - norm[0])

    
    if log: 
        diffMatrix = np.log(diffMatrix)
        #for ii in range(len(diffMatrix[0])):
        #    print (ii)
        #    diffMatrix[ii][ii] = diffMatrix[np.isfinite(diffMatrix)].min()-0.1
        #print(diffMatrix)

        thismin = (diffMatrix[np.isfinite(diffMatrix)].min())
    if not ax:
        ax = pl.figure().add_subplot(111)

    if log:
        extrema = [-0.1, 0]
    else:
        extrema = [0, 1.0]

    print (extrema)

    pl.imshow(diffMatrix, interpolation = 'nearest', cmap='bone',
              vmin = extrema[0],
              vmax = extrema[1])
    #cbar = pl.colorbar(ticks = np.linspace(-16,0.5, 5),
    #                   orientation='horizontal')
    cbar = pl.colorbar(ticks = extrema,
                       #np.linspace(diffMatrix[np.isfinite(diffMatrix)].min(),
                       #            diffMatrix[np.isfinite(diffMatrix)].max(),
                       #            5),
                       orientation='horizontal')
    cbar.set_ticks([extrema[0]*1.2, extrema[1]*0.8])
    cbar.set_ticklabels(['different', 'same'])    
    #cbar.set_label(r"$\Delta$mag similarity, log scale, day 0-4, 12h intervals", rotation=0)
    pl.xticks(range(diffMatrix.shape[0]), xticks, rotation=30)
    pl.yticks(range(diffMatrix.shape[0]), yticks, rotation=30)
    ax.set_title("%.1f hours"%gap)
#
#        diffMatrix[ii * nEvents] = calcDmag(hours, ph0, mag0)[ -\
#                                   calcDmag(hours, ph1, mag1)
  
#########################################################

if __name__ == '__main__':
    
    wpplot2([0.5, 2, 5, 24])
    fig = pl.figure()
    dm1, xt1 = diffmatrix(0.5)
    dm2, xt2 = diffmatrix(2)
    dm3, xt3 = diffmatrix(5)
    dm4, xt4 = diffmatrix(24)
    
    dmmin = min([dm1[np.isfinite(dm1)].min(),
                 dm2[np.isfinite(dm2)].min(),
                 dm3[np.isfinite(dm3)].min(),
                 dm4[np.isfinite(dm4)].min()])
    dmmax = max([dm1[np.isfinite(dm1)].max(),
                 dm2[np.isfinite(dm2)].max(),
                 dm3[np.isfinite(dm3)].max(),
                 dm4[np.isfinite(dm4)].max()])

    ax = fig.add_subplot(141)
    plotDiffMatrix(dm1, xt1, (dmmin, dmmax), 0.5, ax=ax, log=False)

    ax = fig.add_subplot(142)
    plotDiffMatrix(dm2, xt2, (dmmin, dmmax), 2, ax=ax, log=False)  

    ax = fig.add_subplot(143)
    plotDiffMatrix(dm3, xt3, (dmmin, dmmax), 5, ax=ax, log=False)  

    ax = fig.add_subplot(144)
    plotDiffMatrix(dm4, xt4, (dmmin, dmmax), 24, ax=ax, log=False)
    #pl.show()

    fig = pl.figure()
 
    ax = fig.add_subplot(141)
    plotDiffMatrix(dm1, xt1, (dmmin, dmmax), 0.5, ax=ax, log=True)

    ax = fig.add_subplot(142)
    plotDiffMatrix(dm2, xt2, (dmmin, dmmax), 2, ax=ax, log=True)  

    ax = fig.add_subplot(143)
    plotDiffMatrix(dm3, xt3, (dmmin, dmmax), 5, ax=ax, log=True)  

    ax = fig.add_subplot(144)
    plotDiffMatrix(dm4, xt4, (dmmin, dmmax), 24, ax=ax, log=True)
    pl.show()
