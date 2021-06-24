import os
import re
from astropy.io import fits 
import scipy
import scipy.signal
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


include_path='/home/simon/common/python/include/'
sys.path.append(include_path)
import ImUtils.Resamp as Resamp
import ImUtils.Cube2Im as Cube2Im


def spider(filetag='',r1=0.3,r2=3.,nrlevs=20,nphilevs=40,fileout='spider.pdf'):

    nplotsx=1
    nplotsy=1
    iplotpos=1
    
    ax = plt.subplot(nplotsy, nplotsx, iplotpos)
    
    VisibleYaxis=False
    VisibleXaxis=False
    
    plt.setp(ax.get_xticklabels(),visible=VisibleXaxis)
    plt.setp(ax.get_yticklabels(),visible=VisibleYaxis)
    
    ax.tick_params(axis='both',length = 5, width=1., color = 'grey',direction='in',left=True, right=True,bottom=True, top=True)
    
    ax.spines['right'].set_color('grey')
    ax.spines['left'].set_color('grey')
    ax.spines['top'].set_color('grey')
    ax.spines['bottom'].set_color('grey')
    
    filename_grey=filetag+'rrs_top_sky.fits'
    f = fits.open(filename_grey)
    im_rrs = f[0].data
    hdr_rrs= f[0].header
    cdelt=3600.*hdr_rrs['CDELT2']
    side=hdr_rrs['NAXIS2']*cdelt

    filename_grey=filetag+'phis_top_sky.fits'
    f = fits.open(filename_grey)
    im_phis = f[0].data
    hdr_phis= f[0].header


    a0 = side/2.
    a1 = -side/2.
    d0 = -side/2.
    d1 = side/2.

    clevs=(r2-r1)*np.arange(nrlevs)/(nrlevs-1) + r1 
                    
    ax.contour(im_rrs, clevs , origin='lower', linewidths=1.0,
               linestyles = 'solid', 
               extent=[a0,a1,d0,d1], colors='grey')


    clevs=(2.*np.pi)*np.arange(nphilevs)/(nphilevs-1) 
                    
    ax.contour(im_phis, clevs , origin='lower', linewidths=1.0,
               linestyles = 'solid', 
               extent=[a0,a1,d0,d1], colors='grey')

    
    plt.subplots_adjust(hspace=0.)
    plt.subplots_adjust(wspace=0.)
        
    fileout=filetag+fileout
    print( fileout)
    plt.savefig(fileout, bbox_inches='tight') # dpi

