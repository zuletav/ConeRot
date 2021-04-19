import sys
import numpy as np
import scipy  as sp
import os
import os.path
from scipy import ndimage
from astropy.io import fits 
import re
from copy import deepcopy
from astropy.wcs import WCS
from scipy import optimize
import time
from time import gmtime,strftime
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from pylab import *
#import matplotlib.colors as colors

include_path='/home/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import gridding




#import MPolarMaps.Master

import PyVtools.Vtools as Vtools

def func(x, a, b, c):
    return a * np.exp(-b * x) + c



def z_func(r,z0,r0,q):
    taper=(1.-np.tanh(r/r0))
    return (z0 * (r/r0)**q)*taper


def ftanpsi(r,z0,r0,q):
    z=z_func(r,z0,r0,q)
    return z/r

def cartesian2conicpolar(outcoords, inputshape, origin, inc=0.,tanpsi=0.):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    rindex, thetaindex = outcoords
    x0, y0 = origin

    
    theta = float(thetaindex) * 2. * np.pi / float(inputshape[0])
    
    y = rindex*np.cos(theta)

    height = tanpsi * rindex
    #height = -tanpsi * rindex
    x = (rindex*np.sin(theta))/np.cos(inc) + \
             (height - (rindex*np.sin(theta))*np.tan(inc))*np.sin(inc)
        
    ix = -x + x0
    iy = y + float(y0)
    
    return (iy,ix)
 


def cartesian2polar(outcoords, inputshape, origin, inc=0.):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""


    rindex, thetaindex = outcoords
    x0, y0 = origin

    
    theta = float(thetaindex) * 2. * np.pi / float(inputshape[0]-1.)
    
    y = rindex*np.cos(theta)

    x = rindex*np.sin(theta)*np.cos(inc) 
        
    ix = -x + x0
    iy = y + float(y0)
    
    return (iy,ix)
 



def conicpolar2cartesian_ellipse(outcoords, inputshape, origin,inc=0.,tanpsi=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx=inputshape[0]
    ny=inputshape[1]
    x = -float(xindex - x0)
    y = float(yindex - y0)

    tanpsi0=tanpsi


    a=((np.tan(inc) * tanpsi0)**2-1.0)
    b=-2.*x*np.sin(inc)* tanpsi0/(np.cos(inc))**2
    c=y**2+(x**2/(np.cos(inc)**2))
    Delta=b**2-4.*a*c
    rho=(-b-np.sqrt(Delta))/(2.*a)
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

        
    #sintheta=np.sqrt( 1.- costheta**2)
    #if (x<0):
    #    sintheta*=-1
    
    #xp=(rho*sintheta/np.cos(inc)) + (tanpsi0*rho  - rho*sintheta*np.tan(inc)) * np.sin(inc)

    H1=tanpsi0*rho
    num= x - H1 * np.sin(inc)
    denom= rho * ( (1./np.cos(inc))  - np.tan(inc) * np.sin(inc))
    sintheta= num/ denom
    
    theta=np.arccos(costheta)
    
    if sintheta<0:
        theta = 2.*np.pi - theta


    #theta = 2.*np.pi - theta

    #slope=np.tan(inc)#*tanpsi0

    #if (x<0.):
    #    if (x<slope*y):
    #        #print("theta",theta,"x",x,"y",y,"costheta",costheta)


    #if (x<0):
    #    #print("theta",theta,"x",x,"y",y)
    #    theta += np.pi 
    
    thetaindex = (theta * (float(nx)-1.) / (2. * np.pi)) 
    #thetaindex = (theta * float(nx) / (np.pi))


            
    return (rindex,thetaindex)



def polar2cartesian(outcoords, inputshape, origin,inc=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx=inputshape[0]-1
    ny=inputshape[1]-1
    x = -float(xindex - x0)
    y = float(yindex - y0)

    rho=np.sqrt((x**2/(np.cos(inc)**2))+y**2)
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

    # theta=np.arccos(costheta)
    theta = np.arctan2((-x/np.cos(inc)), y) 
    if (theta < 0):
        theta = theta + 2.*np.pi

    thetaindex = (theta * float(nx) / (2. * np.pi))
            
    return (rindex,thetaindex)




# def exec_flatpolar(filename_fullim,hdr2,aresampim,PA,cosi): 
#     rotangle= PA - 180.
#     #    rotangle= PA
#     im1rot = ndimage.rotate(aresampim, rotangle, reshape=False)
#     fileout_rotated=re.sub('fullim.fits', 'rotated.fits', filename_fullim)
#     pf.writeto(fileout_rotated,im1rot, hdr2, overwrite=True)
#     
#     
#     hdr3 = deepcopy(hdr2)
#     hdr3['CDELT1']=hdr3['CDELT1']*cosi
#     fileout_stretched=re.sub('fullim.fits', 'stretched.fits', filename_fullim)
#     im3=np.double(gridding(fileout_rotated,hdr3))
#     pf.writeto(fileout_stretched,im3, hdr2, overwrite=True)
# 
# 
#     nx=hdr2['NAXIS1']
#     ny=hdr2['NAXIS2']
# 
#     aim_polar = sp.ndimage.geometric_transform(im3,cartesian2polar, 
#                                                order=1,
#                                                output_shape = (im3.shape[0], im3.shape[1]),
#                                                extra_keywords = {'inputshape':im3.shape,
#                                                                  'origin':(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)}) 
#     
#     return aim_polar






def carttopolar(im,inc):

    (ny,nx)=im.shape
    (i0,j0)=(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)

    im_polar = sp.ndimage.geometric_transform(im,cartesian2polar, 
                                              order=1,
                                              output_shape = (im.shape[0], im.shape[1]),
                                              extra_keywords = {
                                                  'inputshape':im.shape,
                                                  'inc':inc,
                                                  'origin':(i0,j0)})
    return im_polar


def conicpolartocart(im_polar,inc,tanpsi):

    (ny,nx)=im_polar.shape
    (i0,j0)=(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)

    
    im_cart = sp.ndimage.geometric_transform(im_polar,conicpolar2cartesian_ellipse,
                                            order=1,
                                            output_shape = (im_polar.shape[0], im_polar.shape[1]),
                                            extra_keywords = {'inputshape':im_polar.shape,
                                                  'inc':inc,'tanpsi':tanpsi,
                                                  'origin':(i0,j0)})
    im_cart=np.nan_to_num(im_cart)

    return im_cart



def polartocart(im_polar,inc):

    (ny,nx)=im_polar.shape
    (i0,j0)=(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)

    
    im_cart = sp.ndimage.geometric_transform(im_polar,polar2cartesian,
                                            order=1,
                                            output_shape = (im_polar.shape[0], im_polar.shape[1]),
                                            extra_keywords = {'inputshape':im_polar.shape,
                                                  'inc':inc,
                                                  'origin':(i0,j0)})
    im_cart=np.nan_to_num(im_cart)

    return im_cart





def gen_surface(file_psiprofile,file_canvas,PA=0.,inc=0.,fileouttag='H',ForceTop=False):

    (rregions, psis, psi_ds, psi_us) = np.loadtxt(file_psiprofile,unpack=True)

    Hs=rregions*np.tan(psis*np.pi/180.)
    tanpsis=np.tan(psis*np.pi/180.)
    rmax0=np.max(rregions)

    if (ForceTop):
        Hs=np.fabs(Hs)
        tanpsis=np.fabs(tanpsis)

    f0 = fits.open(file_canvas)
    im_canvas0 = f0[0].data
    hdr_canvas0= f0[0].header

    #Vtools.View(f0)


    zoomfactor=2.
    hdr_canvas = deepcopy(hdr_canvas0)
    hdr_canvas['CDELT1']*=zoomfactor
    hdr_canvas['CDELT2']*=zoomfactor

    f1=gridding(f0,hdr_canvas,ReturnHDUList=True)
    im_canvas = f1[0].data
    


    #Vtools.View(im_canvas)
    

    cdelt=hdr_canvas['CDELT2']*3600.
    
    (ny,nx) = im_canvas.shape 
    x=np.arange(1,nx+1)
    y=np.arange(1,ny+1)
    X, Y = np.meshgrid(x, y)
    
    X0 = np.floor(nx/2)+1
    Y0 = np.floor(ny/2)+1
    dxxs = -cdelt *(X-X0)
    dyys = cdelt *(Y-Y0)
    rrs = np.sqrt( (dxxs)**2 + (dyys)**2)

    rmax=np.max(rrs)
    print("rmax",rmax)
    if (zoomfactor > 1.):
        rregions=np.append([0.],rregions)
        rregions=np.append(rregions,[rmax/2.,rmax])    
        Hs=np.append(0.,Hs)
        Hs=np.append(Hs,[0.,0.])
        tanpsis=np.append(0.,tanpsis)
        tanpsis=np.append(tanpsis,[0.,0.])
        
    else:
            
        rregions=np.append([0.],rregions)
        rregions=np.append(rregions,rmax)    
        Hs=np.append(0.,Hs)
        Hs=np.append(Hs,0.)
        tanpsis=np.append(0.,tanpsis)
        tanpsis=np.append(tanpsis,0.)

        
    print("rregions",rregions)
    print("Hs",Hs)
    

    #fH = interp1d(rregions, Hs, kind='cubic')
    #ftanpsi = interp1d(rregions, tanpsis, kind='cubic')
    ##finvH = interp1d(Hs, rregions,  kind='cubic')
    

    popt, pcov = curve_fit(z_func, rregions, tanpsis, p0=[0.,0.5,1.],bounds=[[-10.,0.,0.1],[10.,10.,4.]])

    #popt, pcov = curve_fit(z_func, rregions, tanpsis, p0=[0.3,0.5,1.])
    
    #popt, pcov = curve_fit(z_func, rregions, tanpsis, p0=[0.,0.5,1.])

    obsprof=np.zeros((len(rregions),2))
    modprof=np.zeros((len(rregions),2))
    obsprof[:,0]=rregions
    obsprof[:,1]=tanpsis

    
    z0=popt[0]
    r0=popt[1]
    q=popt[2]
    print("r0 ",r0,"z0 ",z0, "q", q)

    Retro=False
    H_sign=+1
    
    #if z0 < 0.:
    #    Hsign=-1
    #    Retro=True
        
    modprof[:,0]=rregions
    modprof[:,1]=z_func(rregions,z0,r0,q)
    
    #Vtools.Spec([obsprof,modprof])

    #HHs=fH(rrs)
    HHs=H_sign*z_func(rrs,z0,r0,q)
    
    fileout=fileouttag+'_faceon.fits'
    fits.writeto(fileout,HHs, hdr_canvas, overwrite=True)


    
    ########################################################################
    ## default expansion
    #
    #M.workdir='polarmaps_default/'  # directory for products 
    #M.prep_files()
    #M.polar_expansions()

    nphis=HHs.shape[0]
    print("nphis",nphis)
    nrs=HHs.shape[1]
    print("nrs",nrs)
    
    HHs_polar=carttopolar(HHs,0.)
    hdupolar = fits.PrimaryHDU()
    hdrpolar=hdupolar.header
    hdrpolar['NAXIS1']=nrs
    hdrpolar['NAXIS2']=nphis
    hdrpolar['CRPIX1']=1
    hdrpolar['CRVAL1']=0.
    hdrpolar['CDELT1']=2. * np.pi / nphis
    hdrpolar['CRPIX2']=1
    hdrpolar['CRVAL2']=0.
    hdrpolar['CDELT2']=hdr_canvas['CDELT2']
    hdupolar.header = hdrpolar

    rs =  3600.*(np.arange(hdrpolar['NAXIS2'])-hdrpolar['CRPIX2']+1.0)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']

    rrs_polar=np.zeros(HHs_polar.shape)

    print("nphis",nphis,"nx",nx)

    x=np.arange(1,nx+1)
    y=np.arange(1,ny+1)
    X, Y = np.meshgrid(x, y)
    
    rrs_polar= 3600.*(Y-hdrpolar['CRPIX2']+1.0)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']

    #Vtools.View(rrs_polar)
        
    phis_polar= (X-hdrpolar['CRPIX1']+1.0)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']
                                
    
    #Vtools.View(phis_polar)

    
    HHs_sky_top=np.zeros(HHs.shape)
    HHs_sky_bottom=np.zeros(HHs.shape)

    phis_sky_top=np.zeros(HHs.shape)


    region_domain_top_prev=np.zeros(rrs_polar.shape)
    region_domain_bottom_prev=np.zeros(rrs_polar.shape)

    nmesh=10
    rmesh = (np.arange(nmesh)/(nmesh-1))*rmax

    
    iR2=nmesh-1
    for iregion in list(range(iR2)):
        R1=rmesh[iregion]
        R2=rmesh[iregion+1]
        print("R1 ", R1, "R2: ",R2)
        tanpsi=H_sign*ftanpsi(R2,z0,r0,q)

        region_polar=np.zeros(rrs_polar.shape)
        #if (iregion == (iR2-1)):
        #    print("LAST")
        #    region_polar[(rrs_polar >= R1)]=1.
        #else:
        #    #region_polar[(rrs_polar >= R1) & ( rrs_polar < R2)]=1.
        #    region_polar[( rrs_polar < R2)]=1.

        region_polar[( rrs_polar < R2)]=1.

        print("calling  conicpolartocart with ", inc,tanpsi)

        #plt.imshow(region_polar)
        #plt.show()
        phis_top=conicpolartocart(phis_polar,inc,tanpsi)
        #phis_bottom=conicpolartocart(phis_polar,inc,-tanpsi)
        #Vtools.View(phis_top)
        
        region_domain_top=conicpolartocart(region_polar,inc,tanpsi)
        region_domain_bottom=conicpolartocart(region_polar,inc,-tanpsi)
        

        

        
        #plt.imshow(region_domain_top)
        #plt.show()
        #plt.imshow(phis_top)
        #plt.show()
            
        HHs_sky_region_top=conicpolartocart(HHs_polar,inc,tanpsi)
        HHs_sky_region_bottom=-conicpolartocart(HHs_polar,inc,-tanpsi)
        
        #Vtools.View(HHs_sky_region_top)

        #HHs_sky_top[(region_domain_top > 0.5)]= HHs_sky_region_top[(region_domain_top > 0.5)]
        #HHs_sky_bottom[(region_domain_bottom > 0.5)]= HHs_sky_region_bottom[(region_domain_bottom > 0.5)]

        #nearside= ( (phis_top >= 0.) & (phis_top <= np.pi ))
        #farside=  ( (phis_top < 0.) | (phis_top > np.pi ))
        
        mask= ( (region_domain_top >= 0.5) & (region_domain_top_prev < 0.5) )
        HHs_sky_top[mask]= HHs_sky_region_top[mask]
        phis_sky_top[mask]=phis_top[mask]
        
        #Vtools.View(HHs_sky_top)

        mask= ( (region_domain_bottom > 0.5) & (region_domain_bottom_prev < 0.5))
        HHs_sky_bottom[mask]= HHs_sky_region_bottom[mask]

        region_domain_top_prev=region_domain_top
        region_domain_bottom_prev=region_domain_bottom
        
        #Vtools.View(HHs_sky_top)


    rotangle= -PA

    HHs_sky_rot_top = ndimage.rotate(HHs_sky_top, rotangle, reshape=False)
    phis_sky_rot_top = ndimage.rotate(phis_sky_top, rotangle, reshape=False)
    HHs_sky_rot_bottom = ndimage.rotate(HHs_sky_bottom, rotangle, reshape=False)

    hdutop = fits.PrimaryHDU()
    hdutop.header=hdr_canvas
    hdutop.data=HHs_sky_rot_top
    hdutop0=gridding(hdutop,hdr_canvas0,ReturnHDUList=True)
    fileout=fileouttag+'_top_sky.fits'
    hdutop0.writeto(fileout, overwrite=True)

    hdutop = fits.PrimaryHDU()
    hdutop.header=hdr_canvas
    hdutop.data=phis_sky_rot_top
    hdutop0=gridding(hdutop,hdr_canvas0,ReturnHDUList=True)
    fileout='phis_top_sky.fits'
    hdutop0.writeto(fileout, overwrite=True)

    hdubottom = fits.PrimaryHDU()
    hdubottom.header=hdr_canvas
    hdubottom.data=HHs_sky_rot_bottom
    hdubottom0=gridding(hdubottom,hdr_canvas0,ReturnHDUList=True)
    fileout=fileouttag+'_bottom_sky.fits'
    hdubottom0.writeto(fileout, overwrite=True)

    ##fits.writeto(fileout,HHs_sky_rot_top, hdr_canvas, overwrite=True)
    #fileout=fileouttag+'_bottom_sky.fits'
    #fits.writeto(fileout,HHs_sky_rot_bottom, hdr_canvas, overwrite=True)
    #fileout=fileouttag+'_dif_sky.fits'
    #fits.writeto(fileout,(HHs_sky_rot_top-HHs_sky_rot_bottom), hdr_canvas, overwrite=True)

    
    return


