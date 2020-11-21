import sys
import numpy as np
import scipy  as sp
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
from copy import deepcopy
from astropy.wcs import WCS
from scipy import optimize
import time
from time import gmtime,strftime

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.colors as colors

include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import gridding


if not sys.warnoptions:
    import os, warnings
    #warnings.simplefilter("default") # Change the filter in this process
    warnings.simplefilter("ignore") # Change the filter in this process
    #os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses
    os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

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
        
    theta=np.arccos(costheta)
    thetaindex = (theta * float(nx) / (2. * np.pi))


            
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




def carttoconicpolar(im,inc,tanpsi):

    (ny,nx)=im.shape
    (i0,j0)=(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)

    im_polar = sp.ndimage.geometric_transform(im,cartesian2conicpolar, 
                                              order=1,
                                              output_shape = (im.shape[0], im.shape[1]),
                                              extra_keywords = {
                                                  'inputshape':im.shape,
                                                  'inc':inc,'tanpsi':tanpsi,
                                                  'origin':(i0,j0)})
    return im_polar


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


def exec_prep_files(M):
    
    filename_source=M.filename_source
    workdir=M.workdir
    DumpAllFitsFiles=M.DumpAllFitsFiles

    if (not re.search(r"\/$",workdir)):
        workdir+='/'
        M.workdir=workdir
        print("added trailing back slash to outputdir")



    inbasename=os.path.basename(filename_source)
    filename_fullim=re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim=workdir+filename_fullim
        
    
    if (M.DoErrorMap):
        inbasenameerr=os.path.basename(M.filename_errormap)
        filename_fullimerr=re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr=workdir+filename_fullimerr
        filename_fullimw=re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw=workdir+filename_fullimw        

    if (M.Verbose): #
        print( "BUILDING WORKDIR AND GENERATING CENTERED IMAGE")

    hdu = pf.open(filename_source)
    hdr0=hdu[0].header
    if (hdr0['NAXIS'] > 2):
        if (M.Verbose): #
            print( "running cube2im")
        hdu=cube2im(filename_source,False)
        
    im1=hdu[0].data
    im1=im1*M.unitscale
    print( "applied unit scale factor:",M.unitscale)
    hdr1=hdu[0].header
        
    typicalerror=M.typicalerror
    if (M.InjectNoise):
        print( "INJECTING NOISE")
        im1=im1 + np.random.normal(loc=0.0, scale=typicalerror, size=im1.shape)


    if (DumpAllFitsFiles):
        pf.writeto(filename_fullim,im1, hdr1, overwrite=True)
    
        
    hdu[0].data=im1
    
    hduw=False
    
    if (M.DoErrorMap):
        hduerr = pf.open(M.filename_errormap)
        hdrerr = hduerr[0].header
        if (hdrerr['NAXIS'] > 2):
            hduerr=cube2im(M.filename_errormap,False)


        imerr1=hduerr[0].data
        imerr1=imerr1*M.unitscale

        #typicalerror=np.median(imerr1)
        #imerr1[np.where(imerr1 < (typicalerror/3.))] = typicalerror
        #imerr1[np.where(imerr1 > (100.*typicalerror))] = 1E20



        hdrerr1=hduerr[0].header
            

        imw0=1./imerr1**2
        imw0=np.nan_to_num(imw0)
        imw0[np.where(np.fabs(imw0) > 1E10)] = 0.
        
        if (M.Verbose):
            print( "fullim weights: max",np.max(imw0)," min ",np.min(imw0))

        imerr=np.sqrt(1./imw0)
        imerr=np.nan_to_num(imerr)
        imerr[np.where(np.fabs(imw0) < 1E-30)] = 1E20
        
        typicalerror=np.median(imerr)

        if (M.Verbose):
            print( "resamp errors: max",np.max(imerr)," min ",np.min(imerr))
            print( "typicalerror= ",typicalerror," vs", M.typicalerror)

            
        if (DumpAllFitsFiles):
            pf.writeto(filename_fullimw,imw0, hdrerr1, overwrite=True)

        hduw=pf.PrimaryHDU()

        hduw.data=imw0
        hduw.header=hdrerr1
            


    M.Hdu=hdu
    M.Hduw=hduw
    
    return



def exec_grid_4center(M):

    workdir=M.workdir
    RA=M.RA
    DEC=M.DEC
    x_center=M.x_center
    y_center=M.y_center
    fieldscale=M.fieldscale # shrink radial field of view of polar maps by this factor

    hdu=M.Hdu
    hduw=M.Hduw

    
    filename_source=M.filename_source

    inbasename=os.path.basename(filename_source)
    filename_fullim=re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim=workdir+filename_fullim
    fileout_centered=re.sub('fullim.fits', 'centered.fits', filename_fullim)

    if (M.DoErrorMap):
        inbasenameerr=os.path.basename(M.filename_errormap)
        filename_fullimerr=re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr=workdir+filename_fullimerr
        filename_fullimw=re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw=workdir+filename_fullimw        
        fileout_centerederr=re.sub('fullim.fits', 'centered.fits', filename_fullimerr)
        fileout_centeredw=re.sub('fullim.fits', 'wcentered.fits', filename_fullimerr)


    hdr1=hdu[0].header
        
    if (not (isinstance(y_center,bool))):
        if (not isinstance(RA,float)):
            RA=hdr1['CRVAL1']
            DEC=hdr1['CRVAL2']
            if (M.Verbose):
                print( "using center of coords CRVAL1 CRVAL2")
                            
        #RA=RA+(np.sin(x_center*np.pi/180.)*y_center/3600.)/np.cos(DEC*np.pi/180.)
        RA=RA+((x_center/3600.)/np.cos(DEC*np.pi/180.))
        #DEC=DEC+np.cos(x_center*np.pi/180.)*y_center/3600.
        DEC=DEC+((y_center/3600.)*np.pi/180.)
        if (M.Verbose):
            print( "RA =", RA)
            print( "DEC =",DEC)



    nx=int(hdr1['NAXIS1']/(M.pixscale_factor*M.fieldscale))
    ny=nx


    if ( (nx % 2) == 0):
        nx=nx+1
        ny=ny+1
        
    hdr2 = deepcopy(hdr1)

    hdr2['NAXIS1']=nx
    hdr2['NAXIS2']=ny
    hdr2['CRPIX1']=(nx+1)/2
    hdr2['CRPIX2']=(ny+1)/2
    hdr2['CRVAL1']=RA
    hdr2['CRVAL2']=DEC
    hdr2['CDELT1']=M.pixscale_factor*hdr2['CDELT1']
    hdr2['CDELT2']=M.pixscale_factor*hdr2['CDELT2']

        
    resamp=gridding(hdu,hdr2, fullWCS=False)
    resamp=np.nan_to_num(resamp)

        
    if (M.DumpAllFitsFiles):
        fileout_centered=re.sub('fullim.fits', 'centered.fits', filename_fullim)
        pf.writeto(fileout_centered,resamp, hdr2, overwrite=True)


    
    hducentered=pf.PrimaryHDU()
    hducentered.data=resamp
    hducentered.header=hdr2

    hduwcentered=False
    
    if (M.DoErrorMap):
            resampw=gridding(hduw,hdr2, fullWCS=False)
            resampw=np.nan_to_num(resampw)
            resampw[np.where(resampw < 0.)] =0.

            if (M.Verbose): #
                print( "resamp weights: max",np.max(resampw)," min ",np.min(resampw))

            resamperr=np.sqrt(1./resampw)
                              
            if (M.DumpAllFitsFiles):
                fileout_centeredw=re.sub('fullim.fits', 'wcentered.fits', filename_fullimerr)
                pf.writeto(fileout_centeredw,resampw, hdr2, overwrite=True)
                fileout_centerederr=re.sub('fullim.fits', 'centered.fits', filename_fullimerr)
                pf.writeto(fileout_centerederr,resamperr, hdr2, overwrite=True)

            

            hduwcentered=pf.PrimaryHDU()
            hduwcentered.data=resampw
            hduwcentered.header=hdr2


            
    M.Ncorr=(np.pi/(4.*np.log(2)))*M.bmaj*M.bmin/(hdr2['CDELT2']*3600.)**2

    M.Hducentered=hducentered
    M.Hduwcentered=hduwcentered
    

 

def exec_conicpolar_expansions(M):

    
    filename_source=M.filename_source
    workdir=M.workdir
    inc=M.inc
    tanpsi=M.tanpsi
    RA=M.RA
    DEC=M.DEC
    x_center=M.x_center
    y_center=M.y_center
    DoAzimuthalProfile=M.DoAzimuthalProfile
    PlotRadialProfile=M.PlotRadialProfile
    a_min=M.a_min
    a_max=M.a_max
    typicalerror=M.typicalerror
    InjectNoise=M.InjectNoise
    DumpAllFitsFiles=M.DumpAllFitsFiles
    DoDCone=M.DoDCone
    domain=M.domain
    fieldscale=M.fieldscale # shrink radial field of view of polar maps by this factor
    DoAccr=M.DoAccr
    DoMerid=M.DoMerid
    
    hdu=M.Hducentered
    hduw=M.Hduwcentered

    #typicalerror=ExpectedError


    
    #cosi=np.cos(inc*np.pi/ 180.)


    inbasename=os.path.basename(filename_source)
    filename_fullim=re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim=workdir+filename_fullim
    fileout_centered=re.sub('fullim.fits', 'centered.fits', filename_fullim)


    if (M.DoErrorMap):
        inbasenameerr=os.path.basename(M.filename_errormap)
        filename_fullimerr=re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr=workdir+filename_fullimerr
        filename_fullimw=re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw=workdir+filename_fullimw        
        fileout_centerederr=re.sub('fullim.fits', 'centered.fits', filename_fullimerr)
        fileout_centeredw=re.sub('fullim.fits', 'wcentered.fits', filename_fullimerr)


    resamp=hdu.data
    hdr2=hdu.header
    nx=hdr2['NAXIS1']
    ny=hdr2['NAXIS2']

    if M.Verbose:
        print( "M.Ncorr = ",M.Ncorr)

    if (M.DoErrorMap):
        resampw=hduw.data



    if (M.InheritMumap):
        if (M.mumap is None):
            mumap=np.ones(resamp.shape)
        else:
            mumap=M.mumap
    else:
            mumap=np.ones(resamp.shape)
            
    
    #rotangle= M.PA - 180.
    rotangle= M.PA

    im1rot = ndimage.rotate(resamp, rotangle, reshape=False)
    hdr3 = deepcopy(hdr2)
    im3=np.double(im1rot) 

    if (DumpAllFitsFiles):
        fileout_rotated=re.sub('fullim.fits', 'rotated.fits', filename_fullim)
        pf.writeto(fileout_rotated,im1rot, hdr2, overwrite=True)

        

    if (M.DoErrorMap):
        if (np.any(resampw < 0.)):
            print("min / max:",np.min(resampw),np.max(resampw))
            sys.exit("negative  sky weights!!!!")

        im1rotw = ndimage.rotate(resampw, rotangle, reshape=False,order=0)


        im3w=np.double(im1rotw)
        if (DumpAllFitsFiles):
            fileout_rotatedw=re.sub('fullim.fits', 'wrotated.fits', filename_fullimerr)
            pf.writeto(fileout_rotatedw,im1rotw, hdr2, overwrite=True)

            if (np.any(im1rotw < 0.)):
                print("min / max:",np.min(im1rotw),np.max(im1rotw))
                sys.exit("negative rot sky weights!!!!")
    
        
    # #####################################################################
    # take conic polar transforms


    if (M.Verbose):
        print( "CARTESIAN2CONICPOLAR TRANSFORM START")
        print( "using inc ",inc*np.pi/180.," tanpsi ",tanpsi)

    im_polar = carttoconicpolar(im3,inc,tanpsi)

    if (DoDCone):
        mumap_polarpos = carttoconicpolar(mumap,inc,tanpsi)
        
    nphis,nrs=im_polar.shape
    if ((nphis != nx) or (nrs != ny)):
        sys.exit("bug")

    hdupolar = pf.PrimaryHDU()
    hdupolar.data = im_polar
    hdrpolar=hdupolar.header
    hdrpolar['CRPIX1']=1
    hdrpolar['CRVAL1']=0.
    hdrpolar['CDELT1']=2. * np.pi / nphis
    hdrpolar['CRPIX2']=1
    hdrpolar['CRVAL2']=0.
    hdrpolar['CDELT2']=(hdr3['CDELT2'])
    hdupolar.header = hdrpolar
    
    fileout_polar=re.sub('fullim.fits', 'polar.fits', filename_fullim)
    if (DumpAllFitsFiles):
        hdupolar.writeto(fileout_polar, overwrite=True)
            

    if (M.DoErrorMap):
        

        im_polarw = carttoconicpolar(im3w,inc,tanpsi)

        nphis,nrs=im_polarw.shape
        
        hdupolarw = pf.PrimaryHDU()
        hdupolarw.data = im_polarw
        hdupolarw.header = hdrpolar
        
        fileout_polarw=re.sub('fullim.fits', 'wpolar.fits', filename_fullimerr)
        if (DumpAllFitsFiles):
            hdupolarw.writeto(fileout_polarw, overwrite=True)
    else:
        im_polarw=np.ones(im_polar.shape)/typicalerror**2
        


    ######################################################################
    # take azimuthal averages on polar maps


    weights=im_polarw.copy()
    im_Npolcorr=np.ones(im_polarw.shape)

    if (np.any(weights < 0.)):
        print("min / max:",np.min(weights),np.max(weights))
        sys.exit("negative polarweights!!")
     
    im_polar_av = np.copy(im_polar)

    rrs =  3600.*(np.arange(hdrpolar['NAXIS2'])-hdrpolar['CRPIX2']+1.0)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']
    phis =  (180./np.pi)*((np.arange(hdrpolar['NAXIS1'])-hdrpolar['CRPIX1']+1.0)*hdrpolar['CDELT1']+hdrpolar['CRVAL1'])
    phis_rad =  np.double(((np.arange(hdrpolar['NAXIS1'])-hdrpolar['CRPIX1']+1.0)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']))
    KepAmps= np.double(np.zeros(len(rrs)))
    sKepAmps= np.double(np.zeros(len(rrs)))
    AccrAmps= np.double(np.zeros(len(rrs)))
    sAccrAmps= np.double(np.zeros(len(rrs)))
    MeridAmps= np.double(np.zeros(len(rrs)))
    sMeridAmps= np.double(np.zeros(len(rrs)))

    vsysts=np.zeros(hdrpolar['NAXIS2'])

    if (a_min > 0):
        ia_min = np.argmin(np.abs(rrs - a_min))

    if (a_max > 0):
        ia_max = np.argmin(np.abs(rrs - a_max))

        
    if (M.ComputeSystVelo): 
        if (M.DoErrorMap): 
            for irrs in range(len(rrs)):
                v0_vec=im_polar[irrs,:]
                w_vec=im_polarw[irrs,:]
                av_v0=np.sum(w_vec*v0_vec)/np.sum(w_vec)
                av_cosphi=np.sum(w_vec*np.cos(phis_rad))/np.sum(w_vec)
                KepAmp = np.sum( (v0_vec - av_v0) * w_vec * np.cos(phis_rad)) / np.sum( w_vec * (np.cos(phis_rad))**2 - w_vec * np.cos(phis_rad) * av_cosphi)
                vsysts[irrs]=np.sum(w_vec*(v0_vec - KepAmp*np.cos(phis_rad)))/np.sum(w_vec)
        else:
            for irrs in range(len(rrs)):
                v0_vec=im_polar[irrs,:]  
                KepAmp = np.sum( (v0_vec - np.average(v0_vec)) * np.cos(phis_rad)) / np.sum( (np.cos(phis_rad))**2 - np.cos(phis_rad)* np.average(np.cos(phis_rad)))
                vsysts[irrs]=np.average(v0_vec - KepAmp*np.cos(phis_rad))

        vsyst=np.asscalar(np.median(vsysts[ia_min:ia_max]))
        sigma_vsyst=np.asscalar(np.std(vsysts[ia_min:ia_max]))
        print( "vsyst calculated = ",vsyst,"+-",sigma_vsyst)
        M.vsyst=vsyst

    else:
        vsyst=M.vsyst
        if (M.Verbose):
            print( "vsyst from M = ",vsyst)


    RestrictAvToRadialDomain=False # set to True is faster but may lead to discontinuities in region averages. 
    for irrs in range(len(rrs)):
        
        KepAmps[irrs] = 0.
        sKepAmps[irrs] = 0.
        AccrAmps[irrs] = 0.
        sAccrAmps[irrs] = 0.
        MeridAmps[irrs] = 0.
        sMeridAmps[irrs] = 0.
        im_polar_av[irrs,:] = 0.
        if (( (irrs < ia_min) or (irrs > ia_max) ) and RestrictAvToRadialDomain):
            continue


        
        v0_vec=im_polar[irrs,:]-vsyst
        if (M.DoErrorMap): 
            w_vec=weights[irrs,:]
        else:
            w_vec=np.ones(len(v0_vec))/typicalerror**2

        mask=(w_vec < 1E-10)
        w_vec_nozeros=w_vec
        w_vec_nozeros[mask]=1E-10
        err_vec=(1.0/np.sqrt(w_vec_nozeros))
        err_vec[mask] = 1E20 # np.inf
        w_vec[mask]=0.
                
        thisradius=rrs[irrs] #arcsec
        cosi=np.cos(M.inc)
        #Nind=2.*np.pi*thisradius * np.fabs(cosi) /(M.bmaj) # aprox number of beams at each radius
        Nind=2.*np.pi*thisradius /(M.bmaj) # aprox number of beams at each radius
        Nsum=len(w_vec)
        Npolcorr=Nsum/Nind
        #print("irrs ",irrs,Nsum, Nind, Npolcorr)
        if Npolcorr < 1:
            Npolcorr=1.
        im_Npolcorr[irrs,:]=Npolcorr
        
        if ( (np.sum(w_vec) == 0.) or ( np.sum( w_vec * (np.cos(phis_rad))**2)**2 == 0.)):
            KepAmp=0.
            sKepAmp=1E20
            AccrAmp=0.
            sAccrAmp=1E20
            MeridAmp=0.
            sMeridAmp=1E20
        else:
            if (M.InheritMumap):
                mumap_vec=mumap_polarpos[irrs,:]
                if (DoAccr):
                    sys.exit("Not yet programmed DoAccr with mumap")
                else:
                    KepAmp = np.sum(w_vec* v0_vec * mumap_vec * np.cos(phis_rad)) / np.sum(w_vec * mumap_vec *(np.cos(phis_rad))**2)
                    AccrAmp = 0. 
            else:
                Cramer = True
                if (DoMerid and Cramer):
                    sinphis = np.sin(phis_rad)
                    cosphis = np.cos(phis_rad)
                    a_1 = np.sum(w_vec * cosphis**2)
                    b_1 = np.sum(w_vec * sinphis * cosphis)
                    c_1 = np.sum(w_vec * cosphis)
                    d_1 = np.sum(w_vec * cosphis * v0_vec)
                    vard_1 = np.sum(w_vec * cosphis**2)
                    a_2 = b_1
                    b_2 = np.sum(w_vec * sinphis**2)
                    c_2 = np.sum(w_vec * sinphis)
                    d_2 = np.sum(w_vec * sinphis * v0_vec)
                    vard_2 = np.sum(w_vec * sinphis**2)
                    a_3 = c_1
                    b_3 = c_2
                    c_3 = np.sum(w_vec)
                    d_3 = np.sum(w_vec * v0_vec)
                    vard_3 = np.sum(w_vec)

                    detM = a_1*(b_2*c_3 - b_3*c_2) - b_1*(a_2*c_3 - a_3*c_2) + c_1*(a_2*b_3 - a_3*b_2)
                    if (detM == 0.):
                        #print("singular matrix")
                        continue
                        
                    A= (d_1 * (b_2*c_3 - b_3*c_2) + d_2 * (c_1*b_3 - b_1*c_3) + d_3 * (b_1*c_2-c_1*b_2))/detM
                    sigmaA= np.sqrt( vard_1*(b_2*c_3 - b_3*c_2)**2 + vard_2*(c_1*b_3 - b_1*c_3)**2 + vard_3*(b_1*c_2-c_1*b_2)**2)/detM

                    B= (d_1 * (a_3*c_2 - a_2*c_3) + d_2 *(a_1*c_3 - c_1*a_3) + d_3 * (c_1*a_2 - a_1*c_2))/detM

                    #Bcheck = (a_1 * ( d_2 *c_3 - d_3 * c_2) - d_1 * (a_2 * c_3 - a_3 * c_2) + c_1 * (a_2 * d_3 - a_3 * d_2)) /detM

                              
                    sigmaB= np.sqrt( vard_1*(a_3*c_2 - a_2*c_3)**2 + vard_2*(a_1*c_3 - c_1*a_3)**2 + vard_3*(c_1*a_2-a_1*c_2)**2)/detM

                    C= (d_1 * (a_2*b_3 - a_3*b_2) + d_2 * (b_1*a_3 - a_1*b_3) + d_3 * (a_1*b_2-b_1*a_2))/detM
                    sigmaC= np.sqrt(vard_1 * (a_2*b_3 - a_3*b_2)**2 + vard_2 * (b_1*a_3 - a_1*b_3)**2 + vard_3 * (a_1*b_2-b_1*a_2)**2)/detM

                    #if (np.fabs((Bcheck-B)/B) > 1E-8) :
                    #    print("A ",A," C" ,C)
                    #    print("B ",B,"Bcheck",Bcheck)
                    #    sys.exit("check algebra ")

                    KepAmp = A
                    sKepAmp = sigmaA

                    AccrAmp = B
                    sAccrAmp = sigmaB

                    MeridAmp = C
                    sMeridAmp = sigmaC 
                    
                    sAccrAmp *= np.sqrt(Npolcorr)
                    sKepAmp *= np.sqrt(Npolcorr)
                    sMeridAmp *= np.sqrt(Npolcorr)
                elif (DoMerid):
                    sum_w=np.sum(w_vec)
                    sum_wv0=np.sum(w_vec * v0_vec)
                    sum_wsinphi = np.sum(w_vec*np.sin(phis_rad))
                    sum_wcosphi = np.sum(w_vec*np.cos(phis_rad))
                    sinphis = np.sin(phis_rad)
                    cosphis = np.cos(phis_rad)
                    subA_num= (np.sum( w_vec * sinphis**2 - w_vec * sinphis * sum_wsinphi /sum_w ) /
                               np.sum( w_vec * sinphis * cosphis  - w_vec * cosphis * sum_wsinphi /sum_w ))

                    A_num = ( np.sum(v0_vec * w_vec * sinphis  - w_vec * sinphis * sum_wv0 / sum_w)
                              - subA_num * np.sum( v0_vec * w_vec * cosphis - w_vec * cosphis * sum_wv0 / sum_w))

                    subA_denom1_num1 = np.sum( w_vec * cosphis**2  - w_vec * cosphis * sum_wcosphi /sum_w )
                    subA_denom1_denom1 = np.sum( w_vec * cosphis * sinphis  - w_vec * cosphis * sum_wsinphi /sum_w )

                    A_denom = ( np.sum(  (w_vec * sinphis * cosphis  - w_vec * sinphis * sum_wcosphi /sum_w )
                                - (w_vec *  sinphis**2 - w_vec * sinphis * sum_wsinphi / sum_w) * subA_denom1_num1 / subA_denom1_denom1 ) )

                    A = A_num / A_denom

                    varA_num = np.sum( sinphis**2 * w_vec + w_vec**2 * sinphis**2 / sum_w) + subA_num**2 * np.sum( cosphis**2 * w_vec + w_vec**2 * cosphis**2 / sum_w)
                    sigma_A = np.sqrt(varA_num/A_denom**2)
                                       
                    B_denom = np.sum( w_vec*cosphis *sinphis - w_vec * cosphis * sum_wsinphi / sum_w)
                    subB_num = np.sum(w_vec * cosphis**2 - w_vec * cosphis * sum_wcosphi / sum_w)
                    B =  ( (np.sum( v0_vec * w_vec * cosphis - w_vec * cosphis * sum_wv0 / sum_w ) -
                            A * subB_num) / B_denom )
                    varB =  ( np.sum( w_vec * cosphis**2 + w_vec**2 * cosphis**2 / sum_w)
                                  + sigma_A**2 * subB_num**2 )/B_denom**2
                    sigma_B = np.sqrt(varB)

                    C_denom = sum_w * abs(np.cos(M.inc))
                    C = np.sum( w_vec * (v0_vec - A * cosphis - B * sinphis )) /  C_denom

                    varC = np.sum( w_vec**2 * ( err_vec**2 + sigma_A**2 * cosphis**2 + sigma_B**2 + sinphis**2 )) / C_denom**2
                    sigma_C= np.sqrt(varC)


                    KepAmp = A
                    sKepAmp = sigma_A

                    AccrAmp = B
                    sAccrAmp = sigma_B

                    MeridAmp = C
                    sMeridAmp = sigma_C 
                    
                    sAccrAmp *= np.sqrt(Npolcorr)
                    sKepAmp *= np.sqrt(Npolcorr)
                    sMeridAmp *= np.sqrt(Npolcorr)
                    
                elif (DoAccr):
                    subsum1=np.sum(w_vec * v0_vec  *  np.cos(phis_rad)) / np.sum( w_vec * (np.cos(phis_rad))**2)
                    varsubsum1=np.sum((np.sqrt(w_vec)  *  np.cos(phis_rad)) / np.sum( w_vec * (np.cos(phis_rad))**2)**2)
                    numerator=np.sum(w_vec * np.sin(phis_rad) * (v0_vec - np.cos(phis_rad) * subsum1))
                    #dum1=w_vec**2 * np.sin(phis_rad)**2
                    ##if ((~np.isfinite(varsubsum1)).any()):
                    ##    print("varsubsum1  bad values")                       
                    #dum2=np.cos(phis_rad)**2 * varsubsum1
                    #dum3=((err_vec)**2 + dum2)
                    #varnumerator=np.sum(dum1*dum3)
                    varnumerator=np.sum(w_vec**2 * np.sin(phis_rad)**2 * ((err_vec)**2 + np.cos(phis_rad)**2 * varsubsum1))
                    subsum2=np.sum(w_vec * np.sin(phis_rad) * np.cos(phis_rad)) / np.sum( w_vec * (np.cos(phis_rad))**2)
                    denom=np.sum(w_vec * ( (np.sin(phis_rad))**2  - subsum2))
                    AccrAmp=numerator/denom
                    sAccrAmp=np.sqrt(varnumerator)/denom
                    KepAmp = np.sum(w_vec * (v0_vec - AccrAmp*np.sin(phis_rad)) * np.cos(phis_rad)) / np.sum(w_vec * (np.cos(phis_rad))**2)

                    sKdenom=np.sum(w_vec * (np.cos(phis_rad))**2)
                    varKnum=np.sum( w_vec**2 * (err_vec**2 + sAccrAmp**2 * np.sin(phis_rad)**2) * np.cos(phis_rad)**2)

                    sKnum=np.sqrt(varKnum)

                    sKepAmp = sKnum/sKdenom

                    #print("sKepAmp>",rrs[irrs],KepAmp, varKnum, sKnum,sKdenom, sKepAmp, sAccrAmp)
                    
                    sAccrAmp *= np.sqrt(Npolcorr)
                    sKepAmp *= np.sqrt(Npolcorr)
                    MeridAmp = 0.
                    sMeridAmp = 1E20
                    
                else:
                    denom=np.sum(w_vec * (np.cos(phis_rad))**2)
                    numerator = np.sum(w_vec* v0_vec * np.cos(phis_rad)) 
                    KepAmp = numerator / denom
                    

                    sdenom=np.sqrt(np.sum(w_vec * (np.cos(phis_rad))**2))
                    sKepAmp = 1./sdenom

                    #print( "bmaj = ",M.bmaj, " thisradius  ",thisradius, " inc  ", M.inc ,"Nind ",Nind," Ncorr ",Ncorr, " \n")
                    sKepAmp *= np.sqrt(Npolcorr)
                    
                    AccrAmp = 0. 
                    sAccrAmp = 1E20
                    MeridAmp = 0.
                    sMeridAmp = 1E20
                    
        # else:
        # v0_vec=im_polar[irrs,:]-vsyst
        # KepAmp = np.sum(v0_vec * np.cos(phis_rad)) / np.sum((np.cos(phis_rad))**2)


        if (M.Verbose and (KepAmp < 0.) and (irrs == int(len(rrs)/4.))):
            print( "KepAmp negative, wrong PA: KepAmp=",KepAmp," PA: ",M.PA," inc",M.inc*np.pi/180.,"deg  tanpsi",tanpsi )
            print((np.array(zip(v0_vec,np.cos(phis_rad)))))

        KepAmps[irrs] = KepAmp
        sKepAmps[irrs] = sKepAmp
        AccrAmps[irrs] = AccrAmp
        sAccrAmps[irrs] = sAccrAmp
        MeridAmps[irrs] = MeridAmp
        sMeridAmps[irrs] = sMeridAmp

        if (DoMerid):
            # v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(phis_rad) * MeridAmp*np.cos(M.inc)
            v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(phis_rad) +  MeridAmp
            im_polar_av[irrs,:] = v0_vec_av + vsyst           
        elif (DoAccr):
            v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(phis_rad)
            im_polar_av[irrs,:] = v0_vec_av + vsyst
        else:
            v0_vec_av = KepAmp * np.cos(phis_rad)
            im_polar_av[irrs,:] = v0_vec_av + vsyst


    # SIGNS CALIBRATED ON THE RT TRIALS WIGGLERT
    # /strelka_ssd/simon/wiggleRT/

    v_Phi_prof = KepAmps / np.sin(M.inc)
    sv_Phi_prof = sKepAmps / np.sin(M.inc)
    v_Phi_prof = np.nan_to_num(v_Phi_prof)
    sv_Phi_prof = np.nan_to_num(sv_Phi_prof)

    #v_Phi_prof[np.isnan(v_Phi_prof)]=0.
    #sv_Phi_prof[np.isnan(sv_Phi_prof)]=0.
    
    v_R_prof = AccrAmps / np.sin(M.inc)
    sv_R_prof = sAccrAmps / np.sin(M.inc)
    v_R_prof = np.nan_to_num(v_R_prof)
    sv_R_prof = np.nan_to_num(sv_R_prof) 
    
    v_z_prof =  - MeridAmps / np.cos(M.inc)
    sv_z_prof = sMeridAmps / np.cos(M.inc)
    v_z_prof = np.nan_to_num(v_z_prof)
    sv_z_prof = np.nan_to_num(sv_z_prof) 

    
    #bmaj=hdr2['BMAJ']
    #print( "bmaj = ",bmaj,"\n";)
    #Nind=2.*np.pi*rrs /(cosi * bmaj*3600.)
    #dispv_Phi_prof=sv_Phi_prof.copy()
    #sv_Phi_prof=sv_Phi_prof/np.sqrt(Nind)



        
    ######################################################################
    # now compute chi2 in polar coords

    if (DoDCone):
        
        imazim_far=conicpolartocart(im_polar_av,inc,-tanpsi)
        im_polar_av_far_near=carttoconicpolar(imazim_far,inc,tanpsi)

        for irrs in range(len(rrs)):
            v0_vec=im_polar[irrs,:]-vsyst
            v0_vec_m=im_polar_av[irrs,:]-vsyst
            v0_vec_FN_m=im_polar_av_far_near[irrs,:]-vsyst
            if (M.DoErrorMap): 
                w_vec=weights[irrs,:]
            else:
                w_vec=1./typicalerror**2

            Delta_o=v0_vec -  v0_vec_FN_m
            Delta_m=v0_vec_m - v0_vec_FN_m
            cos_vec = np.fabs(np.cos(phis_rad))

            A =  np.sum(w_vec * cos_vec * Delta_m * (Delta_m - Delta_o)) / np.sum(w_vec * cos_vec**2 * Delta_m**2)

            if (A < 0.):
                A=0.
            # sys.exit("BUG")

            if (A > 0.5):
                #print( "setting A",A,"to 1")
                A=0.5
            
            mu_vec=(1. - A * cos_vec) 

            mumap_polarpos[irrs,:]=mu_vec
            
            # testing for chi2 improvement:
            TestChi2Improv=False
            if TestChi2Improv:
                subchi2 = np.sum(w_vec * (v0_vec   - v0_vec_m)**2)
                subchi2_mumap = np.sum(w_vec * (v0_vec   - (v0_vec_m * mu_vec + v0_vec_FN_m * (1. -mu_vec))   )**2)
                flag=''
                if (subchi2_mumap > subchi2):
                    flag='<<'
                print(("irrs %d irrs r %f A %f chi2 %f chi2_mu %f %s " % (irrs,rrs[irrs],A,subchi2,subchi2_mumap,flag)))

        
        #(nr,nphis)=im_polar_av.shape
        #im_phi=np.tile(phis_rad,(nphis,1))
        #im_cosphi=np.tile(np.cos(phis_rad),(nphis,1))
        

        #zeimage = weights*im_Npolcorr*(im_polar-im_polar_av)**2
        zeimage = weights*(im_polar-im_polar_av)**2/im_Npolcorr
        zeimage = np.nan_to_num(zeimage)
        deltaChi2 =  np.sum(zeimage,axis=1)

        deltaimage=im_polar-im_polar_av
        velodev_med=np.sqrt( np.median(deltaimage[ia_min:ia_max]**2))

        velodev_std=np.std(deltaimage[ia_min:ia_max,:])
        velodev_std_vec=np.std(deltaimage,axis=1)
        velodev_std2=np.std(velodev_std_vec[ia_min:ia_max])


        im_polar_av_DCone = im_polar_av * mumap_polarpos + im_polar_av_far_near * (1. - mumap_polarpos)

        #M.DConeDeltaChi2=False
        zeimage_DCone=weights*(im_polar - im_polar_av_DCone)**2
        deltaChi2_DCone=np.sum(zeimage_DCone,axis=1)

        chi2=sum(deltaChi2[ia_min:ia_max])
        #chi2_DCone=sum(deltaChi2_DCone[ia_min:ia_max])
        if (((np.fabs(chi2 - M.chi2_prev)/chi2) < 1E-8) and not M.DConeDeltaChi2):
            print( "using DConeDeltaChi2 chi2: ",chi2," chi2_prev", M.chi2_prev)
            M.DConeDeltaChi2=True

        if M.DConeDeltaChi2:
            deltaChi2 = deltaChi2_DCone
            
            #zeimage_N=weights*mumap_polarpos*(im_polar-im_polar_av)**2
            #zeimage_N=sp.nan_to_num(zeimage_N)
            #deltaChi2_N =  np.sum(zeimage_N,axis=1)
            #zeimage_F=weights*(1.-mumap_polarpos)*(im_polar-im_polar_av_far_near)**2
            #zeimage_F=sp.nan_to_num(zeimage_F)
            #deltaChi2_F = np.sum(zeimage_F,axis=1)
            #deltaChi2= deltaChi2_N + deltaChi2_F
            
        dispv_Phi_prof = np.sqrt(deltaChi2/np.sum(weights,axis=1))
        dispv_Phi_prof=np.nan_to_num(dispv_Phi_prof)

    else:
        #zeimage = weights*im_Npolcorr*(im_polar-im_polar_av)**2
        zeimage = weights*(im_polar-im_polar_av)**2/im_Npolcorr
 
        #zeimage=weights*(im_polar-im_polar_av)**2
        zeimage=np.nan_to_num(zeimage)
        #deltaChi2 =  np.sum(zeimage,axis=1)

        if (np.any(weights < 0.)):
            sys.exit("negative weights!!!!")
        if (np.any(zeimage < 0.)):
            sys.exit("negative chi2!!!!")

        deltaimage=im_polar-im_polar_av
        velodev_med=np.sqrt( np.median(deltaimage[ia_min:ia_max,:]**2))
        velodev_std=np.std(deltaimage[ia_min:ia_max,:])
        velodev_std_vec=np.std(deltaimage,axis=1)
        velodev_std2=np.std(velodev_std_vec[ia_min:ia_max])



        varim = deltaimage**2 * weights
        varvec = np.sum(varim,axis=1)
        wvec = np.sum(weights,axis=1)
        mask=(wvec < 1E-10)
        vec_w_var=(varvec/wvec)
        vec_w_var[mask]=0.
        vec_median_w=np.median(weights,axis=1)
        vec_typicalerror = np.sqrt(1./vec_median_w)
        deltaChi2 = (vec_w_var/vec_typicalerror**2.)
        
        #colapsed_weights=np.sum(weights,axis=1)
        #dispv_Phi_prof=colapsed_weights.copy()
        #mask=(colapsed_weights > 1E-10)
        #dispv_Phi_prof[mask]=np.sqrt(deltaChi2[mask]/colapsed_weights[mask])  # << weighted dispersion of residuals
        #dispv_Phi_prof[np.invert(mask)]=np.nan
        ## dispv_Phi_prof = np.sqrt(deltaChi2/colapsedweights)

            
        
        #sv_Phi_prof=dispv_Phi_prof.copy()
        #cosi=np.cos(M.inc)
        #Nind=2.*np.pi*rrs * cosi /(M.bmaj) # number of beams at each radius
        #sv_Phi_prof=sv_Phi_prof/np.sqrt(Nind)
        #sv_Phi_prof=np.nan_to_num(sv_Phi_prof)


    if (M.DoMerid):
        M.RadialProfile=[rrs,v_Phi_prof, sv_Phi_prof,  v_R_prof, sv_R_prof, v_z_prof, sv_z_prof]
    elif (M.DoAccr):
        M.RadialProfile=[rrs,v_Phi_prof, sv_Phi_prof,  v_R_prof, sv_R_prof]
    else:
        M.RadialProfile=[rrs,v_Phi_prof,sv_Phi_prof]

    if (DumpAllFitsFiles):
        if (M.DoMerid):
            save_prof = np.zeros((hdrpolar['NAXIS2'],7))
            save_prof[:,0] = rrs
            save_prof[:,1] = v_Phi_prof
            save_prof[:,2] = sv_Phi_prof
            save_prof[:,3] = v_R_prof 
            save_prof[:,4] = sv_R_prof
            save_prof[:,5] = v_z_prof 
            save_prof[:,6] = sv_z_prof
        elif DoAccr:
            save_prof = np.zeros((hdrpolar['NAXIS2'],5))
            save_prof[:,0] = rrs
            save_prof[:,1] = v_Phi_prof
            save_prof[:,2] = sv_Phi_prof
            save_prof[:,3] = v_R_prof 
            save_prof[:,4] = sv_R_prof
        else:
            save_prof = np.zeros((hdrpolar['NAXIS2'],3))
            save_prof[:,0] = rrs
            save_prof[:,1] = v_Phi_prof
            save_prof[:,2] = sv_Phi_prof


            
        fileout_radialprofile=re.sub('fullim.fits', 'radial_profile.dat', filename_fullim)
        np.savetxt(fileout_radialprofile, save_prof)   # x,y,z equal sized 1D arrays


        fileout_polar_Npolcorr=re.sub('fullim.fits', 'Npolcorr.fits', filename_fullim)
        hdupolar.data =im_Npolcorr
        hdupolar.writeto(fileout_polar_Npolcorr, overwrite=True)


        fileout_polar_av=re.sub('fullim.fits', 'polar_av.fits', filename_fullim)
        hdupolar.data = im_polar_av
        hdupolar.writeto(fileout_polar_av, overwrite=True)
            
        im_polar_diff = im_polar - im_polar_av
        fileout_polar_diff=re.sub('fullim.fits', 'polar_diff.fits', filename_fullim)
        hdupolar.data = im_polar_diff
        hdupolar.writeto(fileout_polar_diff, overwrite=True)

        if (DoDCone):
            fileout_polar_av_far_near=re.sub('fullim.fits', 'polar_av_far_near.fits', filename_fullim)
            hdupolar.data = im_polar_av_far_near
            hdupolar.writeto(fileout_polar_av_far_near, overwrite=True)

            fileout_polar_av_far_near=re.sub('fullim.fits', 'polar_av_DCone.fits', filename_fullim)
            hdupolar.data = im_polar_av_DCone
            hdupolar.writeto(fileout_polar_av_far_near, overwrite=True)




    ######################################################################
    # Polar average: back to sky  for diff image

    
    if (M.ComputeSkyImages):
    
        if (M.Verbose):
            print( "CONICPOLAR2CARTESIAN TRANSFORM FOR AZIM AV START im_polar_av")
            print( "using inc ",inc*np.pi/180.,"deg tanpsi ",tanpsi)

        start_time=time.time()

        imazim=conicpolartocart(im_polar_av,inc,tanpsi)

        faceoninc=np.pi

        # if (tanpsi < 0.):
        #    faceoninc=np.pi 

        if (inc > np.pi/2.):
            faceoninc=0.

        imazim_faceon=polartocart(im_polar_av,faceoninc)
        resamp_faceon=polartocart(im_polar,faceoninc)

        im4=imazim #gridding(fileout_stretched_av,hdr3)

        if (M.Verbose):
            print( "CONICPOLAR2CARTESIAN TRANSFORM FOR AZIM AV  DONE in",time.time()-start_time)

        if (DumpAllFitsFiles):
            fileout_azimav=re.sub('fullim.fits', 'azim_av.fits', filename_fullim)
            pf.writeto(fileout_azimav,im4, hdr2, overwrite=True)



        # back to sky - rotate
        im4drot = ndimage.rotate(im4, -rotangle, reshape=False,order=0)

        # diff
        diff_im = resamp-im4drot
        diff_im_faceon = resamp_faceon-imazim_faceon





        im_polar_av_region=im_polar_av.copy()
        im_polar_av_region[0:ia_min]=0.
        im_polar_av_region[ia_min:ia_max]=1.
        im_polar_av_region[ia_max:]=0.
        imazim_region = conicpolartocart(im_polar_av_region,inc,tanpsi)
        imazim_region_drot = ndimage.rotate(imazim_region, -rotangle, reshape=False,order=0)
        imazim_region_faceon = polartocart(im_polar_av_region,faceoninc)

        mask=np.where(imazim_region_drot > 0.9)
        skychi2= np.sum(weights[mask]*diff_im[mask]**2)/M.Ncorr
        M.skychi2=skychi2


        hdudiff=pf.PrimaryHDU()
        hdudiff.data=diff_im
        hdudiff.header=hdr2
        M.Hdudiff=hdudiff

        hdudiff_faceon=pf.PrimaryHDU()
        hdudiff_faceon.data=diff_im_faceon
        hdudiff_faceon.header=hdr2
        M.Hdudiff_faceon=hdudiff_faceon

        hduresamp_faceon=pf.PrimaryHDU()
        hduresamp_faceon.data=resamp_faceon
        hduresamp_faceon.header=hdr2
        M.Hduresamp_faceon=hduresamp_faceon

        hduregion=pf.PrimaryHDU()
        hduregion.data=imazim_region_drot
        hduregion.header=hdr2
        M.Hduregion=hduregion

        hdumoddrot=pf.PrimaryHDU()
        hdumoddrot.data=im4drot
        hdumoddrot.header=hdr2
        M.Hdumoddrot=hdumoddrot

        hduregion_faceon=pf.PrimaryHDU()
        hduregion_faceon.data=imazim_region_faceon
        hduregion_faceon.header=hdr2
        M.Hduregion_faceon=hduregion_faceon
        
        if (DumpAllFitsFiles):
            fileout_drotated=re.sub('fullim.fits', 'azim_av_drot.fits', filename_fullim)
            pf.writeto(fileout_drotated,im4drot, hdr2, overwrite=True)
            fileout_drotated=re.sub('fullim.fits', 'azim_av_drot_diff.fits', filename_fullim)
            pf.writeto(fileout_drotated,diff_im, hdr2, overwrite=True)
            fileout_drotated_region=re.sub('fullim.fits', 'region_drot.fits', filename_fullim)
            pf.writeto(fileout_drotated_region,imazim_region_drot, hdr2, overwrite=True)

            fileout_diff_im_faceon=re.sub('fullim.fits', 'diff_faceon.fits', filename_fullim)
            pf.writeto(fileout_diff_im_faceon,diff_im_faceon, hdr2, overwrite=True)
            
            fileout_im_faceon=re.sub('fullim.fits', 'resamp_faceon.fits', filename_fullim)
            pf.writeto(fileout_im_faceon,resamp_faceon, hdr2, overwrite=True)

            fileout_region_faceon=re.sub('fullim.fits', 'region_faceon.fits', filename_fullim)
            pf.writeto(fileout_region_faceon,imazim_region_faceon, hdr2, overwrite=True)



        ##fileout_proj=re.sub('fullim.fits', 'azim_av_proj.fits', filename_fullim)
        ##if (not OptimOrient):
        ##    pf.writeto(fileout_proj,im4, hdr2, overwrite=True)



        if (DoDCone):


            ######################################################################
            # average +-psi in sky plane
            
            imazim_far=conicpolartocart(im_polar_av,inc,-tanpsi)
            imazim_far_drot = ndimage.rotate(imazim_far, -rotangle, reshape=False,order=0)
            
            mumap0=conicpolartocart(mumap_polarpos,inc,tanpsi)
            mumap=ndimage.rotate(mumap0, -rotangle, reshape=False)

            if M.DConeDeltaChi2:
                M.mumap=mumap

            immod_DCone  = mumap * im4drot  + (1. - mumap) * imazim_far_drot
            diff_im_DCone = resamp - immod_DCone

            # immod_DCone_b_rot = conicpolartocart(im_polar_av_DCone,inc,tanpsi)
            # immod_DCone_b=ndimage.rotate(immod_DCone_b_rot, -rotangle, reshape=False)
            # diff_im_DCone_b = resamp - immod_DCone_b

            if DoDCone:
                hduDConemoddrot=pf.PrimaryHDU()
                hduDConemoddrot.data=immod_DCone
                hduDConemoddrot.header=hdr2
                M.HduDConemoddrot=hduDConemoddrot

                hdudiffDConemoddrot=pf.PrimaryHDU()
                hdudiffDConemoddrot.data=diff_im_DCone
                hdudiffDConemoddrot.header=hdr2
                M.HdudiffDConemoddrot=hdudiffDConemoddrot

                hdumumap=pf.PrimaryHDU()
                hdumumap.data=mumap
                hdumumap.header=hdr2
                M.Hdumumap=hdumumap

            
            if (DumpAllFitsFiles):
                fileout_mumap=re.sub('fullim.fits', 'mumap.fits', filename_fullim)
                pf.writeto(fileout_mumap,mumap, hdr2, overwrite=True)

                fileout_azimav_far=re.sub('fullim.fits', 'azim_av_far.fits', filename_fullim)
                pf.writeto(fileout_azimav_far,imazim_far, hdr2, overwrite=True)

                fileout_azimav_far_drot=re.sub('fullim.fits', 'azim_av_far_drot.fits', filename_fullim)
                pf.writeto(fileout_azimav_far_drot,imazim_far_drot, hdr2, overwrite=True)

                fileout_immod_DCone=re.sub('fullim.fits', 'immod_DCone.fits', filename_fullim)
                pf.writeto(fileout_immod_DCone,immod_DCone, hdr2, overwrite=True)
                
                # for  some reason _b does not work, there seems to be a background 
                # fileout_immod_DCone_b=re.sub('fullim.fits', 'immod_DCone_b.fits', filename_fullim)
                # pf.writeto(fileout_immod_DCone_b,immod_DCone_b, hdr2, overwrite=True)

                # for  some reason _b does not work, there seems to be a background 
                # fileout_immod_DCone_a_b=re.sub('fullim.fits', 'immod_DCone_a-b.fits', filename_fullim)
                # pf.writeto(fileout_immod_DCone_a_b,immod_DCone-immod_DCone_b, hdr2, overwrite=True)

                fileout_immod_DCone_diff=re.sub('fullim.fits', 'diff_DCone.fits', filename_fullim)
                pf.writeto(fileout_immod_DCone_diff,diff_im_DCone, hdr2, overwrite=True)

                # for  some reason _b does not work, there seems to be a background 
                # fileout_immod_DCone_diff_b=re.sub('fullim.fits', 'diff_DCone_b.fits', filename_fullim)
                # pf.writeto(fileout_immod_DCone_diff_b,diff_im_DCone_b, hdr2, overwrite=True)

                
    
    
    ######################################################################
    # CROSS CHECK INVERSE TRANSFORM

    if (M.XCheckInv):
        
        
        print( "CONICPOLAR2CARTESIAN TRANSFORM FOR XCHECK INV  START XCheckInv")
        start_time=time.time()
        im_x=conicpolartocart(im_polar,inc,tanpsi)

        
        print( "CONICPOLAR2CARTESIAN TRANSFORM FOR XCHECK INV   DONE in",time.time()-start_time)


        fileout_stretched_x=re.sub('fullim.fits', 'stretched_x.fits', filename_fullim)
        pf.writeto(fileout_stretched_x,im_x, hdr2, overwrite=True)
        
        #hdr3 = deepcopy(hdr2)
        fileout_proj_x=re.sub('fullim.fits', 'x_proj.fits', filename_fullim)
        #im4_x=gridding(fileout_stretched_x,hdr3)
        pf.writeto(fileout_proj_x,im_x, hdr2, overwrite=True)
        
        im4_x_drot = ndimage.rotate(im_x, -rotangle, reshape=False)

        fileout_drotated_x=re.sub('fullim.fits', 'x_drot.fits', filename_fullim)
        pf.writeto(fileout_drotated_x,im4_x_drot, hdr2, overwrite=True)
    

        fileout_skydiff=re.sub('fullim.fits', 'x_diff.fits', filename_fullim)
        pf.writeto(fileout_skydiff,resamp-im4_x_drot, hdr2, overwrite=True)
    

    
    
    ######################################################################
    # profile plotting

    if (PlotRadialProfile and DumpAllFitsFiles):
        
        
        # -----------------------------------------------------------
        # nice fonts
        # -----------------------------------------------------------
        matplotlib.rc('font', family='sans-serif') 
        matplotlib.rcParams.update({'font.size': 9})
        
        
        plt.figure(figsize=(10, 8))
        axprofile=plt.subplot(111)
        
        rmax=np.max(rrs)
        
        plt.setp(axprofile.get_xticklabels(),visible=True) #, fontsize=6)
        plt.setp(axprofile.get_yticklabels(),visible=True) #, fontsize=6)
        plt.xlim(0.,rmax)
        #plt.ylim(np.min(v_Phi_prof),1.1*np.max(v_Phi_prof))

        plt.ylim(-0.1*np.max(v_Phi_prof[ia_min:ia_max]),1.1*np.max(v_Phi_prof[ia_min:ia_max]))
        plt.plot(rrs,v_Phi_prof,color='grey',linewidth=0.1,linestyle='solid')
        plt.fill_between(rrs, v_Phi_prof+sv_Phi_prof, v_Phi_prof-sv_Phi_prof, lw=0.1,color='r', alpha=0.3, interpolate=True, step='mid')

        if (DoMerid):
            plt.plot(rrs,v_z_prof,color='orange',linewidth=1.,linestyle='solid',alpha=0.5,label='v_z')
            plt.plot(rrs,v_R_prof,color='green',linewidth=1.,linestyle='solid',alpha=0.5,label='v_R')
        elif (DoAccr):
            plt.plot(rrs,v_R_prof,color='green',linewidth=1.,linestyle='solid',alpha=0.5,label='v_R')
            

            
        #print( np.min(v_Phi_prof))
        #print( np.max(v_Phi_prof))
        
        plt.plot(rrs,rrs*0.,color='black',linewidth=0.1,linestyle='solid')
        
        #plt.plot(rrs,Icutav,color='blue',linewidth=1,linestyle='solid')
        
        
        #plt.fill_between(rrs, v_Phi_prof+dispv_Phi_prof, v_Phi_prof-dispv_Phi_prof, lw=0.1,color='b', alpha=0.3, interpolate=True)
        
        plt.ylabel(r'$\langle |v_{\circ}(r)| \rangle$')
    
        plt.xlabel(r'$r$ / arcsec')
        
    
    
        fileout_fig=re.sub('fullim.fits', 'fig_profile.pdf', filename_fullim)
        plt.savefig(fileout_fig, bbox_inches='tight')



    ######################################################################

    chi2=sum(deltaChi2[ia_min:ia_max])
    M.chi2_prev=chi2
    M.velodev_med=velodev_med
    M.velodev_std=velodev_std
    M.velodev_std2=velodev_std2
    
    #print( "chi2: ",chi2," typical error ",typicalerror," velodev_std ",velodev_std," velodev_std2 ",velodev_std2,"velodev_med",velodev_med)

    # print( "chi2: ",chi2," typical error ",typicalerror," velodev_med ",velodev_med)
    

    M.polarchi2=chi2
    retchi2=chi2  #  / M.Ncorr ## very aprox correction for correlated pixels because retchi2 is in the polar domain, so NCorr is variable over the polar map. 
    return retchi2


    
