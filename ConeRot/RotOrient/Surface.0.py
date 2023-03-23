import sys
import numpy as np
import scipy  as sp
import os
import os.path
from scipy import ndimage
from astropy.io import fits 
from copy import deepcopy

from scipy.optimize import curve_fit



from multiprocessing import Pool
from tqdm import tqdm


if not sys.warnoptions:
    import os, warnings
    #warnings.simplefilter("default") # Change the filter in this process
    warnings.simplefilter("ignore") # Change the filter in this process
    #os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses
    os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

from ImUtils.Resamp import gridding




#import MPolarMaps.Master

import PyVtools.Vtools as Vtools



def ftaper_gap(r,r1):
    taper=(0.1+ ((np.tanh((r-r1)/r1))**5)*(1.-0.1))
    #taper=r/r
    return taper

def ftaper(r,r1):
    taper=(1.-np.tanh((r-r1)/r1))**3.
    return taper

def z_func_gap(r,z0,r0,q,r1,r2):
    taper_trunc=1.
    powerdisk=1.
    if (isinstance(r,np.ndarray)):
        taper_trunc=np.ones(r.shape)
        alltaper=ftaper(r,r2)
        taper_trunc[(r >r2)] = alltaper[(r>r2)]
    else:
        if (r>r2):
            taper_trunc= ftaper(r,r2)
    if (isinstance(r,np.ndarray)):
        powerdisk=np.ones(r.shape)
        powerdisk[(r <r0)] = (r[(r<r0)]/r0)**q
    else:
        if (r<r0):
            powerdisk=(r/r0)**q


            
    if (isinstance(r,np.ndarray)):
        taper_gap=np.ones(r.shape)
        alltaper=ftaper_gap(r,r1)
        taper_gap[(r < r1)] = 0.1
        taper_gap[(r >= r1)] = alltaper[(r>=r1)]
    else:
        if (r<r1):
           taper_gap= 0.1
        else:
           taper_gap= ftaper_gap(r,r1)
            
    return z0 * powerdisk*taper_gap*taper_trunc


def z_func(r,z0,r0,r1,q):
    taper=1.
    if (isinstance(r,np.ndarray)):
        taper=np.ones(r.shape)
        alltaper=ftaper(r,r1)
        taper[(r >r1)] = alltaper[(r>r1)]
    else:
        if (r>r1):
            taper= ftaper(r,r1)
    return (z0 * (r/r0)**q)*taper


#def ftanpsi(r,z0,r0,r1,q):
#    z=z_func(r,z0,r0,r1,q)
#    return z/r
#
#def ftanpsi_gap(r,z0,r0,q,r1,r2):
#    z=z_func_gap(r,z0,r0,q,r1,r2)
#    return z/r

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
    rho_m=(-b-np.sqrt(Delta))/(2.*a)
    rho_p=(-b+np.sqrt(Delta))/(2.*a)
    if (rho_p > 0.):
        print("rho_p > 0",rho_p, rho_m,x,y,inc*180./np.pi,tanpsi,np.arctan(tanpsi)*180./np.pi)

    rho=rho_m
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

    H1=tanpsi0*rho
    num= x - H1 * np.sin(inc)
    denom= rho * ( (1./np.cos(inc))  - np.tan(inc) * np.sin(inc))
    sintheta= num/ denom
    
    theta=np.arccos(costheta)
    
    if sintheta<0:
        theta = 2.*np.pi - theta


    
    thetaindex = (theta * (float(nx)-1.) / (2. * np.pi)) 


            
    return (rindex,thetaindex)


def cartesian2offsetpolar(outcoords, inputshape, origin, inc=0., Hpix=0.):

    rindex, phiindex = outcoords
    x0, y0 = origin

    side=float(inputshape[0])
     
    phi = float(phiindex) * 2. * np.pi / (side-1.)
    phi = 2.*np.pi - phi
    
    y = rindex*np.cos(phi)

    x = rindex*np.sin(phi)*np.cos(inc)  + Hpix *np.sin(inc) 
        
    ix = -x + x0
    iy = y + float(y0)
    
    return (iy,ix)


def offsetpolar2cartesian(outcoords, inputshape, origin,inc=0.,Hpix=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx=inputshape[0]-1
    ny=inputshape[1]-1
    x = -float(xindex - x0)
    y = float(yindex - y0)

    side=float(inputshape[0])
    
    rho=np.sqrt(((x-Hpix*np.sin(inc))**2/(np.cos(inc)**2))+y**2)
    rindex = rho
    if (rho == 0.):
        cosphi = 0.
    else:
        cosphi = y / rho

    # theta=np.arccos(costheta)
    phi = np.arctan2((-(x-Hpix*np.sin(inc))/np.cos(inc)), y) 
    if (phi < 0):
        phi = phi + 2.*np.pi

    phi = 2.*np.pi - phi
    phiindex = (phi * float(nx) / (2. * np.pi))
            
    return (rindex,phiindex)


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


def offsetpolartocart(im_polar,inc,Hpix):

    (ny,nx)=im_polar.shape
    (i0,j0)=(((float(nx)+1.)/2.)-1.,((float(ny)+1.)/2.)-1.)

    
    im_cart = sp.ndimage.geometric_transform(im_polar,offsetpolar2cartesian,
                                            order=1,
                                            output_shape = (im_polar.shape[0], im_polar.shape[1]),
                                            extra_keywords = {'inputshape':im_polar.shape,
                                                  'inc':inc,'Hpix':Hpix,
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




def proc_1region(regionparams):
    
    Rmesh1=regionparams['Rmesh1']
    Rmesh2=regionparams['Rmesh2']
    inc=regionparams['inc']
    tanpsi=regionparams['tanpsi']
    rrs_polar=regionparams['rrs_polar']
    HHs_polar=regionparams['HHs_polar']
    phis_polar=regionparams['phis_polar']
    pixscale=regionparams['pixscale']

    Hpix=tanpsi*Rmesh1/pixscale

    domain_polar=np.zeros(rrs_polar.shape)
    domain_polar[( rrs_polar < Rmesh2)]=1.

    region_polar=np.zeros(rrs_polar.shape)
    region_polar[( rrs_polar >= Rmesh1) & ( rrs_polar < Rmesh2)]=1.
    
    if Debug:
        print("calling  conicpolartocart with ", inc,tanpsi)


    if ConicPolar:
    
        phis_sky_domain_top=conicpolartocart(phis_polar,inc,tanpsi)
        phis_sky_domain_bottom=conicpolartocart(phis_polar,inc,-tanpsi)
        rrs_sky_domain_top=conicpolartocart(rrs_polar,inc,tanpsi)
        rrs_sky_domain_bottom=conicpolartocart(rrs_polar,inc,-tanpsi)
        HHs_sky_domain_top=conicpolartocart(HHs_polar,inc,tanpsi)
        HHs_sky_domain_bottom=-conicpolartocart(HHs_polar,inc,-tanpsi)
        
    elif OffsetPolar:

        phis_sky_domain_top=offsetpolartocart(phis_polar,inc,Hpix)
        phis_sky_domain_bottom=offsetpolartocart(phis_polar,inc,-Hpix)
        rrs_sky_domain_top=offsetpolartocart(rrs_polar,inc,Hpix)
        rrs_sky_domain_bottom=offsetpolartocart(rrs_polar,inc,-Hpix)
        HHs_sky_domain_top=offsetpolartocart(HHs_polar,inc,Hpix)
        HHs_sky_domain_bottom=-offsetpolartocart(HHs_polar,inc,-Hpix)

    else:
        sys.exit("choose transform")
        
    passout={
        'HHs_sky_domain_top':HHs_sky_domain_top,
        'HHs_sky_domain_bottom':HHs_sky_domain_bottom,
        'rrs_sky_domain_top':rrs_sky_domain_top,
        'rrs_sky_domain_bottom':rrs_sky_domain_bottom,
        'phis_sky_domain_top':phis_sky_domain_top,
        'phis_sky_domain_bottom':phis_sky_domain_bottom}
             
    return passout


def punch_skymap(im_sky,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='H_top_sky.fits'):
    
    rotangle= -PA
    im_sky_rot = ndimage.rotate(im_sky, rotangle, reshape=False)
    hdu = fits.PrimaryHDU()
    hdu.header=hdr_canvas
    hdu.data=im_sky_rot
    hdu0=gridding(hdu,hdr_canvas0,ReturnHDUList=True)
    fileout=fileouttag+fileout_basename
    hdu0.writeto(fileout, overwrite=True)
    return hdu0


def gen_surface(file_canvas,file_psiprofile=False,PA=0.,inc=0.,fileouttag='H',ForceTop=False,ncores=4,Verbose=False,RunTransforms=False,nrmesh=36,z0=0.4,r0=1.,q=1.,r1=2.,r2=4,zoomfactor=1.,DoConicPolar=False,DoOffsetPolar=True):

    global Debug
    Debug=Verbose
    global ConicPolar
    global OffsetPolar
    ConicPolar=DoConicPolar
    OffsetPolar=DoOffsetPolar
    if ConicPolar:
        OffsetPolar=False

    if file_psiprofile:
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


    hdr_canvas = deepcopy(hdr_canvas0)
    hdr_canvas['CDELT1']*=zoomfactor
    hdr_canvas['CDELT2']*=zoomfactor

    f1=gridding(f0,hdr_canvas,ReturnHDUList=True)
    im_canvas = f1[0].data
    


    #Vtools.View(im_canvas)
    

    pixscale=hdr_canvas['CDELT2']*3600.
    
    (ny,nx) = im_canvas.shape 
    x=np.arange(1,nx+1)
    y=np.arange(1,ny+1)
    X, Y = np.meshgrid(x, y)
    
    X0 = np.floor(nx/2)+1
    Y0 = np.floor(ny/2)+1
    dxxs = -pixscale *(X-X0)
    dyys = pixscale *(Y-Y0)
    rrs = np.sqrt( (dxxs)**2 + (dyys)**2)

    rmax=np.max(rrs)
    if Debug:
        print("rmax",rmax)

    if file_psiprofile:
        
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

        if Debug:
            print("rregions",rregions)
            print("Hs",Hs)
    

        #fH = interp1d(rregions, Hs, kind='cubic')
        #ftanpsi = interp1d(rregions, tanpsis, kind='cubic')
        ##finvH = interp1d(Hs, rregions,  kind='cubic')
    

        popt, pcov = curve_fit(z_func, rregions, tanpsis, p0=[0.,0.5,1.,1.],bounds=[[-10.,0.,0.,0.1],[10.,10.,10.,4.]])


        obsprof=np.zeros((len(rregions),2))
        modprof=np.zeros((len(rregions),2))
        obsprof[:,0]=rregions
        obsprof[:,1]=tanpsis

    
        z0=popt[0]
        r0=popt[1]
        r1=popt[2]
        q=popt[3]

        if Debug:
            print("z0 ",z0, "r0 ",r0, "r1 ",r1,"q", q)

            
        Retro=False
    
        #if z0 < 0.:
        #    Hsign=-1
        #    Retro=True
        
        modprof[:,0]=rregions
        modprof[:,1]=z_func(rregions,z0,r0,r1,q)
        if Debug:
            Vtools.Spec([obsprof,modprof])


    #HHs=fH(rrs)
    #HHs=z_func(rrs,z0,r0,r1,q)
    #r2=4.
    #r1=2.
    
    HHs=z_func_gap(rrs,z0,r0,q,r1,r2)

    master_Hsign=np.sign(z0)
    
    fileout=fileouttag+'H_faceon.fits'
    fits.writeto(fileout,HHs, hdr_canvas, overwrite=True)


    
    ########################################################################
    ## default expansion
    #
    #M.workdir='polarmaps_default/'  # directory for products 
    #M.prep_files()
    #M.polar_expansions()

    nphis=HHs.shape[0]
    if Debug:
        print("nphis",nphis)
    nrs=HHs.shape[1]
    if Debug:
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

    

    if Debug and not file_psiprofile:
        modprof=np.zeros((len(rs),2))

        modprof[:,0]=rs
        modprof[:,1]=z_func_gap(rs,z0,r0,q,r1,r2)

        #modprof[:,1]=z_func(rs,z0,r0,r1,q)
        print("model H(R)")
        Vtools.Spec([modprof])

        #modprof[:,1]=ftanpsi(rs,z0,r0,r1,q)
        #modprof[:,1]=ftanpsi_gap(rs,z0,r0,q,r1,r2)
        modprof[:,1]=z_func_gap(rs,z0,r0,q,r1,r2)/rs
        print("model h(R)")
        Vtools.Spec([modprof,])

        #modprof[:,1]=np.arctan(ftanpsi(rs,z0,r0,r1,q))*180./np.pi
        modprof[:,1]=np.arctan(z_func_gap(rs,z0,r0,q,r1,r2)/rs)*180./np.pi
        print("model psi(R)")
        Vtools.Spec([modprof,])

    rrs_polar=np.zeros(HHs_polar.shape)

    print("nphis",nphis,"nx",nx)

    x=np.arange(1,nx+1)
    y=np.arange(1,ny+1)
    X, Y = np.meshgrid(x, y)
    
    rrs_polar= 3600.*(Y-hdrpolar['CRPIX2']+1.0)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']

        
    phis_polar= (X-hdrpolar['CRPIX1']+1.0)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']
                                
    if Debug:
        Vtools.View(phis_polar)


    zero_offset=hdr_canvas['CDELT2']*3600.
    rmesh = zero_offset+ (np.arange(nrmesh)/(nrmesh-1))*rmax
    tasks=[]

    for iregion in list(range(nrmesh-1)):
        Rmesh1=rmesh[iregion]
        Rmesh2=rmesh[iregion+1]
        tanpsi=z_func_gap(Rmesh1,z0,r0,q,r1,r2)/Rmesh1
        psi_deg_mod = np.fabs(np.arctan(tanpsi)) *180. / np.pi
        inc_deg_mod = inc * 180. / np.pi
        if (inc_deg_mod > 90.):
            inc_deg_mod=180.-inc_deg
        if Debug:
            print("Rmesh1: ", Rmesh1, "Rmesh2: ",Rmesh2)
            print("tanpsi",tanpsi,"psi_deg_mod",psi_deg_mod,"inc_deg_mod",inc_deg_mod)
        if (ConicPolar & (psi_deg_mod  >  (90. - inc_deg_mod))):
            print("psi_deg_mod",psi_deg_mod,"inc_deg_mod",inc_deg_mod)
            sys.exit("opening angle too steep, no solutions via single-valued conic transforms -> develop bi-valued transforms")

        regionparams={'Rmesh1':Rmesh1,'Rmesh2':Rmesh2,'inc':inc,'tanpsi':tanpsi,'rrs_polar':rrs_polar,'HHs_polar':HHs_polar,'phis_polar':phis_polar,'pixscale':pixscale}
        tasks.append(regionparams)

    domain_polar=np.zeros(rrs_polar.shape)

    datafile=fileouttag+'binfile_Pooloutput.npy'
    if RunTransforms:
        with Pool(ncores) as pool:
            Pooloutput = list(tqdm(pool.imap(proc_1region, tasks), total=len(tasks)))
            pool.close()
            pool.join()
            np.save(datafile,Pooloutput)
    else:
        Pooloutput=np.load(datafile,allow_pickle=True)
        
            
    HHs_sky_top=np.zeros(HHs.shape)    
    HHs_sky_bottom=np.zeros(HHs.shape)
    phis_sky_top=np.zeros(HHs.shape)    
    phis_sky_bottom=np.zeros(HHs.shape)
    rrs_sky_top=np.zeros(HHs.shape)    
    rrs_sky_bottom=np.zeros(HHs.shape)

    for iregion,aregion in enumerate(Pooloutput):
        regionparams=tasks[iregion]
        Rmesh1=regionparams['Rmesh1']
        Rmesh2=regionparams['Rmesh2']
        
        HHs_sky_domain_top=aregion['HHs_sky_domain_top']
        HHs_sky_domain_bottom=aregion['HHs_sky_domain_bottom']  
        rrs_sky_domain_top=aregion['rrs_sky_domain_top']            
        rrs_sky_domain_bottom=aregion['rrs_sky_domain_bottom']         
        phis_sky_domain_top=aregion['phis_sky_domain_top']           
        phis_sky_domain_bottom=aregion['phis_sky_domain_bottom']

        
        inc_deg = inc * 180. / np.pi
        tanpsi=tasks[iregion]['tanpsi']
        psi_deg = np.arctan(tanpsi) *180. / np.pi
        
        if Debug:
            print("region tanpsi",tanpsi,"psi_deg",psi_deg,"inc_deg",inc_deg)
            

        maskN=( ((rrs_sky_domain_top >= Rmesh1) | (master_Hsign*np.sign(HHs_sky_top) <= 0. ) ) &  (rrs_sky_domain_top < Rmesh2) & (phis_sky_domain_top > np.pi) )

        HHs_sky_top[maskN]= HHs_sky_domain_top[maskN]
        phis_sky_top[maskN]=phis_sky_domain_top[maskN]
        rrs_sky_top[maskN]=rrs_sky_domain_top[maskN]

        maskF=( ((rrs_sky_domain_top >= Rmesh1) | (master_Hsign*np.sign(HHs_sky_top) <= 0. ) ) &  (rrs_sky_domain_top < Rmesh2) & (phis_sky_domain_top <= np.pi) & (master_Hsign*np.sign(HHs_sky_top) <= 0.))

        #maskF=( ((rrs_sky_domain_top >= Rmesh1) ) &  (rrs_sky_domain_top < Rmesh2) & (phis_sky_domain_top <= np.pi) & (master_Hsign*np.sign(HHs_sky_top) <= 0.))

        
        HHs_sky_top[maskF]= HHs_sky_domain_top[maskF]
        phis_sky_top[maskF]=phis_sky_domain_top[maskF]
        rrs_sky_top[maskF]=rrs_sky_domain_top[maskF]
        
        if ((Rmesh1 > r1) & Debug):
            print("r2max",r2max)
            Vtools.View(phis_sky_domain_top)
            print("Mask N")
            arrmaskN=np.zeros(rrs_sky_top.shape)
            arrmaskN[maskN]=1.
            Vtools.View(arrmaskN)
            print("Mask F")
            arrmaskF=np.zeros(rrs_sky_top.shape)
            arrmaskF[maskF]=1.
            Vtools.View(arrmaskF)
            print("combine H")
            Vtools.View(HHs_sky_top)
               

        maskN=( ((rrs_sky_domain_bottom >= Rmesh1) | (-master_Hsign*np.sign(HHs_sky_bottom) <= 0. ) ) &  (rrs_sky_domain_bottom < Rmesh2) & (phis_sky_domain_bottom > np.pi) )
        

        HHs_sky_bottom[maskN]= HHs_sky_domain_bottom[maskN]
        phis_sky_bottom[maskN]=phis_sky_domain_bottom[maskN]
        rrs_sky_bottom[maskN]=rrs_sky_domain_bottom[maskN]

        maskF=( ((rrs_sky_domain_bottom >= Rmesh1) | (-master_Hsign*np.sign(HHs_sky_bottom) <= 0. ) ) &  (rrs_sky_domain_bottom < Rmesh2) & (phis_sky_domain_bottom <= np.pi) & (-master_Hsign*np.sign(HHs_sky_bottom) <= 0.))
        
        HHs_sky_bottom[maskF]= HHs_sky_domain_bottom[maskF]
        phis_sky_bottom[maskF]=phis_sky_domain_bottom[maskF]
        rrs_sky_bottom[maskF]=rrs_sky_domain_bottom[maskF]
 

        
    hdu_H_top_sky=punch_skymap(HHs_sky_top,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='H_top_sky.fits')

    hdu_phis_top_sky=punch_skymap(phis_sky_top,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='phis_top_sky.fits')

    hdu_rrs_top_sky=punch_skymap(rrs_sky_top,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='rrs_top_sky.fits')

    hdu_H_bottom_sky=punch_skymap(HHs_sky_bottom,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='H_bottom_sky.fits')

    hdu_phis_bottom_sky=punch_skymap(phis_sky_bottom,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='phis_bottom_sky.fits')

    hdu_rrs_bottom_sky=punch_skymap(rrs_sky_bottom,hdr_canvas,hdr_canvas0,PA,fileouttag,fileout_basename='rrs_bottom_sky.fits')


    #punch_skymap(-HHs_sky_top,hdr_canvas,hdr_canvas0,180.,fileouttag,fileout_basename='H_top_xcheck_180degrot.fits')
    
    return


