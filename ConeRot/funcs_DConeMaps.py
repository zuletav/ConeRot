import sys
import numpy as np
import scipy as sp
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
from copy import deepcopy
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.colors as colors

# include_path = '/home/simon/common/python/include/'
# sys.path.append(include_path)
from ImUtils.Resamp import gridding
from ImUtils.Cube2Im import slice0
import ConeRot.TakeAzAv as TakeAzAv
import ConeRot.ConicTransforms_numba as ConicTransforms
# import PyVtools.Vtools as Vtools

if not sys.warnoptions:
    import os, warnings
    #warnings.simplefilter("default") # Change the filter in this process
    warnings.simplefilter("ignore")  # Change the filter in this process
    #os.environ["PYTHONWARNINGS"] = "default" # Also affect subprocesses
    os.environ["PYTHONWARNINGS"] = "ignore"  # Also affect subprocesses


def cartesian2conicpolar(outcoords, inputshape, origin, inc=0., tanpsi=0.):
    """Coordinate transform for converting a conic polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    rindex, thetaindex = outcoords
    x0, y0 = origin

    #theta = float(thetaindex) * 2. * np.pi / float(inputshape[0])
    theta = thetaindex * 2. * np.pi / inputshape[0]

    y = rindex * np.cos(theta)

    height = tanpsi * rindex
    #height = -tanpsi * rindex
    x = (rindex*np.sin(theta))/np.cos(inc) + \
             (height - (rindex*np.sin(theta))*np.tan(inc))*np.sin(inc)

    ix = -x + x0
    iy = y + y0

    return (iy, ix)


def cartesian2polar(outcoords, inputshape, origin, inc=0.):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    rindex, thetaindex = outcoords
    x0, y0 = origin

    #theta = float(thetaindex) * 2. * np.pi / float(inputshape[0] - 1.)
    theta = thetaindex * 2. * np.pi / (inputshape[0] - 1.)

    y = rindex * np.cos(theta)

    x = rindex * np.sin(theta) * np.cos(inc)

    ix = -x + x0
    iy = y + y0

    return (iy, ix)


def conicpolar2cartesian_ellipse_reflectedazimuths(outcoords,
                                                   inputshape,
                                                   origin,
                                                   inc=0.,
                                                   tanpsi=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx = inputshape[0]
    ny = inputshape[1]
    x = -float(xindex - x0)
    y = float(yindex - y0)

    tanpsi0 = tanpsi

    a = ((np.tan(inc) * tanpsi0)**2 - 1.0)
    b = -2. * x * np.sin(inc) * tanpsi0 / (np.cos(inc))**2
    c = y**2 + (x**2 / (np.cos(inc)**2))
    Delta = b**2 - 4. * a * c
    rho = (-b - np.sqrt(Delta)) / (2. * a)
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

    theta = np.arccos(costheta)
    thetaindex = (theta * float(nx) / (2. * np.pi))

    return (rindex, thetaindex)


def conicpolar2cartesian_ellipse(outcoords,
                                 inputshape,
                                 origin,
                                 inc=0.,
                                 tanpsi=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx = inputshape[0]
    ny = inputshape[1]
    x = -float(xindex - x0)
    y = float(yindex - y0)

    tanpsi0 = tanpsi

    a = ((np.tan(inc) * tanpsi0)**2 - 1.0)
    b = -2. * x * np.sin(inc) * tanpsi0 / (np.cos(inc))**2
    c = y**2 + (x**2 / (np.cos(inc)**2))
    Delta = b**2 - 4. * a * c
    rho_m = (-b - np.sqrt(Delta)) / (2. * a)
    rho_p = (-b + np.sqrt(Delta)) / (2. * a)
    if (rho_p > 0.):
        print("rho_p > 0", rho_p, rho_m, x, y, inc * 180. / np.pi, tanpsi,
              np.arctan(tanpsi) * 180. / np.pi)

    rho = rho_m
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

    H1 = tanpsi0 * rho
    num = x - H1 * np.sin(inc)
    denom = rho * ((1. / np.cos(inc)) - np.tan(inc) * np.sin(inc))
    sintheta = num / denom

    theta = np.arccos(costheta)

    if sintheta < 0:
        theta = 2. * np.pi - theta

    thetaindex = (theta * (float(nx) - 1.) / (2. * np.pi))

    return (rindex, thetaindex)


def polar2cartesian(outcoords, inputshape, origin, inc=0.):
    yindex, xindex = outcoords
    x0, y0 = origin
    nx = inputshape[0]
    ny = inputshape[1]
    x = -float(xindex - x0)
    y = float(yindex - y0)

    rho = np.sqrt((x**2 / (np.cos(inc)**2)) + y**2)
    rindex = rho
    if (rho == 0.):
        costheta = 0.
    else:
        costheta = y / rho

    # theta=np.arccos(costheta)
    theta = np.arctan2((-x / np.cos(inc)), y)
    if (theta < 0):
        theta = theta + 2. * np.pi

    thetaindex = (theta * (float(nx) - 1.) / (2. * np.pi))

    return (rindex, thetaindex)


def carttoconicpolar_nonumba(im, inc, tanpsi):

    (ny, nx) = im.shape
    (i0, j0) = (((float(nx) + 1.) / 2.) - 1., ((float(ny) + 1.) / 2.) - 1.)

    im_polar = sp.ndimage.geometric_transform(im,
                                              cartesian2conicpolar,
                                              order=1,
                                              output_shape=(im.shape[0],
                                                            im.shape[1]),
                                              extra_keywords={
                                                  'inputshape': im.shape,
                                                  'inc': inc,
                                                  'tanpsi': tanpsi,
                                                  'origin': (i0, j0)
                                              })
    return im_polar


def carttoconicpolar(im, inc, tanpsi):

    #(ny, nx) = im.shape
    #im = np.float32(im)
    #im_polar = np.zeros(im.shape, dtype=float32)
    #xoffset_polar = np.zeros(im.shape, dtype=float32)
    #yoffset_polar = np.zeros(im.shape, dtype=float32)

    im_polar = np.zeros(im.shape)
    xoffset_polar = np.zeros(im.shape)
    yoffset_polar = np.zeros(im.shape)

    ConicTransforms.cart2conicpolar_matrix(im,
                                           im_polar,
                                           xoffset_polar,
                                           yoffset_polar,
                                           inc=inc,
                                           tanpsi=tanpsi)

    return im_polar


def carttopolar(im, inc):

    (ny, nx) = im.shape
    (i0, j0) = (((float(nx) + 1.) / 2.) - 1., ((float(ny) + 1.) / 2.) - 1.)

    im_polar = sp.ndimage.geometric_transform(im,
                                              cartesian2polar,
                                              order=1,
                                              output_shape=(im.shape[0],
                                                            im.shape[1]),
                                              extra_keywords={
                                                  'inputshape': im.shape,
                                                  'inc': inc,
                                                  'origin': (i0, j0)
                                              })
    return im_polar


def conicpolartocart(im_polar, inc, tanpsi):

    (ny, nx) = im_polar.shape
    (i0, j0) = (((float(nx) + 1.) / 2.) - 1., ((float(ny) + 1.) / 2.) - 1.)

    im_cart = sp.ndimage.geometric_transform(im_polar,
                                             conicpolar2cartesian_ellipse,
                                             order=1,
                                             output_shape=(im_polar.shape[0],
                                                           im_polar.shape[1]),
                                             extra_keywords={
                                                 'inputshape': im_polar.shape,
                                                 'inc': inc,
                                                 'tanpsi': tanpsi,
                                                 'origin': (i0, j0)
                                             })
    im_cart = np.nan_to_num(im_cart)

    return im_cart


def polartocart(im_polar, inc):

    (ny, nx) = im_polar.shape
    (i0, j0) = (((float(nx) + 1.) / 2.) - 1., ((float(ny) + 1.) / 2.) - 1.)

    im_cart = sp.ndimage.geometric_transform(im_polar,
                                             polar2cartesian,
                                             order=1,
                                             output_shape=(im_polar.shape[0],
                                                           im_polar.shape[1]),
                                             extra_keywords={
                                                 'inputshape': im_polar.shape,
                                                 'inc': inc,
                                                 'origin': (i0, j0)
                                             })
    im_cart = np.nan_to_num(im_cart)

    return im_cart


def exec_prep_files(M):
    """Load the input FITS files prepare the input FITS files
    """
    filename_source = M.filename_source
    workdir = M.workdir
    DumpAllFitsFiles = M.DumpAllFitsFiles

    if (not re.search(r"\/$", workdir)):
        workdir += '/'
        M.workdir = workdir
        print("added trailing back slash to outputdir")

    inbasename = os.path.basename(filename_source)
    filename_fullim = re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim = workdir + filename_fullim

    if (M.DoErrorMap):
        inbasenameerr = os.path.basename(M.filename_errormap)
        filename_fullimerr = re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr = workdir + filename_fullimerr
        filename_fullimw = re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw = workdir + filename_fullimw

    if (M.Verbose):  #
        print("BUILDING WORKDIR AND GENERATING CENTERED IMAGE")

    hdu = pf.open(filename_source)
    hdr0 = hdu[0].header
    if (hdr0['NAXIS'] > 2):
        if (M.Verbose):  #
            print("running cube2im")
        hdu = slice0(filename_source, False)

    im1 = hdu[0].data
    im1 = im1 * M.unitscale
    print("applied unit scale factor:", M.unitscale)
    hdr1 = hdu[0].header

    typicalerror = M.typicalerror
    if (M.InjectNoise):
        print("INJECTING NOISE")
        im1 = im1 + np.random.normal(
            loc=0.0, scale=typicalerror, size=im1.shape)

    if (DumpAllFitsFiles):
        pf.writeto(filename_fullim, im1, hdr1, overwrite=True)

    hdu[0].data = im1

    hduw = False
    if (M.DoErrorMap):
        hduerr = pf.open(M.filename_errormap)
        hdrerr = hduerr[0].header
        if (hdrerr['NAXIS'] > 2):
            hduerr = slice0(M.filename_errormap, False)

        imerr1 = hduerr[0].data
        imerr1 = imerr1 * M.unitscale

        #typicalerror=np.median(imerr1)
        #imerr1[np.where(imerr1 < (typicalerror/3.))] = typicalerror
        #imerr1[np.where(imerr1 > (100.*typicalerror))] = 1E20

        hdrerr1 = hduerr[0].header

        imw0 = 1. / imerr1**2
        imw0 = np.nan_to_num(imw0)
        imw0[np.where(np.fabs(imw0) > 1E10)] = 0.

        if (M.Verbose):
            print("fullim weights: max", np.max(imw0), " min ", np.min(imw0))

        imerr = np.sqrt(1. / imw0)
        imerr = np.nan_to_num(imerr)
        imerr[np.where(np.fabs(imw0) < 1E-30)] = 1E20

        typicalerror = np.median(imerr)

        if (M.Verbose):
            print("resamp errors: max", np.max(imerr), " min ", np.min(imerr))
            print("typicalerror= ", typicalerror, " vs", M.typicalerror)

        if (DumpAllFitsFiles):
            pf.writeto(filename_fullimw, imw0, hdrerr1, overwrite=True)

        hduw = pf.PrimaryHDU()

        hduw.data = imw0
        hduw.header = hdrerr1

    M.Hdu = hdu
    M.Hduw = hduw

    return


def exec_grid_4center(M):
    """prepare the input FITS files: zoom, resample, center, etc...
    """

    workdir = M.workdir
    RA = M.RA
    DEC = M.DEC
    x_center = M.x_center
    y_center = M.y_center
    fieldscale = M.fieldscale  # shrink radial field of view of polar maps by this factor

    hdu = M.Hdu
    hduw = M.Hduw

    filename_source = M.filename_source

    inbasename = os.path.basename(filename_source)
    filename_fullim = re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim = workdir + filename_fullim
    fileout_centered = re.sub('fullim.fits', 'centered.fits', filename_fullim)

    if (M.DoErrorMap):
        inbasenameerr = os.path.basename(M.filename_errormap)
        filename_fullimerr = re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr = workdir + filename_fullimerr
        filename_fullimw = re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw = workdir + filename_fullimw
        fileout_centerederr = re.sub('fullim.fits', 'centered.fits',
                                     filename_fullimerr)
        fileout_centeredw = re.sub('fullim.fits', 'wcentered.fits',
                                   filename_fullimerr)

    hdr1 = hdu[0].header

    if (not (isinstance(y_center, bool))):
        if (not isinstance(RA, float)):
            RA = hdr1['CRVAL1']
            DEC = hdr1['CRVAL2']
            if (M.Verbose):
                print("using center of coords CRVAL1 CRVAL2")

        #RA=RA+(np.sin(x_center*np.pi/180.)*y_center/3600.)/np.cos(DEC*np.pi/180.)
        RA = RA + ((x_center / 3600.) / np.cos(DEC * np.pi / 180.))
        #DEC=DEC+np.cos(x_center*np.pi/180.)*y_center/3600.
        DEC = DEC + ((y_center / 3600.) * np.pi / 180.)
        if (M.Verbose):
            print("RA =", RA)
            print("DEC =", DEC)

    nx = int(hdr1['NAXIS1'] / (M.pixscale_factor * M.fieldscale))
    ny = nx

    if ((nx % 2) == 0):
        nx = nx + 1
        ny = ny + 1

    hdr2 = deepcopy(hdr1)

    hdr2['NAXIS1'] = nx
    hdr2['NAXIS2'] = ny
    hdr2['CRPIX1'] = (nx + 1) / 2
    hdr2['CRPIX2'] = (ny + 1) / 2
    hdr2['CRVAL1'] = RA
    hdr2['CRVAL2'] = DEC
    hdr2['CDELT1'] = M.pixscale_factor * hdr2['CDELT1']
    hdr2['CDELT2'] = M.pixscale_factor * hdr2['CDELT2']

    resamp = gridding(hdu, hdr2, fullWCS=False)
    resamp = np.nan_to_num(resamp)

    if (M.DumpAllFitsFiles):
        fileout_centered = re.sub('fullim.fits', 'centered.fits',
                                  filename_fullim)
        pf.writeto(fileout_centered, resamp, hdr2, overwrite=True)

    hducentered = pf.PrimaryHDU()
    hducentered.data = resamp
    hducentered.header = hdr2

    hduwcentered = False

    if (M.DoErrorMap):
        resampw = gridding(hduw, hdr2, fullWCS=False)
        resampw = np.nan_to_num(resampw)
        resampw[np.where(resampw < 0.)] = 0.

        if (M.Verbose):  #
            print("resamp weights: max", np.max(resampw), " min ",
                  np.min(resampw))

        resamperr = np.sqrt(1. / resampw)

        if (M.DumpAllFitsFiles):
            fileout_centeredw = re.sub('fullim.fits', 'wcentered.fits',
                                       filename_fullimerr)
            pf.writeto(fileout_centeredw, resampw, hdr2, overwrite=True)
            fileout_centerederr = re.sub('fullim.fits', 'centered.fits',
                                         filename_fullimerr)
            pf.writeto(fileout_centerederr, resamperr, hdr2, overwrite=True)

        hduwcentered = pf.PrimaryHDU()
        hduwcentered.data = resampw
        hduwcentered.header = hdr2

    M.Ncorr = (np.pi / (4. * np.log(2))) * M.bmaj * M.bmin / (hdr2['CDELT2'] *
                                                              3600.)**2

    M.Hducentered = hducentered
    M.Hduwcentered = hduwcentered


def exec_conicpolar_expansions(M):

    filename_source = M.filename_source
    workdir = M.workdir
    inc = M.inc
    tanpsi = M.tanpsi
    PlotRadialProfile = M.PlotRadialProfile
    a_min = M.a_min
    a_max = M.a_max
    typicalerror = M.typicalerror
    DumpAllFitsFiles = M.DumpAllFitsFiles
    DoAccr = M.DoAccr
    DoMerid = M.DoMerid
    RestrictAvToRadialDomain = M.RestrictAvToRadialDomain  # set to True is faster but may lead to discontinuities in region averages.

    DoFarSideOnly = M.DoFarSideOnly

    hdu = M.Hducentered
    hduw = M.Hduwcentered

    #typicalerror=ExpectedError

    #cosi=np.cos(inc*np.pi/ 180.)

    inbasename = os.path.basename(filename_source)
    filename_fullim = re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim = workdir + filename_fullim

    if (M.DoErrorMap):
        inbasenameerr = os.path.basename(M.filename_errormap)
        filename_fullimerr = re.sub('.fits', '_fullim.fits', inbasenameerr)
        filename_fullimerr = workdir + filename_fullimerr
        filename_fullimw = re.sub('.fits', '_fullimw.fits', inbasenameerr)
        filename_fullimw = workdir + filename_fullimw

    resamp = hdu.data
    hdr2 = hdu.header
    nx = hdr2['NAXIS1']
    ny = hdr2['NAXIS2']

    if M.Verbose:
        print("M.Ncorr = ", M.Ncorr)

    if (M.DoErrorMap):
        resampw = hduw.data

    if (M.InheritMumap):
        if (M.mumap is None):
            mumap = np.ones(resamp.shape)
        else:
            mumap = M.mumap
    else:
        mumap = np.ones(resamp.shape)

    rotangle = M.PA

    im1rot = ndimage.rotate(resamp, rotangle, reshape=False)
    hdr3 = deepcopy(hdr2)
    im3 = np.double(im1rot)

    if (DumpAllFitsFiles):
        fileout_rotated = re.sub('fullim.fits', 'rotated.fits',
                                 filename_fullim)
        pf.writeto(fileout_rotated, im1rot, hdr2, overwrite=True)

    if (M.DoErrorMap):
        if (np.any(resampw < 0.)):
            print("min / max:", np.min(resampw), np.max(resampw))
            sys.exit("negative  sky weights!!!!")

        im1rotw = ndimage.rotate(resampw, rotangle, reshape=False, order=0)

        im3w = np.double(im1rotw)
        if (DumpAllFitsFiles):
            fileout_rotatedw = re.sub('fullim.fits', 'wrotated.fits',
                                      filename_fullimerr)
            pf.writeto(fileout_rotatedw, im1rotw, hdr2, overwrite=True)

            if (np.any(im1rotw < 0.)):
                print("min / max:", np.min(im1rotw), np.max(im1rotw))
                sys.exit("negative rot sky weights!!!!")

    # #####################################################################
    # take conic polar transforms

    if (M.Verbose):
        print("CARTESIAN2CONICPOLAR TRANSFORM START")
        print("using inc ", inc * np.pi / 180., " tanpsi ", tanpsi)

    im_polar = carttoconicpolar(im3, inc, tanpsi)

    nphis, nrs = im_polar.shape
    if ((nphis != nx) or (nrs != ny)):
        sys.exit("bug")

    hdupolar = pf.PrimaryHDU()
    hdupolar.data = im_polar
    hdrpolar = hdupolar.header
    hdrpolar['CRPIX1'] = 1
    hdrpolar['CRVAL1'] = 0.
    hdrpolar['CDELT1'] = 2. * np.pi / nphis
    hdrpolar['CRPIX2'] = 1
    hdrpolar['CRVAL2'] = 0.
    hdrpolar['CDELT2'] = (hdr3['CDELT2'])
    hdupolar.header = hdrpolar

    fileout_polar = re.sub('fullim.fits', 'polar.fits', filename_fullim)
    if (DumpAllFitsFiles):
        hdupolar.writeto(fileout_polar, overwrite=True)

    if (M.DoErrorMap):

        im_polarw = carttoconicpolar(im3w, inc, tanpsi)

        nphis, nrs = im_polarw.shape

        hdupolarw = pf.PrimaryHDU()
        hdupolarw.data = im_polarw
        hdupolarw.header = hdrpolar

        fileout_polarw = re.sub('fullim.fits', 'wpolar.fits',
                                filename_fullimerr)
        if (DumpAllFitsFiles):
            hdupolarw.writeto(fileout_polarw, overwrite=True)
    else:
        # im_polarw = np.ones(im_polar.shape, dtype=float32) / typicalerror**2
        im_polarw = np.ones(im_polar.shape) / typicalerror**2

    ######################################################################
    # take azimuthal averages on polar maps

    weights = im_polarw.copy()
    #im_Npolcorr = np.ones(im_polarw.shape, dtype=float32)
    im_Npolcorr = np.ones(im_polarw.shape)

    if (np.any(weights < 0.)):
        print("min / max:", np.min(weights), np.max(weights))
        sys.exit("negative polarweights!!")

    im_polar_av = np.copy(im_polar)
    im_polar_rrs = np.zeros(im_polar.shape)
    im_polar_phis = np.zeros(im_polar.shape)

    rrs = 3600. * (np.arange(hdrpolar['NAXIS2']) - hdrpolar['CRPIX2'] +
                   1.0) * hdrpolar['CDELT2'] + hdrpolar['CRVAL2']
    phis = (180. / np.pi) * (
        (np.arange(hdrpolar['NAXIS1']) - hdrpolar['CRPIX1'] + 1.0) *
        hdrpolar['CDELT1'] + hdrpolar['CRVAL1'])
    phis_rad = np.double(
        ((np.arange(hdrpolar['NAXIS1']) - hdrpolar['CRPIX1'] + 1.0) *
         hdrpolar['CDELT1'] + hdrpolar['CRVAL1']))
    KepAmps = np.double(np.zeros(len(rrs)))
    sKepAmps = np.double(np.zeros(len(rrs)))
    AccrAmps = np.double(np.zeros(len(rrs)))
    sAccrAmps = np.double(np.zeros(len(rrs)))
    MeridAmps = np.double(np.zeros(len(rrs)))
    sMeridAmps = np.double(np.zeros(len(rrs)))

    vsysts = np.zeros(hdrpolar['NAXIS2'])

    if (a_min > 0):
        ia_min = np.argmin(np.abs(rrs - a_min))

    if (a_max > 0):
        ia_max = np.argmin(np.abs(rrs - a_max))

    if (M.ComputeSystVelo):
        if (M.DoErrorMap):
            for irrs in range(len(rrs)):
                v0_vec = im_polar[irrs, :]
                w_vec = im_polarw[irrs, :]
                av_v0 = np.sum(w_vec * v0_vec) / np.sum(w_vec)
                av_cosphi = np.sum(w_vec * np.cos(phis_rad)) / np.sum(w_vec)
                KepAmp = np.sum(
                    (v0_vec - av_v0) * w_vec *
                    np.cos(phis_rad)) / np.sum(w_vec *
                                               (np.cos(phis_rad))**2 - w_vec *
                                               np.cos(phis_rad) * av_cosphi)
                vsysts[irrs] = np.sum(
                    w_vec *
                    (v0_vec - KepAmp * np.cos(phis_rad))) / np.sum(w_vec)
        else:
            for irrs in range(len(rrs)):
                v0_vec = im_polar[irrs, :]
                KepAmp = np.sum(
                    (v0_vec - np.average(v0_vec)) * np.cos(phis_rad)) / np.sum(
                        (np.cos(phis_rad))**2 -
                        np.cos(phis_rad) * np.average(np.cos(phis_rad)))
                vsysts[irrs] = np.average(v0_vec - KepAmp * np.cos(phis_rad))

        vsyst = np.asscalar(np.median(vsysts[ia_min:ia_max]))
        sigma_vsyst = np.asscalar(np.std(vsysts[ia_min:ia_max]))
        print("vsyst calculated = ", vsyst, "+-", sigma_vsyst)
        M.vsyst = vsyst
        M.sigma_vsyst = sigma_vsyst
    else:
        vsyst = M.vsyst
        if (M.Verbose):
            print("vsyst from M = ", vsyst)

    TakeAzAv.exec_av(M.DoErrorMap,
                     M.bmaj,
                     M.InheritMumap,
                     M.Verbose,
                     M.PA,
                     M.inc,
                     M.tanpsi,
                     rrs,
                     phis_rad,
                     im_polar,
                     KepAmps,
                     sKepAmps,
                     AccrAmps,
                     sAccrAmps,
                     MeridAmps,
                     sMeridAmps,
                     im_polar_av,
                     im_polar_rrs,
                     im_polar_phis,
                     ia_min,
                     ia_max,
                     im_Npolcorr,
                     vsyst=vsyst,
                     typicalerror=typicalerror,
                     weights=weights,
                     RestrictAvToRadialDomain=RestrictAvToRadialDomain,
                     DoAccr=DoAccr,
                     DoMerid=DoMerid,
                     DoFarSideOnly=DoFarSideOnly,
                     mumap_polarpos=None)

    # SIGNS CALIBRATED ON THE RT TRIALS WIGGLERT
    # beware when flipping across the sky as observer sees the cone with psi < 0, with z<0, where a wind would have v_z < 0 in disk coordinates.
    # /strelka_ssd/simon/wiggleRT/
    sini = np.sin(M.inc)
    cosi = np.cos(M.inc)
    if (np.fabs(M.inc) > (np.pi / 2.)):
        # this fabs to report values for the upper side z>0 of the disk even when retrograde
        # cosi=np.fabs(np.cos(M.inc))
        # best not to do change the results, and instead take care when reporting in RotOrient
        cosi = np.cos(M.inc)

    v_Phi_prof = KepAmps / sini
    sv_Phi_prof = sKepAmps / np.fabs(sini)
    v_Phi_prof = np.nan_to_num(v_Phi_prof)
    sv_Phi_prof = np.nan_to_num(sv_Phi_prof)

    v_R_prof = AccrAmps / sini
    sv_R_prof = sAccrAmps / np.fabs(sini)
    v_R_prof = np.nan_to_num(v_R_prof)
    sv_R_prof = np.nan_to_num(sv_R_prof)

    v_z_prof = -MeridAmps / cosi
    sv_z_prof = sMeridAmps / np.fabs(cosi)
    v_z_prof = np.nan_to_num(v_z_prof)
    sv_z_prof = np.nan_to_num(sv_z_prof)

    ######################################################################
    # now compute chi2 in polar coords

    # im_Npolcorr**2 anticipates correlated pixels in the radial direction
    chi2image = weights * (im_polar - im_polar_av)**2 / im_Npolcorr**2
    chi2image = np.nan_to_num(chi2image)
    deltaChi2 = np.sum(chi2image, axis=1)

    if (np.any(weights < 0.)):
        sys.exit("negative weights!!!!")
    if (np.any(chi2image < 0.)):
        sys.exit("negative chi2!!!!")

    deltaimage = im_polar - im_polar_av
    velodev_med = np.sqrt(np.median(deltaimage[ia_min:ia_max, :]**2))
    velodev_std = np.std(deltaimage[ia_min:ia_max, :])
    velodev_std_vec = np.std(deltaimage, axis=1)
    velodev_std2 = np.std(velodev_std_vec[ia_min:ia_max])

    #varim = deltaimage**2 * weights
    #varvec = np.sum(varim, axis=1)
    #wvec = np.sum(weights, axis=1)
    #mask = (wvec < 1E-10)
    #vec_w_var = varvec  # / wvec)
    #vec_w_var[mask] = 0.
    # #vec_median_w = np.median(weights, axis=1)
    # #vec_typicalerror = np.sqrt(1. / vec_median_w)
    # deltaChi2 = vec_w_var  # / vec_typicalerror**2.)

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
        M.RadialProfile = [
            rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof, v_z_prof,
            sv_z_prof
        ]
    elif (M.DoAccr):
        M.RadialProfile = [rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof]
    else:
        M.RadialProfile = [rrs, v_Phi_prof, sv_Phi_prof]

    if (DumpAllFitsFiles):
        if (M.DoMerid):
            save_prof = np.zeros((hdrpolar['NAXIS2'], 7))
            save_prof[:, 0] = rrs
            save_prof[:, 1] = v_Phi_prof
            save_prof[:, 2] = sv_Phi_prof
            save_prof[:, 3] = v_R_prof
            save_prof[:, 4] = sv_R_prof
            save_prof[:, 5] = v_z_prof
            save_prof[:, 6] = sv_z_prof
        elif DoAccr:
            save_prof = np.zeros((hdrpolar['NAXIS2'], 5))
            save_prof[:, 0] = rrs
            save_prof[:, 1] = v_Phi_prof
            save_prof[:, 2] = sv_Phi_prof
            save_prof[:, 3] = v_R_prof
            save_prof[:, 4] = sv_R_prof
        else:
            save_prof = np.zeros((hdrpolar['NAXIS2'], 3))
            save_prof[:, 0] = rrs
            save_prof[:, 1] = v_Phi_prof
            save_prof[:, 2] = sv_Phi_prof

        fileout_radialprofile = re.sub('fullim.fits', 'radial_profile.dat',
                                       filename_fullim)
        np.savetxt(fileout_radialprofile,
                   save_prof)  # x,y,z equal sized 1D arrays

        fileout_polar_Npolcorr = re.sub('fullim.fits', 'Npolcorr.fits',
                                        filename_fullim)
        hdupolar.data = im_Npolcorr
        hdupolar.writeto(fileout_polar_Npolcorr, overwrite=True)

        fileout_polar_av = re.sub('fullim.fits', 'polar_av.fits',
                                  filename_fullim)
        hdupolar.data = im_polar_av
        hdupolar.writeto(fileout_polar_av, overwrite=True)

        im_polar_diff = im_polar - im_polar_av
        fileout_polar_diff = re.sub('fullim.fits', 'polar_diff.fits',
                                    filename_fullim)
        hdupolar.data = im_polar_diff
        hdupolar.writeto(fileout_polar_diff, overwrite=True)

    ######################################################################
    # Polar average: back to sky  for diff image

    if (M.ComputeSkyImages):

        if (M.Verbose):
            print(
                "CONICPOLAR2CARTESIAN TRANSFORM FOR AZIM AV START im_polar_av")
            print("using inc ", inc * np.pi / 180., "deg tanpsi ", tanpsi)

        #(ny,nx) = im_polar.shape
        #x=np.arange(0,nx)
        #y=np.arange(0,ny)
        #X, Y = np.meshgrid(x, y)
        #rrs_polar= 3600.*(Y-hdrpolar['CRPIX2']+1.0)*hdrpolar['CDELT2']+hdrpolar['CRVAL2']
        #phis_polar= (X-hdrpolar['CRPIX1']+1.0)*hdrpolar['CDELT1']+hdrpolar['CRVAL1']
        #HHs_polar = rrs_polar * tanpsi
        #phis_sky_domain_top=conicpolartocart(phis_polar,inc,tanpsi)
        #rrs_sky_domain_top=conicpolartocart(rrs_polar,inc,tanpsi)
        #HHs_sky_domain_top=conicpolartocart(HHs_polar,inc,tanpsi)
        #phis_sky_domain_top_drot = ndimage.rotate(phis_sky_domain_top, -rotangle, reshape=False,order=0)
        #rrs_sky_domain_top_drot = ndimage.rotate(rrs_sky_domain_top, -rotangle, reshape=False,order=0)
        #HHs_domain_top_drot = ndimage.rotate(HHs_sky_domain_top, -rotangle, reshape=False,order=0)
        #M.diskgeometry={
        #    'HHs_sky_domain_top':HHs_sky_domain_top_drot,
        #    'rrs_sky_domain_top':rrs_sky_domain_top_drot,
        #    'phis_sky_domain_top':phis_sky_domain_top_drot}

        imazim = conicpolartocart(im_polar_av, inc, tanpsi)

        faceoninc = np.pi

        # if (tanpsi < 0.):
        #    faceoninc=np.pi

        if (inc > np.pi / 2.):
            faceoninc = 0.

        imazim_faceon = polartocart(im_polar_av, faceoninc)
        resamp_faceon = polartocart(im_polar, faceoninc)

        im4 = imazim  #gridding(fileout_stretched_av,hdr3)

        #if (M.Verbose):
        #    print("CONICPOLAR2CARTESIAN TRANSFORM FOR AZIM AV  DONE in",
        #          time.time() - start_time)

        if (DumpAllFitsFiles):
            fileout_azimav = re.sub('fullim.fits', 'azim_av.fits',
                                    filename_fullim)
            pf.writeto(fileout_azimav, im4, hdr2, overwrite=True)

        # back to sky - rotate
        im4drot = ndimage.rotate(im4, -rotangle, reshape=False, order=0)

        # diff
        diff_im = resamp - im4drot
        diff_im_faceon = resamp_faceon - imazim_faceon

        im_polar_av_region = im_polar_av.copy()

        im_polar_av_region[0:ia_min] = 0.
        im_polar_av_region[ia_min:ia_max] = 1.
        im_polar_av_region[ia_max:] = 0.

        if M.ExtendRegions:
            if M.iregion == 0:
                print("extending inwards domain of inner region")
                im_polar_av_region[0:ia_min] = 1.
            if M.iregion == (M.n_abins - 2):
                print("extending outwards domain of inner region")
                im_polar_av_region[ia_max:] = 1.

        imazim_region = conicpolartocart(im_polar_av_region, inc, tanpsi)
        imazim_region_drot = ndimage.rotate(imazim_region,
                                            -rotangle,
                                            reshape=False,
                                            order=0)
        imazim_region_faceon = polartocart(im_polar_av_region, faceoninc)

        imazim_rrs = conicpolartocart(im_polar_rrs, inc, tanpsi)
        imazim_rrs_drot = ndimage.rotate(imazim_rrs,
                                         -rotangle,
                                         reshape=False,
                                         order=0)
        imazim_rrs_faceon = polartocart(im_polar_rrs, faceoninc)

        imazim_phis = conicpolartocart(im_polar_phis, inc, tanpsi)
        imazim_phis_drot = ndimage.rotate(imazim_phis,
                                          -rotangle,
                                          reshape=False,
                                          order=0)
        imazim_phis_faceon = polartocart(im_polar_phis, faceoninc)

        mask = np.where(imazim_region_drot > 0.9)
        skychi2 = np.sum(weights[mask] * diff_im[mask]**2) / M.Ncorr
        M.skychi2 = skychi2

        hdudiff = pf.PrimaryHDU()
        hdudiff.data = diff_im
        hdudiff.header = hdr2
        M.Hdudiff = hdudiff

        hdurrs = pf.PrimaryHDU()
        hdurrs.data = imazim_rrs_drot
        hdurrs.header = hdr2
        M.Hdurrs = hdurrs

        hdurrs_faceon = pf.PrimaryHDU()
        hdurrs_faceon.data = imazim_rrs_faceon
        hdurrs_faceon.header = hdr2
        M.Hdurrs_faceon = hdurrs_faceon

        hduphis_faceon = pf.PrimaryHDU()
        hduphis_faceon.data = imazim_phis_faceon
        hduphis_faceon.header = hdr2
        M.Hduphis_faceon = hduphis_faceon

        hduphis = pf.PrimaryHDU()
        hduphis.data = imazim_phis_drot
        hduphis.header = hdr2
        M.Hduphis = hduphis

        hdudiff_faceon = pf.PrimaryHDU()
        hdudiff_faceon.data = diff_im_faceon
        hdudiff_faceon.header = hdr2
        M.Hdudiff_faceon = hdudiff_faceon

        hduresamp_faceon = pf.PrimaryHDU()
        hduresamp_faceon.data = resamp_faceon
        hduresamp_faceon.header = hdr2
        M.Hduresamp_faceon = hduresamp_faceon

        hduregion = pf.PrimaryHDU()
        hduregion.data = imazim_region_drot
        hduregion.header = hdr2
        M.Hduregion = hduregion

        hdumoddrot = pf.PrimaryHDU()
        hdumoddrot.data = im4drot
        hdumoddrot.header = hdr2
        M.Hdumoddrot = hdumoddrot

        hduregion_faceon = pf.PrimaryHDU()
        hduregion_faceon.data = imazim_region_faceon
        hduregion_faceon.header = hdr2
        M.Hduregion_faceon = hduregion_faceon

        if (DumpAllFitsFiles):
            fileout_drotated = re.sub('fullim.fits', 'azim_av_drot.fits',
                                      filename_fullim)
            pf.writeto(fileout_drotated, im4drot, hdr2, overwrite=True)
            fileout_drotated = re.sub('fullim.fits', 'azim_av_drot_diff.fits',
                                      filename_fullim)
            pf.writeto(fileout_drotated, diff_im, hdr2, overwrite=True)
            fileout_drotated_region = re.sub('fullim.fits', 'region_drot.fits',
                                             filename_fullim)
            pf.writeto(fileout_drotated_region,
                       imazim_region_drot,
                       hdr2,
                       overwrite=True)

            fileout_diff_im_faceon = re.sub('fullim.fits', 'diff_faceon.fits',
                                            filename_fullim)
            pf.writeto(fileout_diff_im_faceon,
                       diff_im_faceon,
                       hdr2,
                       overwrite=True)

            fileout_im_faceon = re.sub('fullim.fits', 'resamp_faceon.fits',
                                       filename_fullim)
            pf.writeto(fileout_im_faceon, resamp_faceon, hdr2, overwrite=True)

            fileout_region_faceon = re.sub('fullim.fits', 'region_faceon.fits',
                                           filename_fullim)
            pf.writeto(fileout_region_faceon,
                       imazim_region_faceon,
                       hdr2,
                       overwrite=True)

        ##fileout_proj=re.sub('fullim.fits', 'azim_av_proj.fits', filename_fullim)
        ##if (not OptimOrient):
        ##    pf.writeto(fileout_proj,im4, hdr2, overwrite=True)

    ######################################################################
    # CROSS CHECK INVERSE TRANSFORM

    if (M.XCheckInv):

        print("CONICPOLAR2CARTESIAN TRANSFORM FOR XCHECK INV  START XCheckInv")
        start_time = time.time()
        im_x = conicpolartocart(im_polar, inc, tanpsi)

        print("CONICPOLAR2CARTESIAN TRANSFORM FOR XCHECK INV   DONE in",
              time.time() - start_time)

        fileout_stretched_x = re.sub('fullim.fits', 'stretched_x.fits',
                                     filename_fullim)
        pf.writeto(fileout_stretched_x, im_x, hdr2, overwrite=True)

        #hdr3 = deepcopy(hdr2)
        fileout_proj_x = re.sub('fullim.fits', 'x_proj.fits', filename_fullim)
        #im4_x=gridding(fileout_stretched_x,hdr3)
        pf.writeto(fileout_proj_x, im_x, hdr2, overwrite=True)

        im4_x_drot = ndimage.rotate(im_x, -rotangle, reshape=False)

        fileout_drotated_x = re.sub('fullim.fits', 'x_drot.fits',
                                    filename_fullim)
        pf.writeto(fileout_drotated_x, im4_x_drot, hdr2, overwrite=True)

        fileout_skydiff = re.sub('fullim.fits', 'x_diff.fits', filename_fullim)
        pf.writeto(fileout_skydiff, resamp - im4_x_drot, hdr2, overwrite=True)

    ######################################################################
    # profile plotting

    if (PlotRadialProfile and DumpAllFitsFiles):

        # -----------------------------------------------------------
        # nice fonts
        # -----------------------------------------------------------
        matplotlib.rc('font', family='sans-serif')
        matplotlib.rcParams.update({'font.size': 9})

        plt.figure(figsize=(10, 8))
        axprofile = plt.subplot(111)

        rmax = np.max(rrs)

        plt.setp(axprofile.get_xticklabels(), visible=True)  #, fontsize=6)
        plt.setp(axprofile.get_yticklabels(), visible=True)  #, fontsize=6)
        plt.xlim(0., rmax)
        #plt.ylim(np.min(v_Phi_prof),1.1*np.max(v_Phi_prof))

        plt.ylim(-0.1 * np.max(v_Phi_prof[ia_min:ia_max]),
                 1.1 * np.max(v_Phi_prof[ia_min:ia_max]))
        plt.plot(rrs,
                 v_Phi_prof,
                 color='grey',
                 linewidth=0.1,
                 linestyle='solid')
        plt.fill_between(rrs,
                         v_Phi_prof + sv_Phi_prof,
                         v_Phi_prof - sv_Phi_prof,
                         lw=0.1,
                         color='r',
                         alpha=0.3,
                         interpolate=True,
                         step='mid')

        if (DoMerid):
            plt.plot(rrs,
                     v_z_prof,
                     color='orange',
                     linewidth=1.,
                     linestyle='solid',
                     alpha=0.5,
                     label='v_z')
            plt.plot(rrs,
                     v_R_prof,
                     color='green',
                     linewidth=1.,
                     linestyle='solid',
                     alpha=0.5,
                     label='v_R')
        elif (DoAccr):
            plt.plot(rrs,
                     v_R_prof,
                     color='green',
                     linewidth=1.,
                     linestyle='solid',
                     alpha=0.5,
                     label='v_R')

        #print( np.min(v_Phi_prof))
        #print( np.max(v_Phi_prof))

        plt.plot(rrs,
                 rrs * 0.,
                 color='black',
                 linewidth=0.1,
                 linestyle='solid')

        #plt.plot(rrs,Icutav,color='blue',linewidth=1,linestyle='solid')

        #plt.fill_between(rrs, v_Phi_prof+dispv_Phi_prof, v_Phi_prof-dispv_Phi_prof, lw=0.1,color='b', alpha=0.3, interpolate=True)

        plt.ylabel(r'$\langle |v_{\circ}(r)| \rangle$')

        plt.xlabel(r'$r$ / arcsec')

        fileout_fig = re.sub('fullim.fits', 'fig_profile.pdf', filename_fullim)
        plt.savefig(fileout_fig, bbox_inches='tight')

    ######################################################################

    chi2 = sum(deltaChi2[ia_min:ia_max])
    M.chi2_prev = chi2
    M.velodev_med = velodev_med
    M.velodev_std = velodev_std
    M.velodev_std2 = velodev_std2

    #print( "chi2: ",chi2," typical error ",typicalerror," velodev_std ",velodev_std," velodev_std2 ",velodev_std2,"velodev_med",velodev_med)

    # print( "chi2: ",chi2," typical error ",typicalerror," velodev_med ",velodev_med)

    M.polarchi2 = chi2
    retchi2 = chi2  #  / M.Ncorr ## very aprox correction for correlated pixels because retchi2 is in the polar domain, so NCorr is variable over the polar map.
    return retchi2
