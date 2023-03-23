import math as ma
from numba import jit

@jit(nopython=True)
def cartesian2conicpolar(rindex,
                         thetaindex,
                         nx,
                         ny,
                         origin,
                         inc=0.,
                         tanpsi=0.):

    x0, y0 = origin
    theta = thetaindex * 2. * ma.pi / nx
    y = rindex * ma.cos(theta)
    height = tanpsi * rindex
    x = (rindex*ma.sin(theta))/ma.cos(inc) + \
             (height - (rindex*ma.sin(theta))*ma.tan(inc))*ma.sin(inc)

    ix = -x + x0
    iy = y + y0

    #if ix < 0.:
    #    ix = 0
    #if ix > nx-1:
    #    ix = nx-1.
    #if iy < 0.:
    #    iy = 0
    #if iy > ny-1:
    #    iy = ny-1.

    return (iy, ix)


@jit(nopython=True)
def cart2conicpolar_matrix(im,
                           im_polar,
                           xindex_polar,
                           yindex_polar,
                           inc=10.,
                           tanpsi=0.1):

    nx = im.shape[0]
    ny = im.shape[1]
    (i0, j0) = (((float(nx) + 1.) / 2.) - 1., ((float(ny) + 1.) / 2.) - 1.)
    origin = (i0, j0)
    for ir in range(nx):
        for iphi in range(ny):
            yindex, xindex = cartesian2conicpolar(ir,
                                                  iphi,
                                                  nx,
                                                  ny,
                                                  origin,
                                                  inc=inc,
                                                  tanpsi=tanpsi)
            xindex_polar[ir, iphi] = xindex
            yindex_polar[ir, iphi] = yindex

    Bilin = True
    for ir in range(nx):
        for iphi in range(ny):
            ix = xindex_polar[ir, iphi]
            iy = yindex_polar[ir, iphi]
            i = ma.floor(ix)
            j = ma.floor(iy)
            OutOfBounds = False
            if (i < 0):
                i = 0
                OutOfBounds = True
            if (j < 0):
                j = 0
                OutOfBounds = True
            if i >= (nx - 2):
                i = nx - 2
                OutOfBounds = True
            if j >= (ny - 2):
                j = ny - 2
                OutOfBounds = True
            if OutOfBounds:
                f00 = 0.
            else:
                f00 = im[j, i]

            im_polar[ir, iphi] = f00

            if Bilin:
                dx = ix - i
                dy = iy - j
                ip1 = i + 1
                jp1 = j + 1
                if i >= (nx - 1):
                    ip1 = nx - 1
                if i >= (ny - 1):
                    jp1 = ny - 1

                f00 = im[j, i]
                f10 = im[j, ip1]
                f01 = im[jp1, i]
                f11 = im[jp1, ip1]
                if OutOfBounds:
                    fbilin = 0.
                else:
                    fbilin = f00 * (1. - dx) * (1. - dy) + f10 * dx * (
                        1. - dy) + f01 * (1. - dx) * dy + f11 * dx * dy

                im_polar[ir, iphi] = fbilin
