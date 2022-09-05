import sys
import numpy as np
from numba import jit


@jit(nopython=True)
def exec_av(DoErrorMap,
            bmaj,
            InheritMumap,
            Verbose,
            PA,
            inc,
            tanpsi,
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
            ia_min,
            ia_max,
            im_Npolcorr,
            vsyst=0.,
            typicalerror=1.,
            weights=None,
            RestrictAvToRadialDomain=False,
            DoAccr=False,
            DoMerid=False,
            DoFarSideOnly=False,
            mumap_polarpos=None):

    for irrs in range(len(rrs)):

        KepAmps[irrs] = 0.
        sKepAmps[irrs] = 0.
        AccrAmps[irrs] = 0.
        sAccrAmps[irrs] = 0.
        MeridAmps[irrs] = 0.
        sMeridAmps[irrs] = 0.
        im_polar_av[irrs, :] = 0.
        if (((irrs < ia_min) or (irrs > ia_max)) and RestrictAvToRadialDomain):
            continue

        v0_vec = im_polar[irrs, :] - vsyst
        # print("type v0_vec", v0_vec.dtype) ## DEV 
        if (DoErrorMap):
            w_vec = weights[irrs, :]
        else:
            w_vec = np.ones(len(v0_vec), dtype='float32') / typicalerror**2

        # print("type w_vec", w_vec.dtype) ## DEV
        mask = (w_vec < 1E-10)
        w_vec_nozeros = w_vec
        w_vec_nozeros[mask] = 1E-10
        err_vec = (1.0 / np.sqrt(w_vec_nozeros))
        err_vec[mask] = 1E20  # np.inf
        w_vec[mask] = 0.

        thisradius = rrs[irrs]  #arcsec
        # cosi=np.cos(inc)
        #Nind=2.*np.pi*thisradius * np.fabs(cosi) /(bmaj) # aprox number of beams at each radius
        Nind = 2. * np.pi * thisradius / (
            bmaj)  # aprox number of beams at each radius
        if (Nind < 1.):
            Nind = 1.

        Nsum = len(w_vec)
        Npolcorr = Nsum / Nind
        #print("irrs ",irrs,Nsum, Nind, Npolcorr)
        if Npolcorr < 1:
            Npolcorr = 1.
        im_Npolcorr[irrs, :] = Npolcorr

        if ((np.sum(w_vec) == 0.)
                or (np.sum(w_vec * (np.cos(phis_rad))**2)**2 == 0.)):
            KepAmp = 0.
            sKepAmp = 1E20
            AccrAmp = 0.
            sAccrAmp = 1E20
            MeridAmp = 0.
            sMeridAmp = 1E20
        else:
            if (InheritMumap):
                print("coucou")  # DEV DEV
                # mumap_vec = mumap_polarpos[irrs, :]
                # if (DoAccr):
                #     sys.exit("Not yet programmed DoAccr with mumap")
                # else:
                #     KepAmp = np.sum(w_vec * v0_vec * mumap_vec *
                #                     np.cos(phis_rad)) / np.sum(
                #                         w_vec * mumap_vec *
                #                         (np.cos(phis_rad))**2)
                #     AccrAmp = 0.
            else:
                Cramer = True

                if (DoMerid and Cramer):
                    if (np.any(phis_rad < 0.)):
                        print("BUG BUG phi should not be negative")
                        # sys.exit("negative phis_rad")

                    maskphis = (phis_rad >= 0.)
                    if DoFarSideOnly:
                        maskphis = ((phis_rad > 0.) & (phis_rad < np.pi))
                    subw_vec = w_vec[maskphis]
                    subv0_vec = v0_vec[maskphis]
                    subphis_rad = phis_rad[maskphis]

                    sinphis = np.sin(subphis_rad)
                    cosphis = np.cos(subphis_rad)
                    a_1 = np.sum(subw_vec * cosphis**2)
                    b_1 = np.sum(subw_vec * sinphis * cosphis)
                    c_1 = np.sum(subw_vec * cosphis)
                    d_1 = np.sum(subw_vec * cosphis * subv0_vec)
                    vard_1 = np.sum(subw_vec * cosphis**2)
                    a_2 = b_1
                    b_2 = np.sum(subw_vec * sinphis**2)
                    c_2 = np.sum(subw_vec * sinphis)
                    d_2 = np.sum(subw_vec * sinphis * subv0_vec)
                    vard_2 = np.sum(subw_vec * sinphis**2)
                    a_3 = c_1
                    b_3 = c_2
                    c_3 = np.sum(subw_vec)
                    d_3 = np.sum(subw_vec * subv0_vec)
                    vard_3 = np.sum(subw_vec)

                    detM = a_1 * (b_2 * c_3 - b_3 * c_2) - b_1 * (
                        a_2 * c_3 - a_3 * c_2) + c_1 * (a_2 * b_3 - a_3 * b_2)
                    if (detM == 0.):
                        #print("singular matrix")
                        continue

                    A = (d_1 * (b_2 * c_3 - b_3 * c_2) + d_2 *
                         (c_1 * b_3 - b_1 * c_3) + d_3 *
                         (b_1 * c_2 - c_1 * b_2)) / detM
                    sigmaA = np.sqrt(vard_1 *
                                     (b_2 * c_3 - b_3 * c_2)**2 + vard_2 *
                                     (c_1 * b_3 - b_1 * c_3)**2 + vard_3 *
                                     (b_1 * c_2 - c_1 * b_2)**2) / detM

                    B = (d_1 * (a_3 * c_2 - a_2 * c_3) + d_2 *
                         (a_1 * c_3 - c_1 * a_3) + d_3 *
                         (c_1 * a_2 - a_1 * c_2)) / detM

                    #Bcheck = (a_1 * ( d_2 *c_3 - d_3 * c_2) - d_1 * (a_2 * c_3 - a_3 * c_2) + c_1 * (a_2 * d_3 - a_3 * d_2)) /detM

                    sigmaB = np.sqrt(vard_1 *
                                     (a_3 * c_2 - a_2 * c_3)**2 + vard_2 *
                                     (a_1 * c_3 - c_1 * a_3)**2 + vard_3 *
                                     (c_1 * a_2 - a_1 * c_2)**2) / detM

                    C = (d_1 * (a_2 * b_3 - a_3 * b_2) + d_2 *
                         (b_1 * a_3 - a_1 * b_3) + d_3 *
                         (a_1 * b_2 - b_1 * a_2)) / detM
                    sigmaC = np.sqrt(vard_1 *
                                     (a_2 * b_3 - a_3 * b_2)**2 + vard_2 *
                                     (b_1 * a_3 - a_1 * b_3)**2 + vard_3 *
                                     (a_1 * b_2 - b_1 * a_2)**2) / detM

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
                    sum_w = np.sum(w_vec)
                    sum_wv0 = np.sum(w_vec * v0_vec)
                    sum_wsinphi = np.sum(w_vec * np.sin(phis_rad))
                    sum_wcosphi = np.sum(w_vec * np.cos(phis_rad))
                    sinphis = np.sin(phis_rad)
                    cosphis = np.cos(phis_rad)
                    subA_num = (np.sum(w_vec * sinphis**2 -
                                       w_vec * sinphis * sum_wsinphi / sum_w) /
                                np.sum(w_vec * sinphis * cosphis -
                                       w_vec * cosphis * sum_wsinphi / sum_w))

                    A_num = (
                        np.sum(v0_vec * w_vec * sinphis -
                               w_vec * sinphis * sum_wv0 / sum_w) -
                        subA_num * np.sum(v0_vec * w_vec * cosphis -
                                          w_vec * cosphis * sum_wv0 / sum_w))

                    subA_denom1_num1 = np.sum(w_vec * cosphis**2 - w_vec *
                                              cosphis * sum_wcosphi / sum_w)
                    subA_denom1_denom1 = np.sum(w_vec * cosphis * sinphis -
                                                w_vec * cosphis * sum_wsinphi /
                                                sum_w)

                    A_denom = (np.sum((w_vec * sinphis * cosphis -
                                       w_vec * sinphis * sum_wcosphi / sum_w) -
                                      (w_vec * sinphis**2 -
                                       w_vec * sinphis * sum_wsinphi / sum_w) *
                                      subA_denom1_num1 / subA_denom1_denom1))

                    A = A_num / A_denom

                    varA_num = np.sum(
                        sinphis**2 * w_vec +
                        w_vec**2 * sinphis**2 / sum_w) + subA_num**2 * np.sum(
                            cosphis**2 * w_vec + w_vec**2 * cosphis**2 / sum_w)
                    sigma_A = np.sqrt(varA_num / A_denom**2)

                    B_denom = np.sum(w_vec * cosphis * sinphis -
                                     w_vec * cosphis * sum_wsinphi / sum_w)
                    subB_num = np.sum(w_vec * cosphis**2 -
                                      w_vec * cosphis * sum_wcosphi / sum_w)
                    B = ((np.sum(v0_vec * w_vec * cosphis - w_vec * cosphis *
                                 sum_wv0 / sum_w) - A * subB_num) / B_denom)
                    varB = (np.sum(w_vec * cosphis**2 +
                                   w_vec**2 * cosphis**2 / sum_w) +
                            sigma_A**2 * subB_num**2) / B_denom**2
                    sigma_B = np.sqrt(varB)

                    C_denom = sum_w * abs(np.cos(inc))
                    C = np.sum(w_vec *
                               (v0_vec - A * cosphis - B * sinphis)) / C_denom

                    varC = np.sum(w_vec**2 *
                                  (err_vec**2 + sigma_A**2 * cosphis**2 +
                                   sigma_B**2 + sinphis**2)) / C_denom**2
                    sigma_C = np.sqrt(varC)

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
                    subsum1 = np.sum(
                        w_vec * v0_vec * np.cos(phis_rad)) / np.sum(
                            w_vec * (np.cos(phis_rad))**2)
                    varsubsum1 = np.sum(
                        (np.sqrt(w_vec) * np.cos(phis_rad)) /
                        np.sum(w_vec * (np.cos(phis_rad))**2)**2)
                    numerator = np.sum(w_vec * np.sin(phis_rad) *
                                       (v0_vec - np.cos(phis_rad) * subsum1))
                    varnumerator = np.sum(
                        w_vec**2 * np.sin(phis_rad)**2 *
                        ((err_vec)**2 + np.cos(phis_rad)**2 * varsubsum1))
                    subsum2 = np.sum(
                        w_vec * np.sin(phis_rad) * np.cos(phis_rad)) / np.sum(
                            w_vec * (np.cos(phis_rad))**2)
                    denom = np.sum(w_vec * ((np.sin(phis_rad))**2 - subsum2))
                    AccrAmp = numerator / denom
                    sAccrAmp = np.sqrt(varnumerator) / denom
                    KepAmp = np.sum(
                        w_vec * (v0_vec - AccrAmp * np.sin(phis_rad)) *
                        np.cos(phis_rad)) / np.sum(w_vec *
                                                   (np.cos(phis_rad))**2)

                    sKdenom = np.sum(w_vec * (np.cos(phis_rad))**2)
                    varKnum = np.sum(
                        w_vec**2 *
                        (err_vec**2 + sAccrAmp**2 * np.sin(phis_rad)**2) *
                        np.cos(phis_rad)**2)

                    sKnum = np.sqrt(varKnum)

                    sKepAmp = sKnum / sKdenom

                    #print("sKepAmp>",rrs[irrs],KepAmp, varKnum, sKnum,sKdenom, sKepAmp, sAccrAmp)

                    sAccrAmp *= np.sqrt(Npolcorr)
                    sKepAmp *= np.sqrt(Npolcorr)
                    MeridAmp = 0.
                    sMeridAmp = 1E20

                else:

                    denom = np.sum(w_vec * (np.cos(phis_rad))**2)
                    numerator = np.sum(w_vec * v0_vec * np.cos(phis_rad))
                    KepAmp = numerator / denom

                    sdenom = np.sqrt(np.sum(w_vec * (np.cos(phis_rad))**2))
                    sKepAmp = 1. / sdenom

                    #print( "bmaj = ",bmaj, " thisradius  ",thisradius, " inc  ", inc ,"Nind ",Nind," Ncorr ",Ncorr, " \n")
                    sKepAmp *= np.sqrt(Npolcorr)

                    AccrAmp = 0.
                    sAccrAmp = 1E20
                    MeridAmp = 0.
                    sMeridAmp = 1E20

        # else:
        # v0_vec=im_polar[irrs,:]-vsyst
        # KepAmp = np.sum(v0_vec * np.cos(phis_rad)) / np.sum((np.cos(phis_rad))**2)

        if (Verbose and (KepAmp < 0.) and (irrs == int(len(rrs) / 4.))):
            print("KepAmp negative, wrong PA: KepAmp=", KepAmp, " PA: ", PA,
                  " inc", inc * np.pi / 180., "deg  tanpsi", tanpsi)
            # print((np.array(zip(v0_vec, np.cos(phis_rad)))))

        KepAmps[irrs] = KepAmp
        sKepAmps[irrs] = sKepAmp
        AccrAmps[irrs] = AccrAmp
        sAccrAmps[irrs] = sAccrAmp
        MeridAmps[irrs] = MeridAmp
        sMeridAmps[irrs] = sMeridAmp

        if (DoMerid):
            # v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(phis_rad) * MeridAmp*np.cos(inc)
            v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(
                phis_rad) + MeridAmp
            im_polar_av[irrs, :] = v0_vec_av + vsyst
        elif (DoAccr):
            v0_vec_av = KepAmp * np.cos(phis_rad) + AccrAmp * np.sin(phis_rad)
            im_polar_av[irrs, :] = v0_vec_av + vsyst
        else:
            v0_vec_av = KepAmp * np.cos(phis_rad)
            im_polar_av[irrs, :] = v0_vec_av + vsyst
