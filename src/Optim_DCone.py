import os
import os.path
import re
from copy import deepcopy
import numpy as np
from astropy.wcs import WCS
from scipy import optimize
from time import gmtime, strftime
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.colors as colors
import scipy.optimize as op

from ConeRot.src.funcs_DConeMaps import *
from ConeRot.src.DConeMaps import Model
from ConeRot.src.funcs_Optim_DCone import exec_ConjGrad_1region
import ConeRot.src.KineSummary as KineSummary
from ConeRot.src.funcs_Optim_DCone import exec_emcee
from ConeRot.src.funcs_Optim_DCone import exec_emcee
from ConeRot.src.funcs_Optim_DCone import exec_Regions

#def pass_model(Mpass):
#    global M
#    M=Mpass

######################################################################


class OptimModel():

    #def __init__(self,M,PrintOptimStatus=True): #,DoConjGrad=False, RunMCMC=False
    #    self.PrintOptimStatus=M.PrintOptimStatus
    #    #self.DoConjGrad=DoConjGrad
    #    #self.RunMCMC=RunMCMC

    def __init__(
        self,
        M,
        RunMCMC=False,
        Nit=1,  #MCMC iterations
        nwalkers=-1,
        burn_in=50,
        n_cores_MCMC=1,
        BlindMCMC=False,
        #a_min_regions=0.2,
        #a_max_regions=0.3,
        #n_abins=12,
        #StoreRegions=False
        n_cores_regions=4):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            print("OptimModel setting ", a_attribute, " to ",
                  initlocals[a_attribute])
            setattr(self, a_attribute, initlocals[a_attribute])

        # print( "opening log:",M.workdir+M.filelog)
        # fout=open(M.workdir+M.filelog,"w+")
        # M.fout=fout

    def ConjGrad_1region(self, M):
        M.DumpAllFitsFiles = False
        M.Verbose = False
        M.ComputeSkyImages = False
        return exec_ConjGrad_1region(M, self)

    def RecoverConjGrad_1region(self, M):
        names = list(map((lambda x: x[0]), M.domain))
        bnds = list(map((lambda x: x[1]), M.domain))
        nvar = len(names)

        result_ml = np.load(M.workdir + 'result_ml.dat.npy')
        print("result_ml is ", result_ml)
        print("LOG:", M.fout)
        M.fout.write("Global: \n")
        for iparam in range(nvar):
            print(names[iparam], "->", result_ml[iparam])
            setattr(M, names[iparam], result_ml[iparam])
            M.fout.write(names[iparam] + "-> %.6f " % (result_ml[iparam]))

        M.fout.write("\n")
        M.DumpAllFitsFiles = True
        M.prep_files()
        M.grid_4center()
        M.ComputeSkyImages = True

        chi2 = M.conicpolar_expansions()
        print("chi2=", chi2)
        M.fout.write("chi2=%.6e\n" % (chi2))
        inbasename = os.path.basename(M.filename_source)
        inbasename = re.sub('.fits', '', inbasename)
        inbasename = M.workdir + inbasename
        fileout = inbasename + '_fig_summary.pdf'

        nplots = 3
        if (M.DoErrorMap):
            inbasenameerrormap = os.path.basename(M.filename_errormap)
            inbasenameerrormap = re.sub('.fits', '', inbasenameerrormap)
            inbasenameerrormap = M.workdir + inbasenameerrormap
            nplots = 4
        else:
            inbasenameerrormap = False

        KineSummary.exec_summary(inbasename,
                                 fileout,
                                 vsyst=M.vsyst,
                                 basename_errormap=inbasenameerrormap,
                                 nplots=nplots)

    def emcee(self, M):
        M.DumpAllFitsFiles = False
        M.Verbose = False
        M.PrintOptimStatus = False
        M.ComputeSkyImages = False

        result_ml = np.load(M.workdir + 'result_ml.dat.npy')

        retvals = exec_emcee(M, result_ml, True, self)
        return retvals

    def RecoverMCMC(self, M):
        names = list(map((lambda x: x[0]), M.domain))
        bnds = list(map((lambda x: x[1]), M.domain))
        nvar = len(names)

        result_ml = np.load(M.workdir + 'bestparams.dat.npy')

        #self.RunMCMC=False

        # exec_emcee(M,self.Nit,self.nwalkers,result_ml,self.n_cores_MCMC,False,self) #self.RunMCMC
        exec_emcee(M, result_ml, False, self)  #self.RunMCMC
        M.fout.write("Global MCMC: \n")
        print("MCMC best params  is ", result_ml)
        for iparam in list(range(nvar)):
            print(names[iparam], "->", result_ml[iparam])
            setattr(M, names[iparam], result_ml[iparam])
            M.fout.write(names[iparam] + "-> %.6f " % (result_ml[iparam]))

        M.fout.write("\n")

        mcmc_results = np.load(M.workdir + 'mcmc_results.dat.npy')
        mcmc_results = mcmc_results.tolist()

        logstring = "emcee posterior\n"
        for iparam in list(range(nvar)):
            #strparams=names[iparam]+" -> %.6f %.6f %.6f " % mcmc_results[iparam]
            strparams = names[iparam] + " -> {0:.6f} {1:.6f} {2:.6f} ".format(
                *mcmc_results[iparam])
            logstring = logstring + strparams + "\n"

        M.fout.write(logstring)

        M.DumpAllFitsFiles = True
        #M.Optim=False
        M.prep_files()
        M.grid_4center()
        M.ComputeSkyImages = True

        chi2 = M.conicpolar_expansions()
        print("chi2=", chi2)
        print("polar chi2=", M.polarchi2)
        print("sky chi2=", M.skychi2)
        M.fout.write("chi2=%.3e\n" % (M.skychi2))

        inbasename = os.path.basename(M.filename_source)
        inbasename = re.sub('.fits', '', inbasename)
        inbasename = M.workdir + inbasename
        fileout = inbasename + '_emcee_fig_summary.pdf'

        nplots = 3
        if (M.DoErrorMap):
            inbasenameerrormap = os.path.basename(M.filename_errormap)
            inbasenameerrormap = re.sub('.fits', '', inbasenameerrormap)
            inbasenameerrormap = M.workdir + inbasenameerrormap
            nplots = 4
        else:
            inbasenameerrormap = False

        KineSummary.exec_summary(inbasename,
                                 fileout,
                                 vsyst=M.vsyst,
                                 basename_errormap=inbasenameerrormap,
                                 nplots=nplots)

    def Regions(self, M):
        M.DumpAllFitsFiles = False
        M.Verbose = False
        M.PrintOptimStatus = False

        exec_Regions(M, self)  # self.n_cores_regions)
        #exec_Regions(M,self.a_min_regions,self.a_max_regions,self.n_abins,self.n_cores_regions,self.StoreRegions)
