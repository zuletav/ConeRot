import sys
import os
import os.path
import numpy as np
from astropy.io import fits as pf
import re


include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import  ConeRot.DConeMaps  as DConeMaps
from  ImUtils.FixHeads import fixhead_3 
import  ConeRot.Optim_DCone as Optim_DCone


class Setup():
    def __init__(self,
                 filename_source='',
                 filename_errormap='',
                 workdir='',
                 DoErrorMap=False,
                 typicalerror=0.1,
                 ComputeSystVelo=False,
                 vsyst=0.,
                 sigma_vsyst=0., #uncertainty in the inferred value of vsyst
                 fieldscale=1.,
                 pixscale_factor=1.,
                 unitscale=1.,
                 PA=0.,
                 inc=1.,
                 tanpsi=0.1,
                 rangePA=10.,
                 rangeinc=20.*np.pi/180.,
                 rangetanpsi=0.4,
                 ClearWorkDir=True,
                 a_min=1.0,
                 a_max=2.0,
                 DoRegions=False,
                 InheritGlobalInit=False, # to force same initial conditions for conjgrad, in all regions
                 a_min_regions=1.0,
                 a_max_regions=2.0,
                 n_abins=10,
                 DoAccr=False,
                 DoAccr_fixPAinc=False,
                 DoMerid=False,
                 DoMerid_fixPAinc=False,
                 DoExec=True,
                 DoFixOrient=True,
                 DumpAllFitsFiles=False,
                 fout=False,
                 x_center=0.,
                 y_center=0.,
                 bmaj=1.,
                 bmin=1.,
                 DoConjGrad=False,
                 DoMinuit=False,
                 RunMCMC=False,
                 RecoverMCMC=False,
                 n_cores_MCMC=10,
                 Nit=140,
                 burn_in=70,
                 BlindMCMC=False, # if Falses uses conjgrad optim for init
                 nwalkers=-1,  # <= 0 value will use defaults
                 domain=(),
                 RA=False,
                 DEC=False,
                 InjectNoise=False,
                 DoDCone=False,
                 InheritMumap=False,  # pass mumap from a previous orientation - used as weights in KepAmps
                 StoreRegions=False,
                 DoFarSideOnly=False,
                 exec_master_script='exec_master.py'):


        initlocals=locals()
        initlocals.pop('self')
        passargs={}
        for a_attribute in initlocals.keys():
            print( "MasterDConeMaps  setting ",a_attribute," to ",initlocals[a_attribute])
            setattr(self,a_attribute,initlocals[a_attribute])
            passargs['a_attribute']=initlocals[a_attribute]


        self.PA0=self.PA
        self.inc0=self.inc
        self.tanpsi0=self.tanpsi
        self.n_cores_regions=self.n_abins-1 # -1

        

    def Run(self):

        filename_in=self.filename_source
        fixhead_3(self.filename_source,filename_in)
        if self.DoErrorMap:
            filename_in_errormap=self.filename_errormap
            fixhead_3(self.filename_errormap,filename_in_errormap)
        else: 
            filename_in_errormap=False
             
             
             
        if (self.ClearWorkDir):
            os.system("rm -rf  "+self.workdir)

        os.system("mkdir "+self.workdir)

        os.system("rsync -va "+self.exec_master_script+" "+self.workdir)
        os.system("tar cvfz "+self.workdir+"ball_conemaps.tgz "+include_path+"conemaps   ")


        
        M=DConeMaps.Model(VerboseInit=False)

        for a_attribute in M.__dict__.keys():
            if (a_attribute in self.__dict__.keys()):
                print( "setting M ",a_attribute," to ",self.__dict__[a_attribute])
                setattr(M,a_attribute,self.__dict__[a_attribute])

                
        # #####################################################################
                
        # if (not M.fout):
        #    fout=open(M.workdir+M.filelog,"w+")
        fout=open(M.workdir+M.filelog,"a+")
        M.fout=fout
        self.fout=fout
        M.fout.write("Init: \n")
        M.fout.write("PA= "+str(M.PA)+" ")
        M.fout.write("inc= "+str(M.inc*180./np.pi)+" ")
        M.fout.write("tanpsi= "+str(M.tanpsi)+" ")


        OptimM=Optim_DCone.OptimModel(M)

        print( "self.DoConjGrad", self.DoConjGrad)

        if self.DoConjGrad:
            if (self.ComputeSystVelo):
                M.ComputeSystVelo=True
                OptimM.ConjGrad_1region(M)
                M.ComputeSystVelo=False
                print( "Calculated systemic velocity:",M.vsyst)
                M.fout.write("Calculated systemic velocity:%.6f +- %.6f \n" % (M.vsyst, M.sigma_vsyst))
            else:
                M.fout.write("Input systemic velocity:%.6f\n" % (M.vsyst))

            self.vsyst=M.vsyst
            self.sigma_vsyst=M.sigma_vsyst
            

            OptimM.ConjGrad_1region(M)

            print( ">>>  velodev_med "+str(M.velodev_med))
            print( ">>>  velodev_std "+str(M.velodev_std))
            print( ">>>  velodev_std2 "+str(M.velodev_std2))
            
            if ((not self.DoErrorMap) and (M.velodev_med > self.typicalerror*3.)):
                # sys.exit("RERUN with correct typicalerror, recommend to use at least "+str(M.velodev_med))
                print(("RERUN with correct typicalerror, recommend to use at least "+str(M.velodev_med)))

        OptimM.RecoverConjGrad_1region(M)
        
        # #####################################################################
        # emcee
        nparams=len(self.domain)
        print("nparams",nparams)
        if (self.nwalkers <= 0):
            self.nwalkers=int(10*nparams)
        
        if (self.RunMCMC and  (os.path.isdir(self.workdir))):
            OptimM=Optim_DCone.OptimModel(M,RunMCMC=True,Nit=self.Nit,nwalkers=self.nwalkers,n_cores_MCMC=self.n_cores_MCMC,burn_in=self.burn_in)
            print("MasterDConeMaps: calling OptimM.emcee with self.n_cores_MCMC=",self.n_cores_MCMC)
            OptimM.emcee(M)

        if (self.RecoverMCMC):
            OptimM=Optim_DCone.OptimModel(M,RunMCMC=True,Nit=self.Nit,nwalkers=self.nwalkers,burn_in=self.burn_in)
            OptimM.RecoverMCMC(M)

        # #####################################################################
        # REGIONS

        if (self.DoRegions):
            OptimM=Optim_DCone.OptimModel(M,n_cores_regions=self.n_cores_regions,RunMCMC=True,Nit=self.Nit,nwalkers=self.nwalkers,n_cores_MCMC=self.n_cores_MCMC,burn_in=self.burn_in)
            OptimM.Regions(M)

        print("Master.Run closing fout",fout)
        fout.close()

        return


    
    def RunFixOrient(self,ForceGlobalOrient=False, Force_allradsPA=0., Force_allradsinc=0.):

        self.DoAccr=self.DoAccr_fixPAinc
        self.DoMerid=self.DoMerid_fixPAinc
        if (not ForceGlobalOrient):
            file_log=self.workdir+'log_output.txt'
            print("loading file_log",file_log)
            fin= open(file_log,"r")
            log_output=fin.readlines()
            fin.close
            
            AllRads=False
            AllRadsMCMC=False
            emcee_posterior=False

            for aline in log_output:
                if AllRads:
                    matches = re.search("^PA-> (.*) inc-> (.*) tanpsi-> (.*) $",aline)
                    allradsPA=float(matches.group(1))
                    allradsinc=float(matches.group(2)) # *180./np.pi
                    allradstanpsi=float(matches.group(3))
                    AllRads=False
                if AllRadsMCMC:
                    matches = re.search("^PA-> (.*) inc-> (.*) tanpsi-> (.*) $",aline)
                    allradsPAMCMC=float(matches.group(1))
                    allradsincMCMC=float(matches.group(2)) # *180./np.pi
                    allradstanpsiMCMC=float(matches.group(3))

                    allradsPA=allradsPAMCMC
                    allradsinc=allradsincMCMC

                    AllRadsMCMC=False
                if "Global" in aline:
                    AllRads=True
                if "Global MCMC" in aline:
                    AllRadsMCMC=True
                #if "chi2" in aline:
                #    print("chi2 in aline",aline)
                #    matches = re.search("^chi2=(.*)$",aline)
                #    chi2allrads=float(matches.group(1))

                if emcee_posterior:
                    matches = re.search("^PA ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    #matches = re.search("^PA ->\s+(\w+)\s+(\w+)\s+(\w+)\s*$",aline)
                    if (matches):
                        aPA_post=matches.group(1)
                        aPA_upsigma=matches.group(2)
                        aPA_downsigma=matches.group(3)
                        PA_MCMC=[float(aPA_post),float(aPA_upsigma),float(aPA_downsigma)]
                        allradsPA=PA_MCMC[0]
                    matches = re.search("^inc ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    if (matches):
                        ainc_post=matches.group(1)
                        ainc_upsigma=matches.group(2)
                        ainc_downsigma=matches.group(3)
                        inc_MCMC=[float(ainc_post),float(ainc_upsigma),float(ainc_downsigma)]
                        allradsinc=inc_MCMC[0]
                    matches = re.search("^tanpsi ->\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    if (matches):
                        atanpsi_post=matches.group(1)
                        atanpsi_upsigma=matches.group(2)
                        atanpsi_downsigma=matches.group(3)
                        tanpsi_MCMC=[float(atanpsi_post),float(atanpsi_upsigma),float(atanpsi_downsigma)]
                        allradtanpsi=tanpsi_MCMC[0]
                    emcee_posterior=False
                if ("emcee posterior" in aline):
                    emcee_posterior=True
                if ("Regions" in aline):
                    emcee_posterior=False
                    break
                    


            PA=allradsPA
            inc=allradsinc
            tanpsi=allradstanpsi

        else:
            print(">>>>>> Doing Fix Orient, with forced global orientation, at PA=",Force_allradsPA," inc=",Force_allradsinc)

            PA=Force_allradsPA
            inc=Force_allradsinc
            tanpsi=self.tanpsi

        
        # print( "global chi2:",chi2allrads)
        print( "using global PA:",PA)
        print( "using global inc:",inc)
        print( "using global tanpsi:",tanpsi)


        self.PA=PA
        self.inc=inc
        self.tanpsi=tanpsi
        self.nwalkers = 7 # 
        self.domain=( ('tanpsi',(tanpsi-self.rangetanpsi/2.,tanpsi+self.rangetanpsi/2.)),)

        #if self.DoMerid:
        #    self.workdir=re.sub('/$','_Merid/',self.workdir)

        self.workdir=re.sub('/$','_fixPAinc/',self.workdir)

        print("doing fixed orientation, workdir:",self.workdir)
        self.Run()

        # print("Master.RunFixOrient closing M.fout",self.fout)
        # self.fout.close()

        return
