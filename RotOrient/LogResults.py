import numpy as np
import re

def load(workdir,FixPAInc=False,RunMCMC=False):
    
    file_log=workdir+'log_output.txt'

    print( "loading file_log",file_log,"FixPAInc",FixPAInc)
    fin= open(file_log,"r")
    log_output=fin.readlines()
    fin.close
    #print( log_output)
    #print( "entries:",len(log_output))

    AllRads=False
    AllRadsMCMC=False
    Regions=False
    emcee_posterior=False
    
    
    
    if FixPAInc:
        r1s=[]
        r2s=[]

        rregions_fixincPA=[]
        tanpsis_fixincPA=[]
        tanpsis_sigma_fixincPA=[]

        tanpsis_fixincPA_MCMC=[]
        for aline in log_output:
            if AllRads:
                #(allPA, allinc, alltanpsi) = re.search("^PA-> (.*) inc-> (.*) tanpsi-> (.*)$",aline)
                matches = re.search("^tanpsi-> (.*) $",aline)
                #instring=matches.group(0)
                allradstanpsi_fixincPA=float(matches.group(1))
                AllRads=False
            if "Global" in aline:
                AllRads=True
            if (emcee_posterior and (not Regions)):
                matches = re.search("^tanpsi ->\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                if (matches):
                    atanpsi_post=matches.group(1)
                    atanpsi_upsigma=matches.group(2)
                    atanpsi_downsigma=matches.group(3)
                    tanpsi_fixincPA_MCMC=[float(atanpsi_post),float(atanpsi_upsigma),float(atanpsi_downsigma)]
                    emcee_posterior=False
            if ("emcee posterior" in aline):
                emcee_posterior=True
                
            if Regions: 
                matches = re.search("^>>>>> .* from (.*) to (.*)$",aline)
                if (matches):
                    ar1=matches.group(1)
                    ar2=matches.group(2)
                    #print( "a region radiii",ar1,ar2)
                    rregions_fixincPA.append((float(ar1)+float(ar2))/2.)
                    r1s.append(float(ar1))
                    r2s.append(float(ar2))
                matches = re.search("^tanpsi->\s*(\-?\d+\.\d+)\s*\+\-\s*(\d+\.\d+)\s*$",aline)
                if (matches):
                    atanpsi=matches.group(1)        
                    atanpsi_sigma=matches.group(2)        
                    tanpsis_fixincPA.append(float(atanpsi))
                    tanpsis_sigma_fixincPA.append(float(atanpsi_sigma))
                if emcee_posterior:
                    matches = re.search("^tanpsi ->\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    if (matches):
                        atanpsi_post=matches.group(1)
                        atanpsi_upsigma=matches.group(2)
                        atanpsi_downsigma=matches.group(3)
                        tanpsis_fixincPA_MCMC.append([float(atanpsi_post),float(atanpsi_upsigma),float(atanpsi_downsigma)])
                        emcee_posterior=False
                if ("emcee posterior" in aline):
                    emcee_posterior=True

            if "Regions" in aline:
                Regions=True

        rregions_fixincPA=np.array(rregions_fixincPA)
        r1s=np.array(r1s)
        r2s=np.array(r2s)
        tanpsis_fixincPA=np.array(tanpsis_fixincPA)

        psis_fixincPA=np.arctan(tanpsis_fixincPA)*180./np.pi


        allradspsi_fixincPA= np.arctan( tanpsi_fixincPA_MCMC[0]) * 180. / np.pi

        if RunMCMC:

            print(">>>fix>>psis_fixincPA ", psis_fixincPA)
            print(">>>fix>>allradspsi_fixincPA", allradspsi_fixincPA,"tanpsi_fixincPA_MCMC[0]",tanpsi_fixincPA_MCMC[0])
            print(">>>fix>>tanpsis_fixincPA_MCMC",tanpsis_fixincPA_MCMC)
            print(">>>fix>>psis_fixincPA_MCMC",180.*(np.arctan(tanpsis_fixincPA_MCMC))/np.pi)
            
            return (RunMCMC, [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, tanpsis_fixincPA_MCMC])
        else:
            return (RunMCMC, [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA])
            
    else:

        r1s=[]
        r2s=[]
        rregions=[]
        PAs=[]
        PAs_sigma=[]
        incs=[]
        incs_sigma=[]
        tanpsis=[]
        tanpsis_sigma=[]
        PAs_MCMC=[]
        incs_MCMC=[]
        tanpsis_MCMC=[]

    
        
        for aline in log_output:
            if AllRads:
                matches = re.search("^PA-> (.*) inc-> (.*) tanpsi-> (.*) $",aline)
                allradsPA=float(matches.group(1))
                allradsinc=float(matches.group(2))*180./np.pi
                allradstanpsi=float(matches.group(3))
                AllRads=False
            if AllRadsMCMC:
                matches = re.search("^PA-> (.*) inc-> (.*) tanpsi-> (.*) $",aline)
                allradsPAMCMC=float(matches.group(1))
                allradsincMCMC=float(matches.group(2))*180./np.pi
                allradstanpsiMCMC=float(matches.group(3))
                AllRadsMCMC=False
            if "Global" in aline:
                AllRads=True
            if "Global MCMC" in aline:
                AllRadsMCMC=True
                RunMCMC=True

            if (emcee_posterior and (not Regions)):
                matches = re.search("^PA ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                if (matches):
                    aPA_post=matches.group(1)
                    aPA_upsigma=matches.group(2)
                    aPA_downsigma=matches.group(3)
                    PA_MCMC=[float(aPA_post),float(aPA_upsigma),float(aPA_downsigma)]
                matches = re.search("^inc ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                if (matches):
                    ainc_post=matches.group(1)
                    ainc_upsigma=matches.group(2)
                    ainc_downsigma=matches.group(3)
                    inc_MCMC=[float(ainc_post),float(ainc_upsigma),float(ainc_downsigma)]
                matches = re.search("^tanpsi ->\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                if (matches):
                    atanpsi_post=matches.group(1)
                    atanpsi_upsigma=matches.group(2)
                    atanpsi_downsigma=matches.group(3)
                    tanpsi_MCMC=[float(atanpsi_post),float(atanpsi_upsigma),float(atanpsi_downsigma)]
                    emcee_posterior=False

            if ("emcee posterior" in aline):
                emcee_posterior=True
                
            if Regions:
                matches = re.search("^>>>>> .* from (.*) to (.*)$",aline)
                if (matches):
                    ar1=matches.group(1)
                    ar2=matches.group(2)
                    rregions.append((float(ar1)+float(ar2))/2.)
                    r1s.append(float(ar1))
                    r2s.append(float(ar2))

                matches = re.search("^PA->\s*(\d+\.\d+)\s*\+\-\s*(\d+\.\d+)\s*inc->\s*(\-?\d+\.\d+)\s*\+\-\s*(\d+\.\d+)\s*tanpsi->\s*(\-?\d+\.\d+)\s*\+\-\s*(\d+\.\d+)\s*$",aline)
                if (matches):
                    aPA=matches.group(1)
                    aPA_sigma=matches.group(2)
                    ainc=matches.group(3)
                    ainc_sigma=matches.group(4)
                    atanpsi=matches.group(5)        
                    atanpsi_sigma=matches.group(6)
                    
                    
                    PAs.append(float(aPA))
                    incs.append(float(ainc))
                    tanpsis.append(float(atanpsi))
                    PAs_sigma.append(float(aPA_sigma))
                    incs_sigma.append(float(ainc_sigma))
                    tanpsis_sigma.append(float(atanpsi_sigma))
                    
                if emcee_posterior:
                    matches = re.search("^PA ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    #matches = re.search("^PA ->\s+(\w+)\s+(\w+)\s+(\w+)\s*$",aline)
                    if (matches):
                        aPA_post=matches.group(1)
                        aPA_upsigma=matches.group(2)
                        aPA_downsigma=matches.group(3)
                        PAs_MCMC.append([float(aPA_post),float(aPA_upsigma),float(aPA_downsigma)]) 
                    matches = re.search("^inc ->\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    if (matches):
                        ainc_post=matches.group(1)
                        ainc_upsigma=matches.group(2)
                        ainc_downsigma=matches.group(3)
                        incs_MCMC.append([float(ainc_post),float(ainc_upsigma),float(ainc_downsigma)])
                    matches = re.search("^tanpsi ->\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$",aline)
                    if (matches):
                        atanpsi_post=matches.group(1)
                        atanpsi_upsigma=matches.group(2)
                        atanpsi_downsigma=matches.group(3)
                        print("found atanpsi_post ",atanpsi_post)
                        tanpsis_MCMC.append([float(atanpsi_post),float(atanpsi_upsigma),float(atanpsi_downsigma)])
                        emcee_posterior=False
                if ("emcee posterior" in aline):
                    emcee_posterior=True

            if "Regions" in aline:
                Regions=True


        rregions=np.array(rregions)
        r1s=np.array(r1s)
        r2s=np.array(r2s)
        incs=180.*np.array(incs)/np.pi
        tanpsis=np.array(tanpsis)
        PAs=np.array(PAs)
        psis=np.arctan(tanpsis)*180./np.pi
        
        allradspsi= np.arctan( allradstanpsi) * 180. / np.pi


                
        if RunMCMC:

            allradspsiMCMC= np.arctan( tanpsi_MCMC[0]) * 180. / np.pi
            allradsincMCMC= inc_MCMC[0]  * 180. / np.pi
            allradsPAMCMC= PA_MCMC[0]
            
            return (RunMCMC, [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi, PAs_MCMC, tanpsis_MCMC, incs_MCMC, allradsPAMCMC, allradsincMCMC,allradspsiMCMC])


            
            #return (RunMCMC, [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi, PAs_MCMC, tanpsis_MCMC, incs_MCMC, PA_MCMC[0], inc_MCMC[0],tanpsi_MCMC[0]])
        else:
            return (RunMCMC, [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi] )

