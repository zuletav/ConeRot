import sys
import os
import os.path
from astropy.io import fits as pf

import scipy.optimize as op
from multiprocessing import Pool
from iminuit import Minuit
import matplotlib.pyplot as plt


include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
from  ConeRot.DConeMaps import *
import ConeRot.KineSummary as KineSummary



import time
#from time import gmtime, strftime
t_i = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())



    
def pass_model(Mpass,OptimMpass):
    global M
    global OptimM
    M=Mpass
    OptimM=OptimMpass
    


def neglnlike4Minuit(*args):
    theta=arange(len(args))
    for iparam in range(len(args)):
        theta[iparam]=args[iparam]
    print( 'theta',theta)
    retval = -1.*lnlike(theta)
    print( "retval",retval)
    return retval

    

def lnlike(theta):

    nvar=len(theta)

    names = list(map( (lambda x: x[0]),M.domain))
    
    for iparam in range(nvar):
        setattr(M,names[iparam],theta[iparam])

    #aPA=theta[0]
    #ainc=theta[1]
    #apsi=theta[2]
    #acosi=np.cos(np.pi*ainc/180.)
    #atanpsi=np.tan(np.pi*apsi/180.)

    chi2=M.conicpolar_expansions()
        
    statusstring=''
        
    for iparam in range(nvar):
        statusstring=statusstring+names[iparam]+" "+str(theta[iparam])+" "

    #statusstring=statusstring+" -> "+str(chi2)+" "+str(M.velodev_med)
    statusstring=statusstring+" -> "+str(chi2)

    if (M.PrintOptimStatus):
        print( statusstring)
    
    return -0.5*chi2



def lnprior(theta):
    inside=1
    bnds = list(map( (lambda x: x[1]),M.domain))
    
    for iparam in list(range(len(theta))):
        if (bnds[iparam][0] < theta[iparam] < bnds[iparam][1]):
            inside *=1
        else:
            inside *=0

    if (inside): 
        return 0.0
    else:
        return -np.inf

    
#def lnprob(theta, bnds):
#    lp = lnprior(theta,bnds)
#    if not np.isfinite(lp):
#        return -np.inf
#    return lp + lnlike(theta)

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)


def run_scipy_optimize_minimize(M,OptimM,x,bnds):
    pass_model(M,OptimM)
    print( "starting op.minimize")
    start_time=time.time()
    nll = lambda *args: -lnlike(*args)
    print( "domain: ",M.domain)
    ftol=0.00001 # 1e-10 too small leads to abnormal termination
    result = op.minimize(nll, x, tol=ftol,bounds=bnds,options={'eps':1E-4})
    print( "result",result)
    result_ml  = result["x"]
    print( "Optim done in (elapsed time):",   time.time()-start_time)
    print( "computing errors with Hessian")
    tmp_i = np.zeros(len(result_ml))
    errors_ml= np.zeros(len(result_ml))
    for i in list(range(len(result_ml))):
        tmp_i[i] = 1.0
        uncertainty_i = np.sqrt(result.hess_inv(tmp_i)[i])
        errors_ml[i]=uncertainty_i
        tmp_i[i] = 0.0
        print(('{0:12.4e} +- {1:.1e}'.format(result.x[i], uncertainty_i)))
    return (result_ml,errors_ml)


def run_Minuit(M,OptimM,x,bnds,names):
    pass_model(M,OptimM)
    print( "starting op.minimize")
    start_time=time.time()

    # f = lambda *args: -lnlike(*args)
    f = lambda *args: -neglnlike4Minuit(*args)

    # a,mu,sigma,a2,mu2,sigma2,base_a,base_b: chi2_2gauss_wbase(selected_velocities, a, mu, sigma, a2, mu2, sigma2, signal_a, rmsnoise, baseparams=[base_a,base_b])

    
    minuitkeyargs={}
    for i in list(range(len(x))):
        minuitkeyargs[names[i]]=x[i]
        minuitkeyargs['limit_'+names[i]]=bnds[i]
        stepsize=abs((bnds[i][1]-bnds[i][0])/100.)
        minuitkeyargs['error_'+names[i]]=stepsize
        

    minuitkeyargs['errordef']=1
    minuitkeyargs['print_level']=0
    minuitkeyargs['pedantic']=0


    print("minuitkeyargs",minuitkeyargs)
    m = Minuit(f, **minuitkeyargs)

    m.migrad()

    print( "Minuit optim done in (elapsed time):",   time.time()-start_time)
    print( "computing errors with Hessian")
    m.hesse() 

    result_ml=np.zeros(len(x))
    errors_ml=np.zeros(len(x))
    for i in list(range(len(x))):
        result_ml[i]=m.values[names[i]]
        errors_ml[i]=m.erros[names[i]]
        print(('{0:12.4e} +- {1:.1e}'.format(result_ml[i], errors_ml[i])))

    return (result_ml,errors_ml)



def exec_ConjGrad_1region(M,OptimM):
    names = list(map( (lambda x: x[0]),M.domain))
    bnds = list(map( (lambda x: x[1]),M.domain))
    nvar=len(list(names))
    sample_theta=list(range(nvar))
    for iparam in list(range(nvar)):
        sample_theta[iparam]=getattr(M,names[iparam])

    x = np.array( sample_theta  ) 

    M.DumpAllFitsFiles=False
    M.Verbose=False
    M.prep_files()
    M.grid_4center()


    if M.DoMinuit:
        (result_ml,errors_ml)=run_Minuit(M,OptimM,x,bnds,names)
    else:
        (result_ml,errors_ml)=run_scipy_optimize_minimize(M,OptimM,x,bnds)
        
    #pass_model(M,OptimM)
    #print( "starting op.minimize")
    #start_time=time.time()
    #nll = lambda *args: -lnlike(*args)
    #print( "domain: ",M.domain)
    ##print( "bnds",bnds)
    #ftol=0.00001 # 1e-10 too small leads to abnormal termination
    #result = op.minimize(nll, x, tol=ftol,bounds=bnds,options={'eps':1E-4})
    #print( "result",result)
    #result_ml  = result["x"]
    #print( "Optim done in (elapsed time):",   time.time()-start_time)
    #print( "computing errors with Hessian")
    #tmp_i = np.zeros(len(result_ml))
    #errors_ml= np.zeros(len(result_ml))
    #for i in list(range(len(result_ml))):
    #    tmp_i[i] = 1.0
    #    #uncertainty_i = np.sqrt(ftol*result.hess_inv(tmp_i)[i])
    #    uncertainty_i = np.sqrt(result.hess_inv(tmp_i)[i])
    #    errors_ml[i]=uncertainty_i
    #    tmp_i[i] = 0.0
    #    print(('{0:12.4e} +- {1:.1e}'.format(result.x[i], uncertainty_i)))
    
    np.save(M.workdir+'result_ml.dat',result_ml)
    np.save(M.workdir+'result_ml_errors.dat',errors_ml)
    return  result_ml


def exec_emcee(M,result_ml,RunMCMC,OptimM):
    Nit=OptimM.Nit
    nwalkers=OptimM.nwalkers
    n_cores=OptimM.n_cores_MCMC
    burn_in=OptimM.burn_in #100


    pass_model(M,OptimM)
    workdir=M.workdir
    names = list(map( (lambda x: x[0]),M.domain))
    bnds = list(map( (lambda x: x[1]),M.domain))
    
    nvar = len(names)
    print( "mcmc with nvar=",nvar)
    
    #ndim =nvar
    ##ndim, nwalkers = nvar, 60
    #pos = [result_ml + 1e-1*np.random.randn(ndim) for i in list(range(nwalkers))]
    
    ranges = list(map( (lambda x: x[1][1]-x[1][0]),M.domain))
        
    allowed_ranges=np.array(ranges)
    print("allowed_ranges ",allowed_ranges)

    
    nvar = len(names)
    print( "mcmc with nvar=",nvar)
    
    ndim =nvar
    #ndim, nwalkers = nvar, 60
    #pos = [result_ml + 1e-1*np.random.randn(ndim) for i in list(range(nwalkers))]
    pos=[]
    for i in list(range(nwalkers)):
        if (np.any(allowed_ranges < 0.)):
            sys.exit("wrong order of bounds in domains")
        awalkerinit=result_ml+(1e-3*np.random.randn(ndim)*allowed_ranges)
        pos.append(awalkerinit)

    print("init for emcee :", result_ml)

    import emcee
    #nit=3000
    print( "in exec_emcee with RunMCMC=",RunMCMC)
    if RunMCMC:
        print( bnds)
        print( "funcs_Optim_DCone:  calling  emcee  with Nit",Nit," nmwalkers",nwalkers," n_cores",n_cores)
        #sampler = emcee.ensemblesampler(nwalkers, ndim, lnprob, args=(bnds))

        #os.environ["OMP_NUM_THREADS"] = "1"
        
        #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=n_cores)
        ##sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
        #sampler.run_mcmc(pos, Nit)



        from multiprocessing import Pool
        with Pool(n_cores) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
            start = time.time()
            sampler.run_mcmc(pos, Nit, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))


        
        print( "************ finish ***************")
        samples = sampler.chain  # chain= array(nwalkers,nit,ndim)
        lnprobs = sampler.lnprobability
        
        ######### save samples
        np.save(workdir+'samples.dat',samples)
        np.save(workdir+'lnprobs.dat',lnprobs)
        # end time
        t_f = time.strftime("%y-%m-%d %h:%m:%s", time.gmtime())
        print( "t_i = "+str(t_i))
        print( "t_f = "+str(t_f))
        
        print(("mean acceptance fraction: {0:.3f} "  .format(np.mean(sampler.acceptance_fraction))))
        f=open(workdir+'acceptance.dat', 'w')
        f.write(str(t_i)+' \n')
        f.write(str(t_f)+' \n')
        f.write("Nit = "+str(Nit)+' \n')
        f.write("nwalkers = "+str(nwalkers)+' \n')
        f.write("ndim = "+str(ndim)+' \n')
        f.write("mean acceptance fraction: {0:.3f}"  .format(np.mean(sampler.acceptance_fraction)) +' \n')
        f.close() 
        
        #autocorr=sampler.get_autocorr_time(c=1, low=1)
        #print( "autocorr\n",autocorr  )
        
    else:
        samples=np.load(workdir+'samples.dat.npy')
        lnprobs=np.load(workdir+'lnprobs.dat.npy')
        



        
    chains=np.zeros(((Nit-burn_in)*nwalkers,ndim))
    chains2=np.zeros((Nit-burn_in, nwalkers,ndim))
    lnpchain=np.zeros(((Nit-burn_in)*nwalkers))
    lnpchain2=np.zeros(((Nit-burn_in), nwalkers))
    


    chains[:,:]=samples[:,burn_in:,:].reshape((nwalkers*(Nit-burn_in), ndim),order='c')
    lnpchain[:]=lnprobs[:,burn_in:].reshape((nwalkers*(Nit-burn_in)),order='c')
    
    ibestparams=np.argmax(lnpchain)
    bestparams=chains[ibestparams,:]
    
    ######### save bestparams
    np.save(workdir+'bestparams.dat',bestparams)
    

    for j in list(range(nwalkers)):
        chains2[:,j,:]=samples[j,burn_in:,:].reshape((Nit-burn_in, ndim),order='c')
        lnpchain2[:,j]=lnprobs[j,burn_in:].reshape(((Nit-burn_in)),order='c')

    fig=plt.figure(figsize=(10,8))
    
    par_labels=names
    ax_lnprob=fig.add_subplot(ndim+1,1,ndim+1)
    for ip in list(range(ndim)):
        ax_chain=fig.add_subplot(ndim+1,1,ip+1)
        for i in list( range(nwalkers)):
            ax_chain.plot(chains2[:,i,ip],alpha=0.1)
            ax_chain.set_ylabel(par_labels[ip])
            
            ax_lnprob.plot(lnpchain2[:,i],alpha=0.1)
            ax_lnprob.set_ylabel('ln(p)')
            
    #plt.show()
    plt.savefig(workdir+'chains.png', bbox_inches='tight')
    plt.close(fig)



    #samples = sampler.chain[:, burn_in:, :].reshape((-1, ndim))

    
    mcmc_results = list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(chains, [16, 50, 84],
                                                axis=0))))


    np.save(workdir+'mcmc_results.dat',mcmc_results)

    mcmc_results_0 = np.zeros(nvar)

    print( "param     distrib     max ")
    for iparam in list(range(nvar)):
        print( names[iparam],mcmc_results[iparam],bestparams[iparam])
        mcmc_results_0[iparam]= mcmc_results[iparam][0]
        

    #print( "mcmc median values:")
    #model_median =  np.array(modelfunk(mcmc_results_0, m))


    
    import corner

    fig=corner.corner(chains,
                      labels=names,
                      quantiles=[0.16, 0.5,0.84],
                      bins=20, truths=bestparams,
                      levels=[0.68, 0.95, 0.997],
                      show_titles=True,
                      title_fmt=".3f",
                      title_kwards={"fontsize": 10}) #, smooth=1.0




    fig.savefig(workdir+M.TriangleFile)

    print( "finished MCMC for region workdir",workdir)
    return [names,mcmc_results]




######################################################################

#Regions
def proc_1region(iregion):
    names = list(map( (lambda x: x[0]),M.domain))
    bnds = list(map( (lambda x: x[1]),M.domain))
    nvar=len(names)
    sample_theta=list(range(nvar))
    for iparam in list(range(nvar)):
        sample_theta[iparam]=getattr(M,names[iparam])
    
        
    x = np.array( sample_theta  ) 

    amesh=M.a_min_regions+ np.arange(M.n_abins+1)*(M.a_max_regions - M.a_min_regions)/M.n_abins
    ameshbis=np.roll(amesh, -1)

    #M.filelog='log_output_region'+str(iregion)+'.txt'
    #print( "opening log for region"+str(iregion)+":",M.workdir+M.filelog)

    M.a_min=amesh[iregion]
    M.a_max=ameshbis[iregion+1]
    print( ">>>>> iregion ",iregion," from ", M.a_min," to ", M.a_max)
    logstring=">>>>> "+str(iregion)+" from %.3f to %.3f \n" % (M.a_min,M.a_max)
    
    M.DumpAllFitsFiles=False
    M.Verbose=False

    masterworkdir=M.workdir

    workdir_region=masterworkdir+'work_region_'+str(iregion)+'/'


    os.system("rm -rf  "+workdir_region)

    os.system("mkdir "+workdir_region)
    M.workdir=workdir_region
    M.ComputeSkyImages=False


    # pass_model(M,OptimM)
    if M.DoConjGrad:
        (result_ml,errors_ml)=run_scipy_optimize_minimize(M,OptimM,x,bnds)


        if (M.StoreRegions):
            np.save(workdir_region+'result_ml_region'+str(iregion)+'.dat',result_ml)
        np.save(masterworkdir+'result_ml_region'+str(iregion)+'.dat',result_ml)
        np.save(M.workdir+'result_ml.dat',result_ml)
        np.save(masterworkdir+'result_ml_errors_region'+str(iregion)+'.dat',errors_ml)
        np.save(M.workdir+'result_ml_errors.dat',errors_ml)

    result_ml=np.load(masterworkdir+'result_ml_region'+str(iregion)+'.dat.npy')
    errors_ml=np.load(masterworkdir+'result_ml_errors_region'+str(iregion)+'.dat.npy')
    print( "result_ml_region is ",result_ml)
    for iparam in list(range(nvar)):
        print( names[iparam],"->",result_ml[iparam])
        setattr(M,names[iparam],result_ml[iparam])
        logstring=logstring+names[iparam]+"-> %.6f +- %.7f " % (result_ml[iparam],errors_ml[iparam])
    logstring=logstring+"\n"


    if (M.RunMCMC):
        M.TriangleFile='triangle_'+str(iregion)+'.png'
        #OptimM=Optim_DCone.OptimModel(M,RunMCMC=True,Nit=Nit,nwalkers=nwalkers,n_cores_MCMC=n_cores_MCMC)
        print( "running emcee for region "+str(iregion))
        print( "M.RunMCMC",M.RunMCMC)
        retvals = OptimM.emcee(M)
        print( "returned from OptimM.emcee")
        names=retvals[0]
        mcmc_results=retvals[1]
        # [names, mcmc_results] = 
        print( "looping over MCMC optim")
        logstring+="emcee posterior\n"
        for iparam in list(range(nvar)):
            strparams=names[iparam]+" -> %.6f %.6f %.6f " % mcmc_results[iparam]
            logstring=logstring+strparams+"\n"

        # OptimM.RecoverMCMC(M)


    
    if (M.StoreRegions):
        M.DumpAllFitsFiles=True

    M.prep_files()
    M.grid_4center()
    M.ComputeSkyImages=True
    chi2=M.conicpolar_expansions()
    print( "chi2=",chi2)

    if (M.StoreRegions):

        inbasename=os.path.basename(M.filename_source)
        inbasename=re.sub('.fits', '', inbasename)
        inbasename=workdir_region+inbasename
        fileout = inbasename+'_fig_summary_region'+str(iregion)+'.pdf'
        exec_summary(inbasename,fileout)

    M.workdir=masterworkdir
    
    # if M.DoDCone:
    #   # return [iregion,M.Hduregion.data,M.Hdudiff.data,M.Hdumoddrot.data,logstring,M.Hdumumap.data,M.HduDConemoddrot.data,M.HdudiffDConemoddrot.data,M.RadialProfile,M.a_min,M.a_max]
    #    return [iregion,M,logstring]
    #else:
    #    # return [iregion,M.Hduregion.data,M.Hdudiff.data,M.Hdumoddrot.data,logstring,M.RadialProfile,M.a_min,M.a_max]
    #    return [iregion,M,logstring]

    

    print( "Done processing region ",iregion)

    if M.DoDCone:
        return [iregion,logstring,M.Hduregion.data,M.Hdudiff.data,M.Hdumoddrot.data,M.RadialProfile,M.a_min,M.a_max,M.Hdudiff_faceon.data,M.Hduresamp_faceon.data,M.Hduregion_faceon.data,M.Hdumumap.data,M.HduDConemoddrot.data,M.HdudiffDConemoddrot.data]
    else:
        return [iregion,logstring,M.Hduregion.data,M.Hdudiff.data,M.Hdumoddrot.data,M.RadialProfile,M.a_min,M.a_max,M.Hdudiff_faceon.data,M.Hduresamp_faceon.data,M.Hduregion_faceon.data]



def exec_Regions(M,OptimM):
    n_cores_regions=OptimM.n_cores_regions
    print( "n_cores_regions = ",n_cores_regions)
    amesh=M.a_min_regions+ np.arange(M.n_abins+1)*(M.a_max_regions - M.a_min_regions)/M.n_abins
    ameshbis=np.roll(amesh, -1)
    print( "amesh ",amesh)
    print( "ameshbis ",ameshbis)
    
    
    im_c=M.Hducentered.data
    hdr_c=M.Hducentered.header
    (ny,nx)=im_c.shape
    #print( "master shape:",(ny,nx))
    cube_regions=np.zeros((M.n_abins-1,ny,nx))
    cube_imdrotdiff=np.zeros((M.n_abins-1,ny,nx))
    cube_immoddrot=np.zeros((M.n_abins-1,ny,nx))

    cube_im_diff_faceon=np.zeros((M.n_abins-1,ny,nx))
    cube_resamp_faceon=np.zeros((M.n_abins-1,ny,nx))
    cube_regions_faceon=np.zeros((M.n_abins-1,ny,nx))

    if M.DoDCone:
        cube_mumap=np.zeros((M.n_abins-1,ny,nx))
        cube_imDConemoddrot=np.zeros((M.n_abins-1,ny,nx))
        cube_diffimDConemoddrot=np.zeros((M.n_abins-1,ny,nx))

    print("Regions M.fout",M.fout)
    M.fout.write("Regions:\n")
    
    M.PrintOptimStatus=False
    workdir=M.workdir

    if M.RunMCMC:
        n_cores_regions=1

    pass_model(M,OptimM)

    if (n_cores_regions > 1):
        p = Pool(n_cores_regions)
        passoutput = p.map(proc_1region, range(M.n_abins-1))
    else:
        passoutput=[]
        for iregion in list(range(M.n_abins-1)):
            passoutput.append(proc_1region(iregion))
            

    # rrs0=passoutput[0][-3][0]
    rrs0=passoutput[0][5][0]
    nrs=len(rrs0)
    stack_v_Phi_profiles=np.zeros((len(passoutput),nrs))
    stack_v_R_profiles=np.zeros((len(passoutput),nrs))
    stack_sv_R_profiles=np.zeros((len(passoutput),nrs))
    stack_v_z_profiles=np.zeros((len(passoutput),nrs))
    stack_sv_z_profiles=np.zeros((len(passoutput),nrs))
    stack_sprofiles=np.zeros((len(passoutput),nrs))
    stack_vecregions=np.zeros((len(passoutput),nrs))
    
    for aregionoutput in passoutput:
        iregion=aregionoutput[0]

        logstring=aregionoutput[1]
        cube_regions[iregion,:,:] = aregionoutput[2]
        cube_imdrotdiff[iregion,:,:] = aregionoutput[3]
        cube_immoddrot[iregion,:,:] = aregionoutput[4]
        radialprofile = aregionoutput[5]
        amin = aregionoutput[6]
        amax = aregionoutput[7]
        cube_im_diff_faceon[iregion,:,:] = aregionoutput[8]
        cube_resamp_faceon[iregion,:,:] = aregionoutput[9]
        cube_regions_faceon[iregion,:,:] = aregionoutput[10]

        M.fout.write(logstring)

        if M.DoDCone:
            cube_mumap[iregion,:,:] = aregionoutput[11]
            cube_imDConemoddrot[iregion,:,:] = aregionoutput[12]
            cube_diffimDConemoddrot[iregion,:,:] = aregionoutput[13]

        rrs=radialprofile[0]
        v_Phi_prof=radialprofile[1]
        sv_Phi_prof=radialprofile[2]
        if (M.DoMerid):
            v_R_prof=radialprofile[3]
            sv_R_prof=radialprofile[4]
            v_z_prof=radialprofile[5]
            sv_z_prof=radialprofile[6]
        elif (M.DoAccr):
            v_R_prof=radialprofile[3]
            sv_R_prof=radialprofile[4]
            
        vecregion=np.zeros(len(rrs))
        vecregion[np.where( (rrs  >= amin) & (rrs <= amax))] = 1.
        stack_v_Phi_profiles[iregion,:] = v_Phi_prof
        stack_sprofiles[iregion,:] = sv_Phi_prof
        stack_vecregions[iregion,:] = vecregion

        if (M.DoMerid):
            stack_v_z_profiles[iregion,:] = v_z_prof
            stack_sv_z_profiles[iregion,:] = sv_z_prof
            stack_v_R_profiles[iregion,:] = v_R_prof
            stack_sv_R_profiles[iregion,:] = sv_R_prof
        elif (M.DoAccr):
            stack_v_R_profiles[iregion,:] = v_R_prof
            stack_sv_R_profiles[iregion,:] = sv_R_prof
            


    print( "Collapsing regions")


    vec_norm=np.sum(stack_vecregions,axis=0)
    mask=(vec_norm < 0.1)
    
    
    allrads_v_Phi_prof=np.sum(stack_v_Phi_profiles*stack_vecregions,axis=0)/vec_norm
    allrads_v_Phi_prof[mask] =0.
    allrads_sv_Phi_prof=np.sum(stack_sprofiles*stack_vecregions,axis=0)/vec_norm
    allrads_sv_Phi_prof[mask] =0.

    if (M.DoMerid):
        allrads_v_R_prof=np.sum(stack_v_R_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_v_R_prof[mask] =0.
        allrads_sv_R_prof=np.sum(stack_sv_R_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_sv_R_prof[mask] =0.   
        allrads_v_z_prof=np.sum(stack_v_z_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_v_z_prof[mask] =0.
        allrads_sv_z_prof=np.sum(stack_sv_z_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_sv_z_prof[mask] =0.   
    elif (M.DoAccr):
        allrads_v_R_prof=np.sum(stack_v_R_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_v_R_prof[mask] =0.
        allrads_sv_R_prof=np.sum(stack_sv_R_profiles*stack_vecregions,axis=0)/vec_norm
        allrads_sv_R_prof[mask] =0.


    inbasename=os.path.basename(M.filename_source)
    filename_fullim=re.sub('.fits', '_fullim.fits', inbasename)
    filename_fullim=workdir+filename_fullim
    fileout_cubediff=re.sub('fullim.fits', 'cube_azim_av_drot_diff.fits', filename_fullim)
    fileout_cuberegions=re.sub('fullim.fits', 'cube_regions.fits', filename_fullim)
    
    pf.writeto(fileout_cubediff,cube_imdrotdiff, hdr_c, overwrite=True)
    pf.writeto(fileout_cuberegions,cube_regions, hdr_c, overwrite=True)
    
    im_norm = np.sum(cube_regions, axis = 0)
    mask= (np.fabs(im_norm) < 0.01)
    im_norm[mask] = 0.01
    
    imdrotdiff= np.sum(cube_imdrotdiff * cube_regions, axis = 0) / im_norm
    imdrotdiff[mask]=0.

    immoddrot= np.sum(cube_immoddrot * cube_regions, axis = 0) / im_norm
    immoddrot[mask]=0.

    imdrotdiff_b=im_c-immoddrot
    imdrotdiff_b[mask]=0.

    im_norm_faceon = np.sum(cube_regions_faceon, axis = 0)
    mask_fon= (np.fabs(im_norm_faceon) < 0.01)
    im_norm_faceon[mask_fon] = 0.01

    imdiff_faceon= np.sum(cube_im_diff_faceon * cube_regions_faceon, axis = 0) / im_norm_faceon
    imdiff_faceon[mask_fon]=0.

    imresamp_faceon= np.sum(cube_resamp_faceon * cube_regions_faceon, axis = 0) / im_norm_faceon
    imresamp_faceon[mask_fon]=0.




    if M.DoDCone:
        imDConemoddrot= np.sum(cube_imDConemoddrot * cube_regions, axis = 0) / im_norm
        imDConemoddrot[mask]=0.
        imdiffDConemoddrot= np.sum(cube_diffimDConemoddrot * cube_regions, axis = 0) / im_norm
        imdiffDConemoddrot[mask]=0.
        immumap= np.sum(cube_mumap * cube_regions, axis = 0) / im_norm
        
    

    fileout_diff=re.sub('fullim.fits', 'allrads_azim_av_drot_diff.fits', filename_fullim)
    pf.writeto(fileout_diff,imdrotdiff, hdr_c, overwrite=True)
    fileout_diff_b=re.sub('fullim.fits', 'allrads_azim_av_drot_diff_b.fits', filename_fullim)
    pf.writeto(fileout_diff_b,imdrotdiff_b, hdr_c, overwrite=True)
    fileout_imnorm=re.sub('fullim.fits', 'imregions.fits', filename_fullim)
    pf.writeto(fileout_imnorm,im_norm, hdr_c, overwrite=True)
    fileout_immoddrot=re.sub('fullim.fits', 'allrads_azim_av_drot.fits', filename_fullim)
    pf.writeto(fileout_immoddrot,immoddrot, hdr_c, overwrite=True)


    fileout_diff_faceon=re.sub('fullim.fits', 'allrads_diff_faceon.fits', filename_fullim)
    pf.writeto(fileout_diff_faceon,imdiff_faceon, hdr_c, overwrite=True)
    fileout_resamp_faceon=re.sub('fullim.fits', 'allrads_resamp_faceon.fits', filename_fullim)
    pf.writeto(fileout_resamp_faceon,imresamp_faceon, hdr_c, overwrite=True)
    fileout_imnorm_faceon=re.sub('fullim.fits', 'imregions_faceon.fits', filename_fullim)
    pf.writeto(fileout_imnorm_faceon,im_norm_faceon, hdr_c, overwrite=True)

    im_c_w=M.Hduwcentered.data # _c -> centered
    regionsmask=np.where(im_norm_faceon > 0.9)
    chi2regions=np.sum( im_c_w[regionsmask]*(imdiff_faceon[regionsmask])**2)/M.Ncorr
    M.fout.write("chi2regions=%.6e\n" % (chi2regions))

    if (M.DoMerid):
        save_prof = np.zeros((nrs,7))
        save_prof[:,0] = rrs0
        save_prof[:,1] = allrads_v_Phi_prof
        save_prof[:,2] = allrads_sv_Phi_prof
        save_prof[:,3] = allrads_v_R_prof
        save_prof[:,4] = allrads_sv_R_prof
        save_prof[:,5] = allrads_v_z_prof
        save_prof[:,6] = allrads_sv_z_prof
    elif (M.DoAccr):
        save_prof = np.zeros((nrs,5))
        save_prof[:,0] = rrs0
        save_prof[:,1] = allrads_v_Phi_prof
        save_prof[:,2] = allrads_sv_Phi_prof
        save_prof[:,3] = allrads_v_R_prof
        save_prof[:,4] = allrads_sv_R_prof
    else:
        save_prof = np.zeros((nrs,3))
        save_prof[:,0] = rrs0
        save_prof[:,1] = allrads_v_Phi_prof
        save_prof[:,2] = allrads_sv_Phi_prof
        
    fileout_allradsradialprofile=re.sub('fullim.fits', 'allrads_radial_profile.dat', filename_fullim)
    np.savetxt(fileout_allradsradialprofile, save_prof)   # x,y,z equal sized 1D arrays


    if M.DoDCone:
        fileout_imDConemoddrot=re.sub('fullim.fits', 'allrads_immod_DCone.fits', filename_fullim)
        pf.writeto(fileout_imDConemoddrot,imDConemoddrot, hdr_c, overwrite=True)

        fileout_imdiffDConemoddrot=re.sub('fullim.fits', 'allrads_diff_DCone.fits', filename_fullim)
        pf.writeto(fileout_imdiffDConemoddrot,imdiffDConemoddrot, hdr_c, overwrite=True)

        fileout_immumap=re.sub('fullim.fits', 'allrads_mumap.fits', filename_fullim)
        pf.writeto(fileout_immumap,immumap, hdr_c, overwrite=True)


    inbasename=os.path.basename(M.filename_source)
    inbasename=re.sub('.fits', '', inbasename)
    inbasename=workdir+inbasename
    fileout = inbasename+'_allrads_fig_summary.pdf'
    
    KineSummary.exec_summary_allrads(inbasename,fileout,vsyst=M.vsyst)
    

    return
    
