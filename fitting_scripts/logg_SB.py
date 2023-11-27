import numpy as np
import math
def logg_SB(Teff,Teff_err,dist,dist_err,G,ag,iter,*args,**kwargs):

    Tsun=5778
    a0=0.06
    a1=6.731e-05
    a2=- 6.647e-08
    a3=2.859e-11
    a4=- 7.197e-15

    Teff_r = np.empty((len(Teff),iter))
    mass = np.empty((len(Teff),iter))
    BC = np.empty((len(Teff),iter))
    dist_r = np.empty((len(Teff),iter))
    G_abs_r = np.empty((len(Teff),iter))
    lum_log_r = np.empty((len(Teff),iter))
    logg_sb_r = np.empty((len(Teff),iter))
    logg_sb = np.empty((len(Teff),1))
    logg_sb_err = np.empty((len(Teff),1))

    for i in range(len(Teff)):

        for k in range(iter):

           #Teff_r[i,k]=np.random.normal(Teff(=[i],Teff_err[i])
            Teff_r[i,k] = np.random.normal( Teff[i],Teff_err[i])
            mass[i,k] = 0.3*np.random.rand(1,1)+ 0.65;
            BC[i,k]=a0 + a1*(Teff_r[i,k] - Tsun) + a2*(Teff_r[i,k] - Tsun) ** 2 + a3*(Teff_r[i,k] - Tsun) ** 3 + a4*(Teff_r[i,k] - Tsun) ** 4
            #dist_r[i,k]=np.random.normal(dist[i],dist_err[i])
            dist_r[i,k]= np.absolute(np.random.normal( dist[i], dist_err[i] ))
            G_abs_r[i,k]=G[i] + 5 - 5*math.log10(dist_r[i,k]*1000) - ag[i]
            lum_log_r[i,k]=(G_abs_r[i,k] + BC[i,k] - 4.74) / (- 2.5)
            logg_sb_r[i,k]=math.log10(mass[i,k]) - (lum_log_r[i,k] - 4*math.log10(Teff_r[i,k] / Tsun)) + 4.437

        logg_sb[i]=np.median(logg_sb_r[i,:])   
        logg_sb_err[i]=np.std(logg_sb_r[i,:])

    return logg_sb,logg_sb_err
