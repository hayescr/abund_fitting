import numpy as np
def MB21(C,met,logg,logg_err,iter,*args, **kwargs):

    logg_r = np.empty((len(met),iter))
    C_r = np.empty((len(met),iter))
    met_r = np.empty((len(met),iter))
    theta_r = np.empty((len(met),iter))
    Teff_mb_r = np.empty((len(met),iter))
    Teff_mb = np.empty((len(met),1))
    Teff_mb_err = np.empty((len(met),1))

    for i in range(len(met)):

        for k in range(iter):

            logg_r[i,k]=np.random.normal(logg[i],logg_err[i])

            if logg_r[i,k] >= 4:
                #For Dwarfs ======
                b0=0.4929
                b1=0.5092
                b2=- 0.0353
                b3=0.0192
                b4=- 0.002
                b5=- 0.0395

            else:
                #For Giants =======
                b0=0.5323
                b1=0.4775
                b2=- 0.0344
                b3=- 0.011
                b4=- 0.002
                b5=- 0.0009


            C_r[i,k]=np.random.normal(C[i],0.05) #C(i)+(np.random.rand(1)-0.5)/0.5*0.1

            met_r[i,k]=met[i] + (np.random.rand(1) - 0.5)
            theta_r[i,k]=b0 + b1*C_r[i,k] + b2*C_r[i,k] ** 2 + b3*met_r[i,k] + b4*met_r[i,k] ** 2 + b5*met_r[i,k]*C_r[i,k]
            Teff_mb_r[i,k]=5040 / theta_r[i,k]

        Teff_mb[i]=np.median(Teff_mb_r[i,:])
        Teff_mb_err[i]=np.std(Teff_mb_r[i,:])


        #Teff_mb(i)=5040/theta(i);

    return Teff_mb,Teff_mb_err
