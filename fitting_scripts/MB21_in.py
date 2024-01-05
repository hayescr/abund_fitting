import numpy as np
def MB21_in(C,met,iter,*args,**kwargs):

    b0= 0.4929
    b1= 0.5092
    b2= -0.0353
    b3= 0.0192
    b4= -0.002
    b5= -0.0395
    B0= 0.5323
    B1= 0.4775
    B2= -0.0344
    B3= -0.011
    B4= -0.002
    B5= -0.0009

    C_r = np.empty((len(C),iter))
    met_r = np.empty((len(C),iter))
    theta_d_r = np.empty((len(C),iter))
    theta_g_r = np.empty((len(C),iter))
    Teff_d_r = np.empty((len(C),iter))
    Teff_g_r = np.empty((len(C),iter))
    Teff_d = np.empty((len(C),1))
    Teff_g = np.empty((len(C),1))
    Teff_mb = np.empty((len(C),1))
    Teff_mb_err = np.empty((len(C),1))

    for i in range(len(C)):

        for k in range(iter):

            C_r[i,k] = np.random.normal(C[i],0.01)
            #met_r[k]=met[i] + (np.random.rand(1) - 0.5) / 0.5*0.5
            met_r[i, k]=met[i] + (np.random.rand(1) - 0.5)
            #met_r = np.append(met_r,met_r[k])
            theta_d_r[i,k]=b0 + b1*C_r[i,k] + b2*C_r[i,k] ** 2 + b3*met_r[i,k] + b4*met_r[i,k] ** 2 + b5*met_r[i,k]*C_r[i,k]

            Teff_d_r[i, k]=5040 / theta_d_r[i,k]
            theta_g_r[i, k]=B0 + B1*C_r[i,k] + B2*C_r[i,k] ** 2 + B3*met_r[i,k] + B4*met_r[i,k] ** 2 + B5*met_r[i,k]*C_r[i,k]

            Teff_g_r[i, k]=5040 / theta_g_r[i,k]

        Teff_d[i]=np.median(Teff_d_r[i])
        #Teff_d_err(i)=std(Teff_d_r)
        Teff_g[i]=np.median(Teff_g_r[i])
        #Teff_g_err(i)=std(Teff_g_r)

        Teff_mb[i]=np.mean([Teff_g[i],Teff_d[i]])
        Teff_mb_err[i]=100

    return Teff_mb, Teff_mb_err
