import numpy as np
import time
import MB21
import logg_SB
def MC_logg_Teff(C,met,logg,logg_err, dist,dist_err,G,ag,iter,*args,**kwargs):

   #First iteration k=0 with logg_sb as an input:
    #Teff_it, Teff_err_it =  MB21.MB21(C,met,logg,logg_err,100)
    #logg_it, logg_err_it = logg_SB.logg_SB(Teff_it, Teff_err_it, dist, dist_err, G, ag, 100)
    logg_it=logg
    logg_err_it=logg_err

    #Building iterative process for k in range (1,iter)
    t0 = time.time()
    k=1
    while k<iter:
        Teff_it, Teff_err_it = MB21.MB21(C,met,logg_it,logg_err_it,100)
        logg_it, logg_err_it = logg_SB.logg_SB(Teff_it, Teff_err_it, dist, dist_err, G, ag, 100)
        k=k+1

    t1 = time.time()
    print('Elapsed time: %s s' % round(t1-t0))
    print('Teff_final:')
    for i in range (len(Teff_it)):
        print(Teff_it[i], '±', Teff_err_it[i])
    print('logg_final: ')
    for i in range (len(Teff_it)):
        print(logg_it[i], '±', logg_err_it[i])

    return Teff_it,logg_it,Teff_err_it,logg_err_it
#diff_T,diff_logg
