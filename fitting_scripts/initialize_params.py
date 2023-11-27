import numpy as np
from analysis_utils import read_json, write_json
import os

sp_file = np.genfromtxt('stellar_params_init.txt', dtype=None, names=True, encoding=None)

for star_name, teff, logg, feh, vmicro, alphafe, cfe, nfe, vmacro in zip(sp_file['STAR'], 
                                                                         sp_file['TEFF'], 
                                                                         sp_file['LOGG'], 
                                                                         sp_file['FEH'], 
                                                                         sp_file['VMICRO'], 
                                                                         sp_file['ALPHA_FE'], 
                                                                         sp_file['C_FE'], 
                                                                         sp_file['N_FE'], 
                                                                         sp_file['VMACRO']):


    filename = f'{star_name}/{star_name}_params.txt'

    if not os.path.isdir(f'{star_name}'):
        os.mkdir(f'{star_name}')
    
    param_dict = {'teff' : teff,
                  'logg' : logg,
                  'feh' : feh,
                  'vmicro' : vmicro,
                  'cfe' : cfe,
                  'nfe' : nfe,
                  'alphafe' : alphafe,
                  'vmacro' : vmacro,
                  'abundances' : {},
                  'isotopes' : {},
                  }

    write_json(filename, param_dict)