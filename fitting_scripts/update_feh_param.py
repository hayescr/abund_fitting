import numpy as np
from analysis_utils import read_json, write_json
from scipy.stats import sigmaclip
import abund_utils as au


star_list = ['ret2_mem', 'ret2_cand']
element = 'Fe'

solar_abundances = au.solar_abund('A09')

for star_name in star_list:
    star = star_name
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)
    res = star_param['vmacro']

    abund_filename = f'{star}/{star}_{element}_abund_fit_convol{res:.1f}.txt'
    abund_file = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)

    good_lines = (abund_file['lim_flag_eqw'].astype('bool')) & (abund_file['lim_flag_chi2'].astype('bool'))
    good_lines &= abund_file['avg_synth_stdev']/(3/(abund_file['snr']*np.sqrt(abund_file['npix']))) > 1.0
    good_lines &= (abund_file['avg_line_depth'] * abund_file['snr'] > 5.)

    a, clip_low, clip_high = sigmaclip(abund_file['wave_shift'], low=3., high=3.)

    good_lines &= (abund_file['wave_shift'] < clip_high) & (abund_file['wave_shift'] > clip_low)

    if len(abund_file[good_lines]) > 0:
            
        star_param['feh'] = np.round(np.mean(abund_file['logeps_chi2'][good_lines]) - solar_abundances['Fe'], 3)
        write_json(param_filename, star_param)