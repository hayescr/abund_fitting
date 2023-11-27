import numpy as np
from analysis_utils import read_json, write_json
from scipy.stats import sigmaclip
import abund_utils as au

def identify_good_lines(abund_file, element, detection_threshold=5.0):

    good_lines = (abund_file['lim_flag_eqw'].astype('bool')) & (
        abund_file['lim_flag_chi2'].astype('bool'))
    good_lines &= abund_file['avg_synth_stdev'] / \
        (1.5/(abund_file['snr'] * np.sqrt(abund_file['npix']))) > 1.0
    good_lines &= (abund_file['avg_line_depth'] * abund_file['snr'] > detection_threshold)

    good_lines &= (np.abs(abund_file['wave_shift']) < 0.1)

    if not np.all(abund_file['wave_shift'] == 0.0):
        a, clip_low, clip_high = sigmaclip(
            abund_file['wave_shift'], low=3., high=3.)

        good_lines &= (abund_file['wave_shift'] < clip_high) & (
            abund_file['wave_shift'] > clip_low)

    return good_lines

#star_list = ['hd122563', 'ret2_mem', 'ret2_cand']
#star_list = ['hd122563', 'ret2_mem']
#star_list = ['hd122563']
#star_list = ['cra2_gdr36912']
star_list = ['sgr2_gdr35584','sgr2_gdr32656','sgr2_gdr39936','aqu2_gdr33776','aqu2_gdr31472']
element = 'C'

inspect=False

detection_threshold = 4.

solar_abundances = au.solar_abund('A09')

for star_name in star_list:
    star = star_name
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)
    res = star_param['vmacro']

    abund_filename = f'{star}/abund_no_velshift/{star}_{element}_abund_fit_convol{res:.1f}.txt'
    print(abund_filename)
    abund_file = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)

    if abund_file.shape == ():
        abund_file = np.array([abund_file])
        
    good_lines = identify_good_lines(abund_file, element, detection_threshold=detection_threshold)

    if element == 'C':
        #abund_file = np.array([abund_file])
        good_lines = np.ones(len(abund_file)).astype(bool)
    if element == 'N':
        if star == 'ret2_cand':
            good_lines = np.ones(len(abund_file)).astype(bool)
            good_lines[abund_file['wavelength'] == 8200.00] = False
        else:
            good_lines = np.zeros(len(abund_file)).astype(bool)

    if element == 'Ca':
        bad_lines = (abund_file['wavelength'] == 8498.02)
        bad_lines |= (abund_file['wavelength'] == 8542.09)
        bad_lines |= (abund_file['wavelength'] == 8662.14)
        good_lines &= ~bad_lines


    if len(abund_file[good_lines]) > 0:
        print(star_name)
        print(abund_file['wavelength'][good_lines])
        print(abund_file['logeps_chi2'][good_lines])
        print(np.mean(abund_file['logeps_chi2'][good_lines]))
        if element != 'Fe':
            print(np.round(np.mean(abund_file['logeps_chi2'][good_lines]) - solar_abundances[element] - star_param['abundances']['Fe'] + solar_abundances['Fe'], 3))
        if not inspect:
            if element == 'C':
                star_param['cfe'] = np.round(np.mean(abund_file['logeps_chi2'][good_lines]) - solar_abundances[element] - star_param['feh'], 3)
            if element == 'Mg':
                star_param['alphafe'] = np.round(np.mean(abund_file['logeps_chi2'][good_lines]) - solar_abundances[element] - star_param['feh'], 3)
            if element == 'N':
                star_param['nfe'] = np.round(np.mean(abund_file['logeps_chi2'][good_lines]) - solar_abundances[element] - star_param['feh'], 3)
            star_param['abundances'][element] = np.mean(abund_file['logeps_chi2'][good_lines])
            write_json(param_filename, star_param)