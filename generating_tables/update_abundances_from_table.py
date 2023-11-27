import numpy as np
from analysis_utils import read_json, write_json
from scipy.stats import sigmaclip
import abund_utils as au
from chem_abund import ChemAbund
from abund_lines import *



def identify_good_lines(abund_file, element):

    good_lines = (abund_file['lim_flag_eqw'].astype('bool')) & (
        abund_file['lim_flag_chi2'].astype('bool'))
    good_lines &= abund_file['avg_synth_stdev'] / \
        (1.5/(abund_file['snr'] * np.sqrt(abund_file['npix']))) > 1.0
    good_lines &= (abund_file['avg_line_depth'] * abund_file['snr'] > 5.)

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

star_list = ['cra2_gdr32672', 'cra2_gdr36912']

inspect=True

solar_abundances = au.solar_abund('A09')

abund_dict = {star: {'abundances': {}, 'errors':{}, 'sigmas': {}, 'nlines': {}} for star in star_list}

element_list = ['Fe', 'C', 'Mg', 'Ti', 'Ca', 'Ni', 'Cr', 'Na', 'K', 'Sc', 'Mn', 'Sr', 'Y', 'Ba', 'Eu']

for star_name in star_list:
    star = star_name
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)
    res = star_param['vmacro']

    for element in element_list:

        abund_filename = f'{star}/abund_no_velshift/{star}_{element}_abund_fit_convol{res:.1f}.txt'
        abund_file = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)
        
        if abund_file.shape == ():
            abund_file = np.array([abund_file])

        good_lines = identify_good_lines(abund_file, element)

        if element == 'C':
            if star in ['ret2_cand', 'ret2_test']:
                good_lines = np.ones(len(abund_file)).astype(bool)
            else:
                good_lines = np.in1d(abund_file['wavelength'], abund_lines['C']['CH'])
        if element == 'N':
            if star in ['ret2_cand', 'ret2_test']:
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
            abund_dict[star]['abundances'][element] = np.mean(abund_file['logeps_chi2'][good_lines])
            abund_dict[star]['sigmas'][element] = np.std(abund_file['logeps_chi2'][good_lines])
            abund_dict[star]['nlines'][element] = len(abund_file['logeps_chi2'][good_lines])
            if len(abund_file['logeps_chi2'][good_lines]) > 5:
                 abund_dict[star]['errors'][element] = np.std(abund_file['logeps_chi2'][good_lines])/np.sqrt(len(abund_file['logeps_chi2'][good_lines]))
            else:
                 #abund_dict[star]['errors'][element] = np.sqrt(np.sum(abund_file['logeps_err_eqw'][good_lines]**2.)/len(abund_file['logeps_chi2'][good_lines]))
                 abund_dict[star]['errors'][element] = abund_dict[star]['sigmas']['Fe']/np.sqrt(len(abund_file['logeps_chi2'][good_lines]))

    abundances = abund_dict[star_name]['abundances']
    star_param['abundances'] = abundances
    write_json(param_filename, star_param)

