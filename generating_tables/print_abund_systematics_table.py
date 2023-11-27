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

    good_lines &= (abund_file['avg_line_depth'] * abund_file['snr'] > 4.)

    good_lines &= (np.abs(abund_file['wave_shift']) < 0.1)

    if not np.all(abund_file['wave_shift'] == 0.0):
        a, clip_low, clip_high = sigmaclip(
            abund_file['wave_shift'], low=3., high=3.)

        good_lines &= (abund_file['wave_shift'] < clip_high) & (
            abund_file['wave_shift'] > clip_low)

    return good_lines

star_list = ['hd122563', 'ret2_mem', 'ret2_test']
#star_list = ['hd122563', 'ret2_mem']
#star_list = ['hd122563']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']

measured_lines_all = np.genfromtxt('cra2_abund_used_line_measurements.txt', encoding=None, dtype=None, names=True)

inspect=True

solar_abundances = au.solar_abund('A09')

element_list = ['Fe', 'C', 'Mg', 'Ti', 'Ca', 'Ni', 'Cr', 'Na', 'K', 'Sc', 'Mn', 'Sr', 'Y', 'Ba', 'Eu']

delta_params = ['abund_no_velshift', 'teff_plus_err', 'logg_plus_err', 'vmicro_plus_err', 'mh_plus_err']

abund_dict_list = []

for delta_param in delta_params:

    abund_dict = {star: {'abundances': {}, 'errors':{}, 'sigmas': {}, 'nlines': {}} for star in star_list}


    for star_name in star_list:
        star = star_name
        param_filename = f'{star_name}/{star_name}_params.txt'
        star_param = read_json(param_filename)
        res = star_param['vmacro']

        for element in element_list:

            abund_filename = f'{star}/{delta_param}/{star}_{element}_abund_fit_convol{res:.1f}.txt'
            abund_file = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)

            if abund_file.shape == ():
                abund_file = np.array([abund_file])

            #good_lines = identify_good_lines(abund_file, element)

            good_lines = np.in1d(abund_file['wavelength'], measured_lines_all['wavelength'][(measured_lines_all['star']==star_name)&(measured_lines_all['element']==element)])


            #if element == 'C':
            #    if (star == 'ret2_cand') or (star == 'ret2_test'):
            #        good_lines = np.ones(len(abund_file)).astype(bool)
            #    else:
            #        good_lines = np.in1d(abund_file['wavelength'], abund_lines['C']['CH'])
            #if element == 'N':
            #    if (star == 'ret2_cand') or (star == 'ret2_test'):
            #        good_lines = np.ones(len(abund_file)).astype(bool)
            #        good_lines[abund_file['wavelength'] == 8200.00] = False
            #    else:
            #        good_lines = np.zeros(len(abund_file)).astype(bool)

            #if element == 'Ca':
            #    bad_lines = (abund_file['wavelength'] == 8498.02)
            #    bad_lines |= (abund_file['wavelength'] == 8542.09)
            #    bad_lines |= (abund_file['wavelength'] == 8662.14)
            #    good_lines &= ~bad_lines

            if len(abund_file[good_lines]) > 0:
                abund_dict[star]['abundances'][element] = {}
                for ion in abund_lines[element]:
                    if (ion == 'convol') or (ion == 'C13'):
                        continue
                    ion_lines = np.in1d(abund_file['wavelength'], abund_lines[element][ion])
                    if len(abund_file['logeps_chi2'][good_lines & ion_lines]) < 1:
                        continue
                    abund_dict[star]['abundances'][element][ion] = np.mean(abund_file['logeps_chi2'][good_lines & ion_lines])


    elem_list = {}
    for element in au.atomic_symbols:
        elem_flag = False
        ion_list = []
        for star_name in star_list:
            if element in abund_dict[star_name]['abundances']:
                elem_flag=True
                for ion in abund_dict[star_name]['abundances'][element]:
                    if ion not in ion_list:
                        ion_list += [ion]
        if elem_flag:
            elem_list[element] = ion_list

    abund_dict_list += [abund_dict]


with open('cra2_systematic_errors.txt', 'w') as file:
    base_abunds = abund_dict_list[0]
    diff_abunds = abund_dict_list[1:]

    print('# Element Species delta_teff delta_logg delta_vmicro delta_mh sigma_sys', file=file)

    for star in star_list:
        print(star, file=file)
        for element in elem_list:
            if element in base_abunds[star]['abundances']:
                for ion in base_abunds[star]['abundances'][element]:
                    print(f'{element+" "+ion:<6s}', end=' ', file=file)
                    total = 0.
                    for i, diff_param in enumerate(['teff', 'logg', 'vmicro', 'mh']):
                        logeps_base = base_abunds[star]['abundances'][element][ion]
                        logeps_diff = diff_abunds[i][star]['abundances'][element][ion]
                        total = np.hypot(total,logeps_diff-logeps_base)
                        print(f'{round(logeps_diff-logeps_base,2):>5.2f}', end=' ', file=file)
                    print(f'{total:>5.2f}', end=' ', file=file)
                    print('', file=file)
