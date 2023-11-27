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

#star_list = ['hd122563', 'ret2_mem', 'ret2_cand', 'ret2_test']
#star_list = ['hd122563', 'ret2_mem']
#star_list = ['hd122563']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']

inspect=True

solar_abundances = au.solar_abund('A09')

abund_dict = {star: {'abundances': {}, 'errors':{}, 'sigmas': {}, 'nlines': {}} for star in star_list}

element_list = ['Fe', 'C', 'Mg', 'Ti', 'Ca', 'Ni', 'Cr', 'Na', 'K', 'Sc', 'Mn', 'Sr', 'Y', 'Ba', 'Eu']

with open('cra2_abund_used_line_measurements.txt', 'w') as output_line_file:
    
    print('# star element wavelength eqw eqw_err logeps_eqw logeps_err_eqw logeps_chi2', file=output_line_file)

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

            if (star == 'ret2_test') & (element != 'Na'):
                abund_filename2 = f'ret2_cand/abund_no_velshift/ret2_cand_{element}_abund_fit_convol{res:.1f}.txt'
                abund_file2 = np.genfromtxt(abund_filename2, dtype=None, encoding=None, names=True)
                good_lines &= np.in1d(abund_file['wavelength'], abund_file2['wavelength'])


            if element == 'C':
                if (star == 'ret2_cand') or (star == 'ret2_test'):
                    good_lines = np.ones(len(abund_file)).astype(bool)
                else:
                    good_lines = np.in1d(abund_file['wavelength'], abund_lines['C']['CH'])
            if element == 'N':
                if (star == 'ret2_cand') or (star == 'ret2_test'):
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
                for abund_file_line in abund_file[good_lines]:
                    print(f'{star:9s} {element:>2s} {abund_file_line["wavelength"]:8.2f} {abund_file_line["eqw"]:8.2f} {abund_file_line["eqw_err"]:6.2f} {abund_file_line["logeps_eqw"]:8.4f} {abund_file_line["logeps_err_eqw"]:.4f} {abund_file_line["logeps_chi2"]:8.4f}',file=output_line_file)
                abund_dict[star]['abundances'][element] = {}
                abund_dict[star]['errors'][element] = {}
                abund_dict[star]['sigmas'][element] = {}
                abund_dict[star]['nlines'][element] = {}
                for ion in abund_lines[element]:
                    if (ion == 'convol') or (ion == 'C13'):
                        continue
                    ion_lines = np.in1d(abund_file['wavelength'], abund_lines[element][ion])
                    if len(abund_file['logeps_chi2'][good_lines & ion_lines]) < 1:
                        continue
                    abund_dict[star]['abundances'][element][ion] = np.mean(abund_file['logeps_chi2'][good_lines & ion_lines])
                    abund_dict[star]['sigmas'][element][ion] = np.std(abund_file['logeps_chi2'][good_lines & ion_lines])
                    abund_dict[star]['nlines'][element][ion] = len(abund_file['logeps_chi2'][good_lines & ion_lines])
                    if len(abund_file['logeps_chi2'][good_lines & ion_lines]) > 5:
                        abund_dict[star]['errors'][element][ion] = np.std(abund_file['logeps_chi2'][good_lines])/np.sqrt(len(abund_file['logeps_chi2'][good_lines & ion_lines]))
                    else:
                        #abund_dict[star]['errors'][element] = np.sqrt(np.sum(abund_file['logeps_err_eqw'][good_lines]**2.)/len(abund_file['logeps_chi2'][good_lines & ion_lines]))
                        abund_dict[star]['errors'][element][ion] = abund_dict[star]['sigmas']['Fe']['I']/np.sqrt(len(abund_file['logeps_chi2'][good_lines & ion_lines]))



for star_name in star_list:
    abund_dict[star_name]['abundances']['Fe'] = {ion : ion_abund - solar_abundances['Fe'] for ion, ion_abund in abund_dict[star_name]['abundances']['Fe'].items()}
    for element in abund_dict[star_name]['abundances']:
        if element == 'Fe':
            continue
        star_fe_ion_dict = {}
        for ion in abund_dict[star_name]['abundances'][element]:
            if ion in abund_dict[star_name]['abundances']['Fe']:
                star_fe_ion_dict[ion] = abund_dict[star_name]['abundances']['Fe'][ion]
            else:
                star_fe_ion_dict[ion] = abund_dict[star_name]['abundances']['Fe']['I']
        #    star_fe_ion_dict[ion] = abund_dict[star_name]['abundances']['Fe']['I']
        abund_dict[star_name]['abundances'][element] = {ion : ion_abund - solar_abundances[element] - star_fe_ion_dict[ion] for ion, ion_abund in abund_dict[star_name]['abundances'][element].items()}

#for star_name in star_list:
#    abundances = abund_dict[star_name]['abundances']
#    abund_obj = ChemAbund(abundances=abundances, reference='logeps', solar_reference='a09')
#    abund_dict[star_name]['abundances'] = abund_obj.format_abundances('Fe')


elem_list = {'Fe' : []}
for star in star_list:
    for ion in abund_dict[star_name]['abundances']['Fe']:
        if ion not in elem_list['Fe']:
            elem_list['Fe'] += [ion]

for element in au.atomic_symbols:
    if element == 'Fe':
        continue
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

print(elem_list)
print(abund_dict)

with open('cra2_abund_error_table_ions.txt', 'w') as file:
    print('# STAR', file=file, end=' ')
    for element, elem_ion_list in elem_list.items():
        if element == 'Fe':
            for ion in elem_ion_list:
                print(f'FE{ion.upper()}_H FE{ion.upper()}_H_ERR FE{ion.upper()}_H_SIGMA FE{ion.upper()}_N_LINES', file=file, end=' ')
        else:
            for ion in elem_ion_list:
                print(f'{element.upper()}{ion.upper()}_FE {element.upper()}{ion.upper()}_FE_ERR {element.upper()}{ion.upper()}_FE_SIGMA {element.upper()}{ion.upper()}_N_LINES', file=file, end=' ')
    print('', file=file)
    for star_name in star_list:
        print(star_name, file=file, end=' ')
        for element, elem_ion_list in elem_list.items():
            for ion in elem_ion_list:
                try:
                    print(f'{abund_dict[star_name]["abundances"][element][ion]:.3f} {abund_dict[star_name]["errors"][element][ion]:.3f} {abund_dict[star_name]["sigmas"][element][ion]:.3f} {abund_dict[star_name]["nlines"][element][ion]}', file=file, end=' ')
                except:
                    print('NaN NaN NaN 0', file=file, end=' ')
        print('', file=file)
