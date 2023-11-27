import numpy as np
from analysis_utils import read_json, write_json
from scipy.stats import sigmaclip
import abund_utils as au
from chem_abund import ChemAbund
from abund_lines import *
from print_nlte_corrections import nlte_corrections



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

def calculate_lte_abunds():

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
        abund_dict[star_name]['logeps'] = dict(abund_dict[star_name]['abundances'])

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
    
    return star_list, elem_list, abund_dict

def calculate_nlte_abunds():
    
    star_list = ['hd122563', 'ret2_mem', 'ret2_cand', 'ret2_test']
    #star_list = ['hd122563', 'ret2_mem']
    #star_list = ['hd122563']
    star_list = ['cra2_gdr32672', 'cra2_gdr36912']


    inspect=True

    solar_abundances = au.solar_abund('A09')

    abund_dict = {star: {'abundances': {}, 'errors':{}, 'sigmas': {}, 'nlines': {}} for star in star_list}

    element_list = ['Fe', 'C', 'Mg', 'Ti', 'Ca', 'Ni', 'Cr', 'Na', 'K', 'Sc', 'Mn', 'Sr', 'Y', 'Ba', 'Eu']

    nlte_corr_dict = nlte_corrections()

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
                        if star in nlte_corr_dict:
                            if element in nlte_corr_dict[star]:
                                if ion in nlte_corr_dict[star][element]:
                                    abund_dict[star]['abundances'][element][ion] += nlte_corr_dict[star][element][ion]

    for star_name in star_list:
        abund_dict[star_name]['logeps'] = dict(abund_dict[star_name]['abundances'])

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

    return star_list, elem_list, abund_dict


def calculate_sys_errors():

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

    sys_err_dict = {star: {'sys_errors': {}} for star in star_list}

    base_abunds = abund_dict_list[0]
    diff_abunds = abund_dict_list[1:]

    for star in star_list:
        for element in elem_list:
            sys_err_dict[star]['sys_errors'][element] = {}
            if element in base_abunds[star]['abundances']:
                for ion in base_abunds[star]['abundances'][element]:
                    total = 0.
                    for i, diff_param in enumerate(['teff', 'logg', 'vmicro', 'mh']):
                        logeps_base = base_abunds[star]['abundances'][element][ion]
                        logeps_diff = diff_abunds[i][star]['abundances'][element][ion]
                        total = np.hypot(total,logeps_diff-logeps_base)
                        sys_err_dict[star]['sys_errors'][element][ion] = total

    return sys_err_dict


lte_star_list, lte_elem_list, lte_abund_dict = calculate_lte_abunds()
nlte_star_list, nlte_elem_list, nlte_abund_dict = calculate_nlte_abunds()
sys_err_dict = calculate_sys_errors()
    

with open('cra2_abund_latex_table.txt', 'w') as file:

    for element, elem_ion_list in lte_elem_list.items():
        for ion in elem_ion_list:
            print(f'\ion[{element}][{ion}]', file=file, end=' ')
            for star_name in lte_star_list:
                try:
                    print(f'& {lte_abund_dict[star_name]["nlines"][element][ion]} & {lte_abund_dict[star_name]["logeps"][element][ion]:.2f} & {lte_abund_dict[star_name]["abundances"][element][ion]:.2f} & {nlte_abund_dict[star_name]["logeps"][element][ion]:.2f} & {nlte_abund_dict[star_name]["abundances"][element][ion]:.2f} & {lte_abund_dict[star_name]["sigmas"][element][ion]:.2f} & {lte_abund_dict[star_name]["errors"][element][ion]:.2f}  & {sys_err_dict[star_name]["sys_errors"][element][ion]:.2f} & {np.hypot(lte_abund_dict[star_name]["errors"][element][ion],sys_err_dict[star_name]["sys_errors"][element][ion]):.2f}', file=file, end=' ')
                except:
                    print('& \\nodata & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata', file=file, end=' ')
            print('\\\\', file=file)

with open('cra2_abund_error_table_ions_all.txt', 'w') as file:
    print('# STAR', file=file, end=' ')
    for element, elem_ion_list in lte_elem_list.items():
        if element == 'Fe':
            for ion in elem_ion_list:
                print(f'FE{ion.upper()}_H FE{ion.upper()}_H_NLTE FE{ion.upper()}_H_ERR_STAT FE{ion.upper()}_H_ERR_SYS FE{ion.upper()}_H_ERR_TOT FE{ion.upper()}_H_SIGMA FE{ion.upper()}_N_LINES', file=file, end=' ')
        else:
            for ion in elem_ion_list:
                print(f'{element.upper()}{ion.upper()}_FE {element.upper()}{ion.upper()}_FE_NLTE {element.upper()}{ion.upper()}_FE_ERR_STAT {element.upper()}{ion.upper()}_FE_ERR_SYS {element.upper()}{ion.upper()}_FE_ERR_TOT {element.upper()}{ion.upper()}_FE_SIGMA {element.upper()}{ion.upper()}_N_LINES', file=file, end=' ')
    print('', file=file)
    for star_name in lte_star_list:
        print(star_name, file=file, end=' ')
        for element, elem_ion_list in lte_elem_list.items():
            for ion in elem_ion_list:
                try:
                    print(f'{lte_abund_dict[star_name]["abundances"][element][ion]:.3f} {nlte_abund_dict[star_name]["abundances"][element][ion]:.3f} {lte_abund_dict[star_name]["errors"][element][ion]:.3f} {sys_err_dict[star_name]["sys_errors"][element][ion]:.3f} {np.hypot(lte_abund_dict[star_name]["errors"][element][ion],sys_err_dict[star_name]["sys_errors"][element][ion]):.3f} {lte_abund_dict[star_name]["sigmas"][element][ion]:.3f} {lte_abund_dict[star_name]["nlines"][element][ion]}', file=file, end=' ')
                except:
                    print('NaN NaN NaN NaN NaN NaN 0', file=file, end=' ')
        print('', file=file)