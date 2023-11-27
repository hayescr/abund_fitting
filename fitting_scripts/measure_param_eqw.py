import numpy as np
from linelist_manager import LinelistManager
from abundometer import Abundometer, Eqwometer
from turbospec_wrapper.turbosynth import TurboSynth
from spec_tools import Spectrum
import abund_utils as au
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.cm as cmx
import os
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from analysis_utils import read_json, write_json
from abund_lines import *
from scipy.stats import sigmaclip, median_abs_deviation
from scipy.optimize import curve_fit
import abund_utils as au
import logg_SB

solar_abundances = au.solar_abund('a09')


def vmicro_mash(teff, logg, feh):
    return 0.14 - (0.08 * feh) + (4.90 * teff/1.0e4) - (0.47 * logg)


def calculate_sb_logg(star_name, teff, teff_err):

    data = np.genfromtxt('stellar_phot_table.txt',
                         names=True, dtype=None, encoding=None)

    data = data[data['name'] == star_name]
    ra = data['ra']
    dec = data['dec']
    ebv = data['ebv']

    av = 3.1*ebv

    ag = av*0.85926
    abp = av*1.06794
    arp = av*0.65199

    dist = data['dist']
    dist_err = data['dist_err']
    G = data['phot_g_mean_mag']
    C = data['bp_rp'] - (abp - arp)  # bp-rp

    teff = np.array([teff])
    teff_err = np.array([teff_err])

    logg_sb, logg_sb_err = logg_SB.logg_SB(
        teff, teff_err, dist, dist_err, G, ag, 1000)

    return logg_sb[0][0], logg_sb_err[0][0]


def initiate_turbosynth(teff, logg, feh, vmicro, cfe, nfe, alphafe, abundances):

    solar_abund = au.solar_abund()
    abundances = abundances
    if 'N' not in abundances:
        abundances['N'] = solar_abund['N'] + feh + nfe
    if 'C' not in abundances:
        abundances['C'] = solar_abund['C'] + feh + cfe

    star = 'test'

    rel_path = '../../Turbospectrum2019/COM-v19.1/'

    linelists = ['DATA/Hlinedata',
                 'new_linelists/kurucz_gfallvac08oct17_300-1100.atoms.turbo.air_original',
                 'new_linelists/12CH_mas14.turbo.air',
                 'new_linelists/13CH_mas14.turbo.air',
                 'new_linelists/24Mg1H_kurucz_all.turbo.air',
                 'new_linelists/16O1H_exomol.turbo.air',
                 'new_linelists/12C14N_exomol.turbo.air',
                 'new_linelists/14N1H_exomol.turbo.air',
                 'new_linelists/12C12C_boltz85_exomol.turbo.air',
                 'new_linelists/28Si1H_exomol.turbo.air',
                 'new_linelists/56Fe1H_exomol.turbo.air',
                 'new_linelists/40Ca1H_exomol.turbo.air',
                 'new_linelists/linemake.euII.mooghfs.atoms.turbo.air'
                 ]

    #linelists = ['DATA/Hlinedata',
    #'bacchus_linelists/atoms_4200-9200.list']

    linelists = [f'{rel_path}{linelist}' for linelist in linelists]

    test = TurboSynth(path='tmp2', turbopath='../../Turbospectrum2019')
    test.init_params(teff, logg, feh, vmicro=vmicro, cfe=cfe, alphafe=alphafe)
    test.init_abunds(abunds=abundances, solar_reference='Asplund2009')
    #test.init_abunds(solar_reference='Asplund2005')
    #test.set_abund('C', test.get_abund('C') + cfe)
    #test.set_abund('N', test.get_abund('N') + cfe)

    test.set_linelists(linelists)
    return test


def write_linelist(ion, line, line_data):
    if ion == '1':
        ts_headers = ["'26.0000             '    1         1",
                      "'Fe I   '"]
    elif ion == '2':
        ts_headers = ["'26.0000             '    2         1",
                      "'Fe II  '"]

    filename = 'linelist.tmp'
    with open(filename, 'w') as tmp_linelist_file:
        for header in ts_headers:
            print(header, file=tmp_linelist_file)
        print(line_data, file=tmp_linelist_file)
    return filename


def measure_synthetic_eqw(mem, star_param, ion, lines, abunds, line_param_data):
    ex_turbosynth2 = initiate_turbosynth(star_param['teff'],
                                         star_param['logg'],
                                         star_param['feh'],
                                         star_param['vmicro'],
                                         star_param['cfe'],
                                         star_param['nfe'],
                                         star_param['alphafe'], star_param['abundances'])

    res = star_param['vmacro']
    if not os.path.isdir(star):
        os.mkdir(star)

    conv_profiles = [2]
    conv_broads = [{'vel': res}]

    synth_line_eqws = []

    for line, abund, line_data in zip(lines, abunds, line_param_data):
        linelist_name = write_linelist(ion, line, line_data)
        ex_turbosynth2.set_linelists([linelist_name])
        ex_turbosynth2.set_abund('Fe', abund)

        eqw_fitter = Eqwometer(ex_turbosynth2, mem, 'Fe', line,
                               conv_profiles, conv_broads,
                               synth_range=1, delta_abunds=[0.])
        line_synth_eqw = eqw_fitter.measure_synth_eqw()[0]
        synth_line_eqws += [line_synth_eqw]

    return np.array(synth_line_eqws)


def measure_eqw_abunds(mem, star_param, ion, lines,
                       abunds, eqws, line_param_data):

    ex_turbosynth2 = initiate_turbosynth(star_param['teff'],
                                         star_param['logg'],
                                         star_param['feh'],
                                         star_param['vmicro'],
                                         star_param['cfe'],
                                         star_param['nfe'],
                                         star_param['alphafe'],
                                         star_param['abundances'])

    if star_name == 'ret2_cand':
        heavy_alphas = ['Si', 'S', 'Ar', 'Ca', 'Ti']
        for alpha in heavy_alphas:
            ex_turbosynth2.set_abund(
                alpha, star_param['feh'] + 0.4 + solar_abundances[alpha])

    res = star_param['vmacro']
    if not os.path.isdir(star):
        os.mkdir(star)

    iterations = 2

    conv_profiles = [2]
    conv_broads = [{'vel': res}]

    synth_eqw_abunds = []

    for line, abund, eqw, line_data in zip(lines, abunds, eqws, line_param_data):
        linelist_name = write_linelist(ion, line, line_data)
        ex_turbosynth2.set_linelists([linelist_name])
        ex_turbosynth2.set_abund('Fe', abund)

        eqw_fitter = Eqwometer(ex_turbosynth2, mem, 'Fe', line,
                               conv_profiles, conv_broads,
                               synth_range=1)
        line_synth_eqw = eqw_fitter.fit(eqw)

        for i in range(iterations-1):
            new_delta_abund = eqw_fitter.fit_abund_eqw - \
                ex_turbosynth2.get_abund('Fe')
            if (new_delta_abund <= eqw_fitter.delta_abunds[1]) | (
                    new_delta_abund >= eqw_fitter.delta_abunds[-2]):

                eqw_fitter = Abundometer(ex_turbosynth2, mem, 'Fe', line,
                                         conv_profiles, conv_broads,
                                         synth_range=1,
                                         delta_abunds=new_delta_abund
                                         + np.arange(-0.6, 0.61, 0.3),
                                         )
                eqw_fitter.fit()

        synth_eqw_abunds += [eqw_fitter.fit_abund_eqw]

    return np.array(synth_eqw_abunds)


fe1_ts_header = ["'26.0000             '    1         1",
                 "'Fe I   '"]

fe2_ts_header = ["'26.0000             '    2         1",
                 "'Fe II  '"]

if not os.path.isdir('tmp'):
    os.mkdir('tmp')

# star_list = ['ret2_cand', 'ret2_mem']
# star_list = ['ret2_cand']
# star_list = ['ret2_mem']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']
# star_list = ['hd122563', 'ret2_mem']
# star_list = ['hd122563']

wavelength_limit = 4000.


method = 'eqw'
clip_abunds = True
fix_vmicro = True
fix_logg = True
update_param = False

params_to_update = ['teff', 'feh', 'vmicro']
# params_to_update = ['logg']
# params_to_update = ['teff', 'logg', 'feh', 'vmicro']
# params_to_update = ['teff']

params_to_update = ['teff', 'logg']
params_to_update = ['feh']


print('logg' in params_to_update)

for star_name in star_list:

    # Reading things in

    star = star_name

    star_data = np.genfromtxt(
        f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
    mem = Spectrum(star_data['wavelength'],
                   star_data['flux'], variance=star_data['variance'])
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)

    teff_i = star_param['teff']
    logg_i = star_param['logg']
    feh_i = star_param['feh']
    vmicro_i = star_param['vmicro']
    res = star_param['vmacro']

    abund_filename = f'{star}/{star}_Fe_abund_fit_convol{res:.1f}.txt'

    line_meas = np.genfromtxt(abund_filename, dtype=None,
                              encoding=None, names=True)
    fe1_line_data = np.genfromtxt('Fe_I_ts_lines.txt', dtype=None,
                                  usecols=(0, 1, 2),
                                  names=('wavelength', 'ep', 'loggf'),
                                  encoding=None)
    fe2_line_data = np.genfromtxt('Fe_II_ts_lines.txt', dtype=None,
                                  usecols=(0, 1, 2),
                                  names=('wavelength', 'ep', 'loggf'),
                                  encoding=None)

    # Cleaning good lines

    fe1_lines = np.array(abund_lines['Fe']['I'])
    fe2_lines = np.array(abund_lines['Fe']['II'])

    good_lines = (line_meas['lim_flag_eqw'].astype('bool')) & (
        line_meas['lim_flag_chi2'].astype('bool'))
    good_lines &= line_meas['avg_synth_stdev'] / \
        (1.5/(line_meas['snr'] * np.sqrt(line_meas['npix']))) > 1.0
    good_lines &= (line_meas['avg_line_depth'] * line_meas['snr'] > 5.)

    #a, clip_low, clip_high = sigmaclip(
    #    line_meas['wave_shift'], low=3., high=3.)

    #good_lines &= (line_meas['wave_shift'] < clip_high) & (
    #    line_meas['wave_shift'] > clip_low)

    line_meas = line_meas[good_lines]

    # Matching line data

    fe1_meas = line_meas[np.in1d(line_meas['wavelength'], fe1_lines)]
    fe2_meas = line_meas[np.in1d(line_meas['wavelength'], fe2_lines)]

    fe1_meas_ep = []
    fe1_line_index = []
    fe2_meas_ep = []
    fe2_line_index = []

    for wave in fe1_meas['wavelength']:
        line_match = np.argmin(np.abs(wave-fe1_line_data['wavelength']))
        fe1_line_index += [line_match]
        fe1_meas_ep += [fe1_line_data['ep'][line_match]]
    fe1_meas_ep = np.array(fe1_meas_ep)

    for wave in fe2_meas['wavelength']:
        line_match = np.argmin(np.abs(wave-fe2_line_data['wavelength']))
        fe2_line_index += [line_match]
        fe2_meas_ep += [fe2_line_data['ep'][line_match]]
    fe2_meas_ep = np.array(fe2_meas_ep)

    fe1_abunds = fe1_meas['logeps_eqw']
    fe2_abunds = fe2_meas['logeps_eqw']

    with open('Fe_I_ts_lines.txt') as lp_file:
        all_line_entries = lp_file.readlines()
        fe1_line_param_data = [
            all_line_entries[line_entry] for line_entry in fe1_line_index]

    with open('Fe_II_ts_lines.txt') as lp_file:
        all_line_entries = lp_file.readlines()
        fe2_line_param_data = [
            all_line_entries[line_entry] for line_entry in fe2_line_index]

    # Getting the synthetic eqws

    fe1_eqw_filename = f'{star_name}/{star_name}_fe1_eqws.txt'
    if os.path.exists(fe1_eqw_filename):
        fe1_eqw_file = np.genfromtxt(fe1_eqw_filename, names=True,
                                     encoding=None, dtype=None)
        fe1_eqws = fe1_eqw_file['eqw']
    else:
        fe1_eqws = measure_synthetic_eqw(mem, star_param, '1', fe1_meas['wavelength'],
                                         fe1_abunds, fe1_line_param_data)
        with open(fe1_eqw_filename, 'w') as eqw_file:
            print('# wavelength eqw', file=eqw_file)
            for wave, eqw in zip(fe1_meas['wavelength'], fe1_eqws):
                print(wave, eqw, file=eqw_file)

    fe2_eqw_filename = f'{star_name}/{star_name}_fe2_eqws.txt'
    if os.path.exists(fe2_eqw_filename):
        fe2_eqw_file = np.genfromtxt(fe2_eqw_filename, names=True,
                                     encoding=None, dtype=None)
        fe2_eqws = fe2_eqw_file['eqw']
    else:
        fe2_eqws = measure_synthetic_eqw(mem, star_param, '2', fe2_meas['wavelength'],
                                         fe2_abunds, fe2_line_param_data)
        with open(fe2_eqw_filename, 'w') as eqw_file:
            print('# wavelength eqw', file=eqw_file)
            for wave, eqw in zip(fe2_meas['wavelength'], fe2_eqws):
                print(wave, eqw, file=eqw_file)

    # Use synthetic eqws to get abundances

    fe1_eqw_abund_filename = f'{star_name}/{star_name}_t{teff_i:.0f}_g{logg_i:2f}_m{feh_i}_v{vmicro_i}_param_fe1_eqw_abunds.txt'
    if os.path.exists(fe1_eqw_abund_filename):
        fe1_eqw_abund_file = np.genfromtxt(fe1_eqw_abund_filename, names=True,
                                           encoding=None, dtype=None)
        fe1_eqw_abunds = fe1_eqw_abund_file['logeps']
    else:
        fe1_eqw_abunds = measure_eqw_abunds(mem, star_param, '1', fe1_meas['wavelength'],
                                            fe1_abunds, fe1_eqws, fe1_line_param_data)
        with open(fe1_eqw_abund_filename, 'w') as eqw_abund_file:
            print('# wavelength eqw logeps', file=eqw_abund_file)
            for wave, eqw, abund in zip(fe1_meas['wavelength'], fe1_eqws, fe1_eqw_abunds):
                print(wave, eqw, abund, file=eqw_abund_file)

    fe2_eqw_abund_filename = f'{star_name}/{star_name}_t{teff_i:.0f}_g{logg_i:2f}_m{feh_i}_v{vmicro_i}_param_fe2_eqw_abunds.txt'

    if os.path.exists(fe2_eqw_abund_filename):
        fe2_eqw_abund_file = np.genfromtxt(fe2_eqw_abund_filename, names=True,
                                           encoding=None, dtype=None)
        fe2_eqw_abunds = fe2_eqw_abund_file['logeps']
    else:
        fe2_eqw_abunds = measure_eqw_abunds(mem, star_param, '2', fe2_meas['wavelength'],
                                            fe2_abunds, fe2_eqws, fe2_line_param_data)
        with open(fe2_eqw_abund_filename, 'w') as eqw_abund_file:
            print('# wavelength eqw logeps', file=eqw_abund_file)
            for wave, eqw, abund in zip(fe2_meas['wavelength'], fe2_eqws, fe2_eqw_abunds):
                print(wave, eqw, abund, file=eqw_abund_file)

    # Setting up inputs

    fe1_norm_eqw = fe1_eqws/fe1_meas['wavelength']
    fe2_norm_eqw = fe2_eqws/fe2_meas['wavelength']

    if clip_abunds:
        a, low_cut, high_cut = sigmaclip(fe1_eqw_abunds, low=3, high=3)
        selection = (fe1_eqw_abunds > low_cut) & (fe1_eqw_abunds < high_cut)
        fe1_eqw_abunds = fe1_eqw_abunds[selection]
        fe1_norm_eqw = fe1_norm_eqw[selection]
        fe1_meas_ep = fe1_meas_ep[selection]

        # a, low_cut, high_cut = sigmaclip(fe2_eqw_abunds, low=3, high=3)
        # print(low_cut, high_cut)
        # selection = (fe2_eqw_abunds > low_cut) & (fe2_eqw_abunds < high_cut)
        # fe2_eqw_abunds = fe2_eqw_abunds[selection]
        # fe2_norm_eqw = fe2_norm_eqw
        # fe2_meas_ep = fe2_meas_ep[selection]

    def line(x, a, b):
        return a*x+b

    # Fitting relations

    ep_par, ep_cov = curve_fit(line, fe1_meas_ep, fe1_eqw_abunds)
    eqw_par, eqw_cov = curve_fit(
        line, fe1_norm_eqw[fe1_norm_eqw < 0.05], fe1_eqw_abunds[fe1_norm_eqw < 0.05])

    ion_abund_diff = np.mean(fe1_eqw_abunds) - np.mean(fe2_eqw_abunds)

    fe1_sigma = np.std(fe1_eqw_abunds)/(len(fe1_eqw_abunds) - 1)
    fe2_sigma = np.std(fe2_eqw_abunds)/(len(fe2_eqw_abunds) - 1)
    if len(fe2_eqw_abunds) < 4:
        ion_abund_sigma = np.std(np.concatenate(
            (fe1_eqw_abunds, fe2_eqw_abunds)))
    else:
        ion_abund_sigma = np.hypot(fe1_sigma, fe2_sigma)

    #ion_abund_diff = np.median(fe2_eqw_abunds) - np.median(fe1_eqw_abunds)

    #fe1_sigma = median_abs_deviation(fe1_eqw_abunds)/(len(fe1_eqw_abunds) -1)
    #fe2_sigma = median_abs_deviation(fe2_eqw_abunds)/(len(fe2_eqw_abunds) -1)
    #ion_abund_sigma = median_abs_deviation(np.concatenate((fe1_eqw_abunds,fe2_eqw_abunds)))
    #ion_abund_sigma = np.hypot(fe1_sigma, fe2_sigma)

    #############
    # Checking parameters
    #############

    # Teff
    print()
    teff_changed = False
    if (ep_par[0] < np.sqrt(ep_cov[0, 0])) and (ep_par[0] > -np.sqrt(ep_cov[0, 0])):
        print(f'Teff good : ds/dx = {ep_par[0]} +/- {np.sqrt(ep_cov[0,0])}')
    else:
        delta_teff = np.sign(ep_par[0]) * \
                             np.minimum(np.abs(ep_par[0]/0.0002), 1000.)
        new_teff = np.round(star_param["teff"]+delta_teff, 2)
        print(
            f'Teff bad : ds/dx = {ep_par[0]:.6f} +/- {np.sqrt(ep_cov[0,0]):.6f}')
        print(f'Delta Teff = {delta_teff:.2f}')
        print(
            f'Old Teff = {star_param["teff"]:.2f} : New Teff = {new_teff:.2f}')
        if 'teff' in params_to_update:
            print('Teff updated.')
            star_param['teff'] = new_teff
            teff_changed = True

    # Logg
    print()
    if (fix_logg):
        print('Using Logg from Stefan-Boltzman')
        new_logg, logg_err = calculate_sb_logg(star, star_param['teff'],
                                               np.sqrt(ep_cov[0, 0]/0.0002))
        print(
            f'Old Logg = {star_param["logg"]:.2f} : New Logg = {new_logg:.2f} : Error = {logg_err:.2f}')
        if teff_changed:
            print('Logg updated')
            star_param['logg'] = new_logg
    elif (ion_abund_diff < ion_abund_sigma) and (ion_abund_diff > -ion_abund_sigma):
        print(
            f'Logg good : delta eps(I - II) = {ion_abund_diff:.4f} +/- {ion_abund_sigma:.4f}')
    else:
        delta_logg = ion_abund_diff / 0.7
        new_logg = np.round(star_param["logg"] + delta_logg, 4)
        print(
            f'Logg bad : delta eps(I - II) = {ion_abund_diff:.4f} +/- {ion_abund_sigma:.4f}')
        print(f'Delta logg = {delta_logg:.2f}')
        print(
            f'Old logg = {star_param["logg"]:.2f} : New logg = {new_logg:.2f}')
        if 'logg' in params_to_update:
            print('logg updated.')
            star_param['logg'] = new_logg

    # [Fe/H]
    print()
    new_feh = np.round(np.mean(fe1_eqw_abunds) - solar_abundances['Fe'], 3)
    print(f'Old [Fe/H] = {star_param["feh"]:.2f} : New [Fe/H] = {new_feh:.3f}')
    if 'feh' in params_to_update:
        print('[Fe/H] updated.')
        star_param['feh'] = new_feh
        star_param['abundances']['Fe'] = np.mean(fe1_eqw_abunds)

    # Vmicro
    print()
    if (fix_vmicro):
        print('Using Vmicro relation from Mashenina et al. 2017')
        new_vmicro = np.round(vmicro_mash(
            star_param['teff'], star_param['logg'], star_param['feh']), 2)
        print(
            f'Old Vmicro = {star_param["vmicro"]:.2f} : New Vmicro = {new_vmicro:.2f}')
        star_param['vmicro'] = new_vmicro
    elif (eqw_par[0] < np.sqrt(eqw_cov[0, 0])) & (eqw_par[0] > -np.sqrt(eqw_cov[0, 0])):
        print(
            f'Vmicro good : ds/deqw = {eqw_par[0]:.6f} +/- {np.sqrt(eqw_cov[0,0]):.6f}')
    else:
        delta_vmicro = eqw_par[0]/42.
        new_vmicro = np.round(star_param["vmicro"]+delta_vmicro, 2)
        print(
            f'Vmicro bad : ds/deqw = {eqw_par[0]:.6f} +/- {np.sqrt(eqw_cov[0,0]):.6f}')
        print(f'Delta Vmicro = {delta_vmicro:.2f}')
        print(
            f'Old Vmicro = {star_param["vmicro"]:.2f} : New Vmicro = {new_vmicro:.2f}')
        if 'vmicro' in params_to_update:
            print('Vmicro updated.')
            star_param['vmicro'] = new_vmicro

    if update_param:
        if (star == 'ret2_cand') & (star_param['feh'] < -2.5):
            star_param['feh'] = -2.5
        write_json(param_filename, star_param)

    ############
    # Plotting
    ############

    fig = plt.figure(figsize=(8, 7))

    ax = fig.add_subplot(211)

    ax.plot(fe1_meas_ep, fe1_eqw_abunds, marker='.', color='k',
            ls='None', label=f'Fe I:  logeps = {np.mean(fe1_eqw_abunds):.3f}')
    ax.plot(fe2_meas_ep, fe2_eqw_abunds, marker='.', color='r',
            ls='None', label=f'Fe II:  logeps = {np.mean(fe2_eqw_abunds):.3f}')

    ax.plot(np.sort(fe1_meas_ep), line(np.sort(fe1_meas_ep), *ep_par), color='k', ls='dashed',
            label=f'ds/dx = {ep_par[0]:.4f} +/- {np.sqrt(ep_cov[0,0]):.4f}')

    leg = ax.legend(loc=8, numpoints=1, framealpha=0.5,
                    prop={'size': 9}, ncol=3)
    leg.get_frame().set_linewidth(0.0)

    ax.set_xlabel(r'Excitation Potential (eV)', fontsize=14)
    ax.set_ylabel('A(Fe)', fontsize=14)

    ax.set_xlim([fe1_meas_ep.min()-0.2, fe1_meas_ep.max()+0.2])
    ax.set_ylim([fe1_eqw_abunds.min()-0.5, fe1_eqw_abunds.max()+0.5])
    ax.tick_params(axis='both', direction='in', labelsize=12)

    ax.set_title(
        f'Initial Params:  Teff = {teff_i}, logg = {logg_i}, [Fe/H] = {feh_i}, Vmicro = {vmicro_i}', fontsize=10)

    ax = fig.add_subplot(212)

    ax.plot(fe1_norm_eqw, fe1_eqw_abunds, marker='.', color='k',
            ls='None', label=f'Fe I:  logeps = {np.mean(fe1_eqw_abunds):.3f}')
    ax.plot(fe2_norm_eqw, fe2_eqw_abunds, marker='.', color='r',
            ls='None', label=f'Fe II:  logeps = {np.mean(fe2_eqw_abunds):.3f}')

    ax.plot(np.sort(fe1_norm_eqw), line(np.sort(fe1_norm_eqw), *eqw_par), color='k', ls='dashed',
            label=f'ds/dx = {eqw_par[0]:.2f} +/- {np.sqrt(eqw_cov[0,0]):.2f}')

    leg = ax.legend(loc=8, numpoints=1, framealpha=0.5,
                    prop={'size': 9}, ncol=3)
    leg.get_frame().set_linewidth(0.0)

    ax.set_xlabel(r'EW / Wavelength', fontsize=14)
    ax.set_ylabel('A(Fe)', fontsize=14)

    ax.set_xlim([fe1_norm_eqw.min()-0.01, fe1_norm_eqw.max()+0.01])
    ax.set_ylim([fe1_eqw_abunds.min()-0.5, fe1_eqw_abunds.max()+0.5])
    ax.tick_params(axis='both', direction='in', labelsize=12)

    fig.subplots_adjust(hspace=0.2)
    fig.savefig(
        f'{star_name}/{star_name}_t{teff_i:.0f}_g{logg_i:2f}_m{feh_i}_v{vmicro_i}_param_eqw.pdf')
    fig.tight_layout()
    plt.close(fig)
