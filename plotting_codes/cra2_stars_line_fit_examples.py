import numpy as np
from linelist_manager import LinelistManager
from abundometer import Abundometer
from turbospec_wrapper.turbosynth import TurboSynth, iso_frac2
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
from scipy.stats import sigmaclip
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.gridspec import GridSpec



nbins = 100
cmap_name = 'hayescr'
synth_colors = ['red', 'orange', 'gold', 'cyan', 'dodgerblue', 'mediumorchid']
cmap = colors.LinearSegmentedColormap.from_list(cmap_name, synth_colors, N=nbins)


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
if not os.path.isdir('tmp2'):
    os.mkdir('tmp2')

pp = PdfPages(f'cra2_stars_fits_example.pdf')

fig = plt.figure(figsize = (16,15))

gs = GridSpec(3,6, figure=fig)

ax = fig.add_subplot(gs[:,:])

ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)

ax.set_ylabel('Normalized Flux + Offset', fontsize=16)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='None', top=False, bottom=False, left=False, right=False)



#star_list = ['hd122563', 'ret2_mem', 'ret2_cand']
#star_list = ['ret2_cand']
#star_list = ['hd122563']
#star_list = ['ret2_mem', 'ret2_cand']
star_list = ['ret2_test']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']

wavelength_limit = 4000.

update_abund = False

element = 'C'
lines = np.array(abund_lines[element]['CH'])
lines = lines[lines > wavelength_limit]
synth_range = 40

add_vel_shift=False
run_selrange=False

solar_abundances = au.solar_abund()

final_abund = np.genfromtxt('cra2_abund_error_table_ions.txt', encoding=None, dtype=None, names=True)

for star_name in star_list:

    logeps_c = final_abund['CCH_FE'] + final_abund['FEI_H'] + solar_abundances['C']
    final_c_fe = final_abund['CCH_FE'][final_abund['STAR']==star_name][0]
    final_c = logeps_c[final_abund['STAR']==star_name][0]

    star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
    mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
    star = f'{star_name}'
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)


    if True:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '25Mg1H', '26Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '12C12C', '28Si1H', 
                   '56Fe1H', '40Ca1H', '12C13C', '13C14N']
    else:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '28Si1H', 
                   '56Fe1H', '40Ca1H']


    ll_helper = LinelistManager(molecules=mp_mols, turbospec_path='../../Turbospectrum2019', source='linemake')

    ex_turbosynth2 = initiate_turbosynth(star_param['teff'], 
                                         star_param['logg'],
                                         star_param['feh'],
                                         star_param['vmicro'],
                                         star_param['cfe'],
                                         star_param['nfe'],
                                         star_param['alphafe'],
                                         star_param['abundances'])
    #ex_turbosynth2.set_isotope({12.024: 0.58, 12.025: 0.22, 12.026: 0.22})
    #c12c13 = 5
    #ex_turbosynth2.set_isotope({6.012: iso_frac2('major', c12c13),
    #                             6.013: iso_frac2('minor', c12c13)})


    res = star_param['vmacro']

    if not os.path.isdir(star):
        os.mkdir(star)

    conv_profiles = [2]
    conv_broads = [{'vel': res}]


    for plot_i, line in enumerate(lines):
        if not os.path.exists(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt'):
            ll_range = [line-15, line+15]
            linelists = ll_helper.build_library(ll_range)
            ex_turbosynth2.set_linelists(linelists)

            abund_fitter = Abundometer(ex_turbosynth2, mem, element, line,
                                       conv_profiles, conv_broads,
                                       synth_range=synth_range,
                                    )
            abund_fitter.fit(add_shift=add_vel_shift, run_selrange=run_selrange)

            abund_fitter.turbosynth_obj.set_abund('C', final_c)

            abund_fitter.synthesize(delta_abunds=[-5., -0.3, 0.0, 0.3])

            obs_spec = abund_fitter.obs_spec
            synths = [synth for synth in abund_fitter.synths]

            with open(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', 'w') as output_file:
                print('# wavelength flux_obs flux_0 flux_1 flux_2 flux_3', file=output_file)
                for wave, obs, s1, s2, s3, s4 in zip(obs_spec.wavelength, obs_spec.flux, synths[0].flux, synths[1].flux, synths[2].flux, synths[3].flux):
                    print(wave, obs, s1, s2, s3, s4, file=output_file)



for plot_i, line in enumerate(lines):

    ax = fig.add_subplot(gs[0, 0:])

    for star_num, star_name in enumerate(star_list):

        plot_offset = 1 - star_num

        logeps_c = final_abund['CCH_FE'] + final_abund['FEI_H'] + solar_abundances['C']
        final_c_fe = final_abund['CCH_FE'][final_abund['STAR']==star_name][0]
        final_c = logeps_c[final_abund['STAR']==star_name][0]

        star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
        mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
        star = f'{star_name}'
        param_filename = f'{star_name}/{star_name}_params.txt'
        star_param = read_json(param_filename)


        c_fit_data = np.genfromtxt(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', dtype=None, encoding=None, names=True)

        obs_spec = Spectrum(c_fit_data['wavelength'], c_fit_data['flux_obs'])
        synths = [Spectrum(c_fit_data['wavelength'], c_fit_data[f'flux_{spec_index}']) for spec_index in range(4)]


        # Make your figure and the axes that you are going to plot on
        #fig = plt.figure(figsize = (16,7))

        if star_num == 0:
            star_label = 'Obs'
        else:
            star_label = None

        ax.plot(obs_spec.wavelength, obs_spec.flux+plot_offset, c='k', lw=1.5, label=star_label)

        c_abund = star_param['abundances']['C']

        print(c_abund)

        plot_colors = ['b', 'None', 'r', 'None']
        plot_styles = ['solid', 'dashed', 'solid', 'dashed']
        plot_widths = [1.0, 1.0, 1.5, 1.0]
        plot_labels = [f'No C', f'[C/Fe] = {final_c_fe-0.3:.2f}', f'[C/Fe] = {final_c_fe:.2f}', f'[C/Fe] = {final_c_fe+0.3:.2f}']

        plot_labels = [f'No C', None, f'best fit C', None]

        for synth, plot_color, plot_style, plot_width, plot_label in zip(synths, plot_colors, plot_styles, plot_widths, plot_labels):
            if (plot_i != 0) or (star_num != 0):
                plot_label= None
            if plot_color == 'None':
                continue
            ax.plot(synth.wavelength, synth.flux+plot_offset, c=plot_color, ls=plot_style, lw=plot_width, label=plot_label)
            #if plot_color == 'b':
            #    ax.plot([0], [0], ls='None', c='None', label=' ')

        upsynth = synths[3]
        losynth = synths[1]

        ax.fill_between(losynth.wavelength, losynth.flux+plot_offset, upsynth.flux+plot_offset, fc='r', alpha=0.4)

        if star_num == 0:
            star_name_label = 'GDR3 2672'
            label_loc = 2.2
        elif star_num == 1:
            star_name_label = 'GDR3 6912'
            label_loc = 1.1

        ax.text(4300.20, label_loc, f'{star_name_label}, [C/Fe] = {final_c_fe:.2f}', ha='left', fontsize=12)

    #ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)

    #if plot_i == 0:
    #    ax.set_ylabel('Normalized Flux + Offset', fontsize=16)
    #else:
    #    ax.set_yticklabels([])

    if plot_i == 0:
       leg = ax.legend(loc=4, numpoints=1, framealpha=0.5, prop={'size':12}, ncol=3)
       leg.get_frame().set_linewidth(0.0)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.tick_params(axis='both', direction='in', labelsize=16, which='both')

    ax.set_xlim([line -10, line+16])
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        #x.xaxis.set_minor_locator(MultipleLocator(1))

    ax.set_ylim([-0.2,2.5])

    ax.text(line+15.65, 2.2, f'CH', ha='right', fontsize=24)












#star_list = ['hd122563', 'ret2_mem', 'ret2_cand']
#star_list = ['ret2_cand']
#star_list = ['hd122563']
#star_list = ['ret2_mem', 'ret2_cand']
star_list = ['ret2_test']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']

wavelength_limit = 4400.

update_abund = False

element = 'Ti'
lines = np.array(abund_lines[element]['II'])
lines = np.array([5129.16, 5188.69, 5226.54])
lines = np.array([4470.85, 5188.69, 5336.79])
lines = np.array([5129.16, 5188.69, 5336.79])
#lines = np.array([5188.69, 5226.54, 5336.79])
synth_range = 5



add_vel_shift=False
run_selrange=False

solar_abundances = au.solar_abund()

final_abund = np.genfromtxt('cra2_abund_error_table_ions.txt', encoding=None, dtype=None, names=True)

for star_name in star_list:

    logeps_ti = final_abund['TIII_FE'] + final_abund['FEI_H'] + solar_abundances[element]
    final_ti_fe = final_abund['TIII_FE'][final_abund['STAR']==star_name][0]
    final_ti = logeps_ti[final_abund['STAR']==star_name][0]

    star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
    mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
    star = f'{star_name}'
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)


    if True:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '25Mg1H', '26Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '12C12C', '28Si1H', 
                   '56Fe1H', '40Ca1H', '12C13C', '13C14N']
    else:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '28Si1H', 
                   '56Fe1H', '40Ca1H']


    ll_helper = LinelistManager(molecules=mp_mols, turbospec_path='../../Turbospectrum2019', source='linemake')

    ex_turbosynth2 = initiate_turbosynth(star_param['teff'], 
                                         star_param['logg'],
                                         star_param['feh'],
                                         star_param['vmicro'],
                                         star_param['cfe'],
                                         star_param['nfe'],
                                         star_param['alphafe'],
                                         star_param['abundances'])
    #ex_turbosynth2.set_isotope({12.024: 0.58, 12.025: 0.22, 12.026: 0.22})
    #c12c13 = 5
    #ex_turbosynth2.set_isotope({6.012: iso_frac2('major', c12c13),
    #                             6.013: iso_frac2('minor', c12c13)})


    res = star_param['vmacro']

    if not os.path.isdir(star):
        os.mkdir(star)

    conv_profiles = [2]
    conv_broads = [{'vel': res}]


    for plot_i, line in enumerate(lines):
        if not os.path.exists(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt'):
            ll_range = [line-15, line+15]
            linelists = ll_helper.build_library(ll_range)
            ex_turbosynth2.set_linelists(linelists)

            abund_fitter = Abundometer(ex_turbosynth2, mem, element, line,
                                       conv_profiles, conv_broads,
                                       synth_range=synth_range,
                                    )
            abund_fitter.fit(add_shift=add_vel_shift, run_selrange=run_selrange)

            abund_fitter.turbosynth_obj.set_abund(element, final_ti)

            abund_fitter.synthesize(delta_abunds=[-5., -0.3, 0.0, 0.3])

            obs_spec = abund_fitter.obs_spec
            synths = [synth for synth in abund_fitter.synths]

            with open(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', 'w') as output_file:
                print('# wavelength flux_obs flux_0 flux_1 flux_2 flux_3', file=output_file)
                for wave, obs, s1, s2, s3, s4 in zip(obs_spec.wavelength, obs_spec.flux, synths[0].flux, synths[1].flux, synths[2].flux, synths[3].flux):
                    print(wave, obs, s1, s2, s3, s4, file=output_file)





for plot_i, line in enumerate(lines):

    ax = fig.add_subplot(gs[1,2*plot_i:2*(plot_i+1)])

    for star_num, star_name in enumerate(star_list):

        plot_offset = 1-star_num

        logeps_ti = final_abund['TIII_FE'] + final_abund['FEI_H'] + solar_abundances[element]
        final_ti_fe = final_abund['TIII_FE'][final_abund['STAR']==star_name][0]
        final_ti = logeps_ti[final_abund['STAR']==star_name][0]

        star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
        mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
        star = f'{star_name}'
        param_filename = f'{star_name}/{star_name}_params.txt'
        star_param = read_json(param_filename)


        c_fit_data = np.genfromtxt(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', dtype=None, encoding=None, names=True)

        obs_spec = Spectrum(c_fit_data['wavelength'], c_fit_data['flux_obs'])
        synths = [Spectrum(c_fit_data['wavelength'], c_fit_data[f'flux_{spec_index}']) for spec_index in range(4)]


        # Make your figure and the axes that you are going to plot on
        #fig = plt.figure(figsize = (16,7))

        if star_num == 0:
            star_label = 'Obs'
        else:
            star_label = None

        ax.plot(obs_spec.wavelength, obs_spec.flux+plot_offset, c='k', lw=2.0, label=star_label)


        plot_colors = ['b', 'None', 'r', 'None']
        plot_styles = ['solid', 'dashed', 'solid', 'dashed']
        plot_widths = [1.0, 1.0, 1.5, 1.0]
        plot_labels = [f'No {element}', f'[{element}/Fe] = {final_ti_fe-0.3:.2f}', f'[{element}/Fe] = {final_ti_fe:.2f}', f'[element/Fe] = {final_ti_fe+0.3:.2f}']

        plot_labels = [f'No {element}', None, f'best fit {element}', None]

        for synth, plot_color, plot_style, plot_width, plot_label in zip(synths, plot_colors, plot_styles, plot_widths, plot_labels):
            if (plot_i != 2) or (star_num != 1):
                plot_label= None
            if plot_color == 'None':
                continue
            ax.plot(synth.wavelength, synth.flux+plot_offset, c=plot_color, ls=plot_style, lw=plot_width, label=plot_label)
            #if plot_color == 'b':
            #    ax.plot([0], [0], ls='None', c='None', label=' ')

        upsynth = synths[3]
        losynth = synths[1]

        ax.fill_between(losynth.wavelength, losynth.flux+plot_offset, upsynth.flux+plot_offset, fc='r', alpha=0.4)

        if star_num == 0:
            star_name_label = 'GDR3 2672'
            label_loc = 2.2
        elif star_num == 1:
            star_name_label = 'GDR3 6912'
            label_loc = 1.10

        if plot_i == 0:
            ax.text(5127.25, label_loc, f'{star_name_label}, [{element}/Fe] = {final_ti_fe:.2f}', ha='left', fontsize=12)

    #ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)

    #if plot_i == 0:
    #    ax.set_ylabel('Normalized Flux + Offset', fontsize=16)
    #else:
    #    ax.set_yticklabels([])

    if plot_i == 2:
       leg = ax.legend(loc=4, numpoints=1, framealpha=0.5, prop={'size':12}, ncol=1)
       leg.get_frame().set_linewidth(0.0)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.tick_params(axis='both', direction='in', labelsize=16, which='both')

    ax.set_xlim([line -2, line+2])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))

    ax.set_ylim([-0.2,2.5])

    ax.text(line+1.85, 2.2, f'Ti II', ha='right', fontsize=24)












#star_list = ['hd122563', 'ret2_mem', 'ret2_cand']
#star_list = ['ret2_cand']
#star_list = ['hd122563']
#star_list = ['ret2_mem', 'ret2_cand']
star_list = ['ret2_test']
star_list = ['cra2_gdr32672', 'cra2_gdr36912']

wavelength_limit = 4400.

update_abund = False

element = 'Ni'
lines = np.array(abund_lines[element]['I'])
lines = np.array([6643.60, 6767.80])
synth_range = 5


add_vel_shift=False
run_selrange=False

solar_abundances = au.solar_abund()

final_abund = np.genfromtxt('cra2_abund_error_table_ions.txt', encoding=None, dtype=None, names=True)

for star_name in star_list:

    logeps_ni = final_abund['NII_FE'] + final_abund['FEI_H'] + solar_abundances['Ni']
    final_ni_fe = final_abund['NII_FE'][final_abund['STAR']==star_name][0]
    final_ni = logeps_ni[final_abund['STAR']==star_name][0]

    star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
    mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
    star = f'{star_name}'
    param_filename = f'{star_name}/{star_name}_params.txt'
    star_param = read_json(param_filename)


    if True:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '25Mg1H', '26Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '12C12C', '28Si1H', 
                   '56Fe1H', '40Ca1H', '12C13C', '13C14N']
    else:
        mp_mols = ['12C1H', '13C1H', 
                   '24Mg1H', '16O1H', 
                   '12C14N', '14N1H', 
                   '28Si1H', 
                   '56Fe1H', '40Ca1H']


    ll_helper = LinelistManager(molecules=mp_mols, turbospec_path='../../Turbospectrum2019', source='linemake')

    ex_turbosynth2 = initiate_turbosynth(star_param['teff'], 
                                         star_param['logg'],
                                         star_param['feh'],
                                         star_param['vmicro'],
                                         star_param['cfe'],
                                         star_param['nfe'],
                                         star_param['alphafe'],
                                         star_param['abundances'])
    #ex_turbosynth2.set_isotope({12.024: 0.58, 12.025: 0.22, 12.026: 0.22})
    #c12c13 = 5
    #ex_turbosynth2.set_isotope({6.012: iso_frac2('major', c12c13),
    #                             6.013: iso_frac2('minor', c12c13)})


    res = star_param['vmacro']

    if not os.path.isdir(star):
        os.mkdir(star)

    conv_profiles = [2]
    conv_broads = [{'vel': res}]


    for plot_i, line in enumerate(lines):
        if not os.path.exists(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt'):
            ll_range = [line-15, line+15]
            linelists = ll_helper.build_library(ll_range)
            ex_turbosynth2.set_linelists(linelists)

            abund_fitter = Abundometer(ex_turbosynth2, mem, element, line,
                                       conv_profiles, conv_broads,
                                       synth_range=synth_range,
                                    )
            abund_fitter.fit(add_shift=add_vel_shift, run_selrange=run_selrange)

            abund_fitter.turbosynth_obj.set_abund('Ni', final_ni)

            abund_fitter.synthesize(delta_abunds=[-5., -0.3, 0.0, 0.3])

            obs_spec = abund_fitter.obs_spec
            synths = [synth for synth in abund_fitter.synths]

            with open(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', 'w') as output_file:
                print('# wavelength flux_obs flux_0 flux_1 flux_2 flux_3', file=output_file)
                for wave, obs, s1, s2, s3, s4 in zip(obs_spec.wavelength, obs_spec.flux, synths[0].flux, synths[1].flux, synths[2].flux, synths[3].flux):
                    print(wave, obs, s1, s2, s3, s4, file=output_file)








for plot_i, line in enumerate(lines):

    ax = fig.add_subplot(gs[2,1+2*plot_i:1+2*(plot_i+1)])

    for star_num, star_name in enumerate(star_list):

        plot_offset = 1 - star_num

        logeps_ni = final_abund['NII_FE'] + final_abund['FEI_H'] + solar_abundances['Ni']
        final_ni_fe = final_abund['NII_FE'][final_abund['STAR']==star_name][0]
        final_ni = logeps_ni[final_abund['STAR']==star_name][0]

        star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
        mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
        star = f'{star_name}'
        param_filename = f'{star_name}/{star_name}_params.txt'
        star_param = read_json(param_filename)


        c_fit_data = np.genfromtxt(f'{star}/abund_no_velshift/{star}_{element}_{line}_synth.txt', dtype=None, encoding=None, names=True)

        obs_spec = Spectrum(c_fit_data['wavelength'], c_fit_data['flux_obs'])
        synths = [Spectrum(c_fit_data['wavelength'], c_fit_data[f'flux_{spec_index}']) for spec_index in range(4)]


        # Make your figure and the axes that you are going to plot on
        #fig = plt.figure(figsize = (16,7))

        if star_num == 0:
            star_label = 'Obs'
        else:
            star_label = None

        ax.plot(obs_spec.wavelength, obs_spec.flux+plot_offset, c='k', lw=2.0, label=star_label)


        plot_colors = ['b', 'None', 'r', 'None']
        plot_styles = ['solid', 'dashed', 'solid', 'dashed']
        plot_widths = [1.0, 1.0, 1.5, 1.0]
        plot_labels = [f'No Ni', f'[Ni/Fe] = {final_ni_fe-0.3:.2f}', f'[Ni/Fe] = {final_ni_fe:.2f}', f'[Ni/Fe] = {final_ni_fe+0.3:.2f}']

        plot_labels = [f'No Ni', None, f'best fit Ni', None]

        for synth, plot_color, plot_style, plot_width, plot_label in zip(synths, plot_colors, plot_styles, plot_widths, plot_labels):
            if (plot_i != 1) or (star_num != 1):
                plot_label= None
            if plot_color == 'None':
                continue
            ax.plot(synth.wavelength, synth.flux+plot_offset, c=plot_color, ls=plot_style, lw=plot_width, label=plot_label)
            #if plot_color == 'b':
            #    ax.plot([0], [0], ls='None', c='None', label=' ')

        upsynth = synths[3]
        losynth = synths[1]

        ax.fill_between(losynth.wavelength, losynth.flux+plot_offset, upsynth.flux+plot_offset, fc='r', alpha=0.4)

        if star_num == 0:
            star_name_label = 'GDR3 2672'
            label_loc = 2.15
        elif star_num == 1:
            star_name_label = 'GDR3 6912'
            label_loc = 1.1

        if plot_i == 0:
            ax.text(6641.70, label_loc, f'{star_name_label}, [Ni/Fe] = {final_ni_fe:.2f}', ha='left', fontsize=12)

    #ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)

    #if plot_i == 0:
    #    ax.set_ylabel('Normalized Flux + Offset', fontsize=16)
    #else:
    #    ax.set_yticklabels([])

    if plot_i == 1:
       leg = ax.legend(loc=4, numpoints=1, framealpha=0.5, prop={'size':12}, ncol=1)
       leg.get_frame().set_linewidth(0.0)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.tick_params(axis='both', direction='in', labelsize=16, which='both')

    ax.set_xlim([line -2, line+2])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))

    ax.set_ylim([-0.2,2.5])

    ax.text(line+1.85, 2.2, f'Ni I', ha='right', fontsize=24)



plt.tight_layout()



pp.savefig()
plt.close(fig)

pp.close()