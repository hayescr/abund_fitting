import numpy as np
from linelist_manager import LinelistManager
from abundometer import Abundometer
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
from scipy.stats import sigmaclip


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

#star_list = ['hd122563', 'ret2_mem', 'ret2_cand']
#star_list = ['ret2_cand']
#star_list = ['ret2_mem', 'ret2_cand']
#star_list = ['hd122563', 'ret2_mem']
star_list = ['sgr2_gdr35584','sgr2_gdr32656','sgr2_gdr39936','aqu2_gdr33776','aqu2_gdr31472']

wavelength_limit = 4500.

update_abund = False

element_list = ['Mg']

for element in element_list:
    
    ions = [key for key in abund_lines[element] if key!='convol']

    lines=[]
    for ion in ions:
        lines += abund_lines[element][ion]
    lines = np.array(lines)
    lines = lines[lines > wavelength_limit]

    synth_range = 5

    add_vel_shift=False

    run_selrange=True
    iterations = 5

    solar_abundances = au.solar_abund()


    for star_name in star_list:

        star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
        mem = Spectrum(star_data['wavelength'], star_data['flux'], variance=star_data['variance'])
        star = f'{star_name}'
        param_filename = f'{star_name}/{star_name}_params.txt'
        star_param = read_json(param_filename)
#         error_file = np.genfromtxt('star_param_errors.txt', dtype=None, encoding=None, names=True)

        if True:
            mp_mols = ['12C1H', '13C1H', 
                       '24Mg1H', '25Mg1H', '26Mg1H', '16O1H', 
                       '12C14N', '14N1H', 
                       '12C12C', '28Si1H', 
                       '56Fe1H', '40Ca1H']
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
        #ex_turbosynth2.set_abund('Dy', ex_turbosynth2.get_abund('Dy')+1.72)
        #ex_turbosynth2.set_abund('Eu', ex_turbosynth2.get_abund('Eu')+1.)

        if star_name=='ret2_cand':
            heavy_alphas = ['Si', 'S', 'Ar', 'Ca', 'Ti']
            for alpha in heavy_alphas:
                if (alpha not in star_param['abundances']):
                    ex_turbosynth2.set_abund(alpha, star_param['feh'] + 0.4 + solar_abundances[alpha])

        res = star_param['vmacro']

        if not os.path.isdir(star):
            os.mkdir(star)

        conv_profiles = [2]
        conv_broads = [{'vel': res}]

        if not os.path.isdir(f'{star}/abund_no_velshift'):
            os.mkdir(f'{star}/abund_no_velshift')


        outfile = open(f'{star}/abund_no_velshift/{star}_{element}_abund_fit_convol{res:.1f}.txt', 'w')
        pp = PdfPages(f'{star}/abund_no_velshift/{star}_{element}_abund_fit_plots_convol{res:.1f}.pdf')

        print('# wavelength snr npix wave_shift avg_line_depth avg_synth_stdev eqw eqw_err',
            'logeps_eqw logeps_err_eqw logeps_lim_eqw lim_flag_eqw',
            'chi2 logeps_chi2 lim_flag_chi2', file=outfile)


        for line in lines:

            try:
                ll_range = [line-15, line+15]
                linelists = ll_helper.build_library(ll_range)
                ex_turbosynth2.set_linelists(linelists)

                abund_fitter = Abundometer(ex_turbosynth2, mem, element, line,
                                           conv_profiles, conv_broads,
                                           synth_range=synth_range,
                                        )
                abund_fitter.fit(add_shift=add_vel_shift, run_selrange=run_selrange)

                for i in range(iterations-1):
                    new_delta_abund = abund_fitter.fit_abund_chi2 - ex_turbosynth2.get_abund(element)
                    if (new_delta_abund <= abund_fitter.delta_abunds[1]) | (
                        new_delta_abund >= abund_fitter.delta_abunds[-2]):

                        abund_fitter = Abundometer(ex_turbosynth2, mem, element, line,
                                                   conv_profiles, conv_broads,
                                                   synth_range=synth_range,
                                                   delta_abunds=new_delta_abund + np.arange(-0.6, 0.61, 0.3),
                                       )
                        abund_fitter.fit(add_shift=add_vel_shift, run_selrange=run_selrange)



                abund_fitter.avg_synth.wavelength, abund_fitter.avg_synth.variance

                synth_sel = (abund_fitter.avg_synth.wavelength > abund_fitter.fit_range[0]) & (
                             abund_fitter.avg_synth.wavelength < abund_fitter.fit_range[1])

                synth_elem_sens = np.mean(abund_fitter.avg_synth.variance[synth_sel])
                fit_range_pix = len(abund_fitter.avg_synth.variance[synth_sel])

                new_delta_abund = abund_fitter.fit_abund_chi2 - ex_turbosynth2.get_abund(element)

                synth_limit_flag = (new_delta_abund > (abund_fitter.delta_abunds[0] + abund_fitter.delta_abunds[1])/2.) & (new_delta_abund < (abund_fitter.delta_abunds[-2] + abund_fitter.delta_abunds[-1])/2.)

                print(f'{line:.2f} {abund_fitter.snr:.2f} {fit_range_pix} {abund_fitter.shift:.4f}',
                      f'{1-abund_fitter.pseudo_line_depth_obs:.4f} {synth_elem_sens:.4f}',
                      f'{abund_fitter.eqw:.1f} {abund_fitter.eqw_err:.1f}',
                      f'{abund_fitter.fit_abund_eqw:.4f} {abund_fitter.fit_abund_eqw_err:.4f}',
                      f'{abund_fitter.eqw_limit:.4f} {abund_fitter.eqw_lim_flag}',
                      f'{abund_fitter.fit_min_chi2:.4f} {abund_fitter.fit_abund_chi2:.4f}',
                      f'{int(synth_limit_flag)}', file=outfile)

                cNorm  = colors.Normalize(vmin=abund_fitter.synth_abunds[0], vmax=abund_fitter.synth_abunds[-1])
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)


                # Make your figure and the axes that you are going to plot on
                fig = plt.figure(figsize = (12,7))

                ax = fig.add_subplot(111)

                ax.plot(abund_fitter.obs_spec.wavelength, abund_fitter.obs_spec.flux, c='k', label='Obs')

                for synth, delta in zip(abund_fitter.synths, abund_fitter.synth_abunds):
                    synth_color = scalarMap.to_rgba(delta)
                    ax.plot(synth.wavelength, synth.flux, c=synth_color, label=f'A({abund_fitter.element}) = {delta:.1f}')

                ax.plot(abund_fitter.continuum.wavelength, abund_fitter.continuum.flux, c='r', marker='.', ms=5, ls='None')
                ax.plot(abund_fitter.avg_synth.wavelength, abund_fitter.avg_synth.variance, c='gray')

                ax.plot(abund_fitter.fit_range, [-0.5, -0.5], c='gray', lw=5)

                ax.text(0.02, 0.96, f'log ϵ = {abund_fitter.fit_abund_chi2:.3f}', ha='left', va='top', 
                        transform=ax.transAxes, size=16)

                ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)
                ax.set_ylabel('Normalized Flux', fontsize=16)

                leg = ax.legend(loc=1, numpoints=1, framealpha=0.5, prop={'size':14}, ncol=2)
                leg.get_frame().set_linewidth(0.0)

                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                ax.get_xaxis().get_major_formatter().set_scientific(False)


                ax.set_xlim([abund_fitter.wave_range[0]+0.15, abund_fitter.wave_range[1]-0.15])







                ax2 = fig.add_axes([0.7, 0.2, 0.15, 0.15]) #left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]

                ax2.plot(abund_fitter.synth_abunds, abund_fitter.synth_chi2s, marker='o', c='r', ls='None')

                new_abunds = np.arange(abund_fitter.synth_abunds.min(),
                                       abund_fitter.synth_abunds.max()+0.0001, 0.0001)

                chi2_interp = interp1d(abund_fitter.synth_abunds, abund_fitter.synth_chi2s, kind='cubic', fill_value='extrapolate')
                new_chi2 = chi2_interp(new_abunds)

                ax2.plot(new_abunds, new_chi2, c='r', ls='dashed')
                ax2.set_xlabel(r'$\rm \log{\rm \ ϵ}$', fontsize=12)
                ax2.set_ylabel(r'$\rm \chi^2$', fontsize=12)

                ax.set_ylim([-1,2])

                pp.savefig()
                plt.close(fig)

            except:
                print(line, -999., 999, 0., -999., -999., -999., -999., -99., -99., -99., 0, -999., -99., 0, file=outfile)



        outfile.close()
        pp.close()

        if update_abund:
            abund_filename = f'{star}/{star}_{element}_abund_fit_convol{res:.1f}.txt'
            abund_file = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)

            good_lines = (abund_file['lim_flag_eqw'].astype('bool')) & (abund_file['lim_flag_chi2'].astype('bool'))
            good_lines &= abund_file['avg_synth_stdev']/(2/(abund_file['snr']*np.sqrt(abund_file['npix']))) > 1.0
            good_lines &= (abund_file['avg_line_depth'] * abund_file['snr'] > 5.)

            if add_vel_shift:

                a, clip_low, clip_high = sigmaclip(abund_file['wave_shift'], low=3., high=3.)

                good_lines &= (abund_file['wave_shift'] < clip_high) & (abund_file['wave_shift'] > clip_low)


            if len(abund_file[good_lines]) > 0:

                star_param['abundances'][element] = np.mean(abund_file['logeps_chi2'][good_lines])
                write_json(param_filename, star_param)