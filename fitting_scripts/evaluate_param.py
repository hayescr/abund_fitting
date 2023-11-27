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

method = 'chi2'
clip_abunds = False
fix_vmicro = True
fix_logg = False
update_param = True

#params_to_update = ['teff', 'feh', 'vmicro']
#params_to_update = ['logg']
#params_to_update = ['teff', 'logg', 'feh', 'vmicro']
params_to_update = ['feh', 'vmicro']
#params_to_update = ['feh', 'teff', 'vmicro', 'logg']

star_list = ['sgr2_gdr35584','sgr2_gdr32656','sgr2_gdr39936','aqu2_gdr33776','aqu2_gdr31472']

star_name = 'sgr2_gdr35584'
star = star_name


#########
# Reading things in
#########

param_filename = f'{star_name}/{star_name}_params.txt'
star_param = read_json(param_filename)

teff_i = star_param['teff']
logg_i = star_param['logg']
feh_i = star_param['feh']
vmicro_i = star_param['vmicro']
res = star_param['vmacro']


abund_filename = f'{star}/{star}_Fe_abund_fit_convol{res:.1f}.txt'

line_meas = np.genfromtxt(abund_filename, dtype=None, encoding=None, names=True)
fe1_line_data = np.genfromtxt('Fe_I_ts_lines.txt', dtype=None, usecols=(0,1,2), 
                              names=('wavelength', 'ep', 'loggf'), encoding=None)
fe2_line_data = np.genfromtxt('Fe_II_ts_lines.txt', dtype=None, usecols=(0,1,2), 
                              names=('wavelength', 'ep', 'loggf'), encoding=None)



#########
# Cleaning good lines
#########

fe1_lines = np.array(abund_lines['Fe']['I'])
fe2_lines = np.array(abund_lines['Fe']['II'])

good_lines = (line_meas['lim_flag_eqw'].astype('bool')) & (line_meas['lim_flag_chi2'].astype('bool'))
good_lines &= line_meas['avg_synth_stdev']/(1.5/(line_meas['snr']*np.sqrt(line_meas['npix']))) > 1.0
good_lines &= (line_meas['avg_line_depth'] * line_meas['snr'] > 4.)

#a, clip_low, clip_high = sigmaclip(line_meas['wave_shift'], low=3., high=3.)

#good_lines &= (line_meas['wave_shift'] < clip_high) & (line_meas['wave_shift'] > clip_low)


line_meas = line_meas[good_lines]

print(line_meas['wavelength'])

##########
# Matching line data
##########

fe1_meas = line_meas[np.in1d(line_meas['wavelength'], fe1_lines)]
fe2_meas = line_meas[np.in1d(line_meas['wavelength'], fe2_lines)]

fe1_meas_ep = []
fe2_meas_ep = []

for wave in fe1_meas['wavelength']:
    line_match = np.argmin(np.abs(wave-fe1_line_data['wavelength']))
    fe1_meas_ep += [fe1_line_data['ep'][line_match]]
fe1_meas_ep = np.array(fe1_meas_ep)
    
for wave in fe2_meas['wavelength']:
    line_match = np.argmin(np.abs(wave-fe2_line_data['wavelength']))
    fe2_meas_ep += [fe2_line_data['ep'][line_match]]
fe2_meas_ep = np.array(fe2_meas_ep)

##########
# Setting up inputs
##########

    
if method == 'eqw':
    fe1_abunds = fe1_meas['logeps_eqw']
    fe2_abunds = fe2_meas['logeps_eqw']
else:
    fe1_abunds = fe1_meas['logeps_chi2']
    fe2_abunds = fe2_meas['logeps_chi2']
    
fe1_norm_eqw = fe1_meas['eqw']/fe1_meas['wavelength']
fe2_norm_eqw = fe2_meas['eqw']/fe2_meas['wavelength']


if clip_abunds:
    a, low_cut, high_cut = sigmaclip(fe1_abunds, low=3, high=3)
    selection = (fe1_abunds > low_cut) & (fe1_abunds < high_cut)
    fe1_abunds = fe1_abunds[selection]
    fe1_norm_eqw = fe1_norm_eqw
    fe1_meas_ep = fe1_meas_ep[selection]

    #a, low_cut, high_cut = sigmaclip(fe2_abunds, low=3, high=3)
    #selection = (fe2_abunds > low_cut) & (fe2_abunds < high_cut)
    #fe2_abunds = fe2_abunds[selection]
    #fe2_norm_eqw = fe2_norm_eqw
    #fe2_meas_ep = fe2_meas_ep[selection]
    
def line(x, a, b):
    return a*x+b


############
# Fitting relations
############

ep_par, ep_cov = curve_fit(line, fe1_meas_ep, fe1_abunds)
eqw_par, eqw_cov = curve_fit(line, fe1_norm_eqw[fe1_norm_eqw < 0.05], fe1_abunds[fe1_norm_eqw < 0.05])

ion_abund_diff = np.mean(fe1_abunds) - np.mean(fe2_abunds)

fe1_sigma = np.std(fe1_abunds)/np.sqrt((len(fe1_abunds) -1))
fe1_disp = np.std(fe1_abunds)
fe2_sigma = np.std(fe2_abunds)/(len(fe2_abunds) -1)
ion_abund_sigma = np.std(np.concatenate((fe1_abunds,fe2_abunds)))
#ion_abund_sigma = np.hypot(fe1_sigma, fe2_sigma)


#ion_abund_diff = np.median(fe2_abunds) - np.median(fe1_abunds)

#fe1_sigma = median_abs_deviation(fe1_abunds)/(len(fe1_abunds) -1)
#fe2_sigma = median_abs_deviation(fe2_abunds)/(len(fe2_abunds) -1)
#ion_abund_sigma = median_abs_deviation(np.concatenate((fe1_abunds,fe2_abunds)))
#ion_abund_sigma = np.hypot(fe1_sigma, fe2_sigma)

#############
# Checking parameters
#############


# Teff
print()
if (ep_par[0] < np.sqrt(ep_cov[0,0])) and (ep_par[0] > -np.sqrt(ep_cov[0,0])):
    print(f'Teff good : ds/dx = {ep_par[0]} +/- {np.sqrt(ep_cov[0,0])}')
else:
    delta_teff = np.sign(ep_par[0])*np.minimum(np.abs(ep_par[0]/0.0002), 1000.)
    new_teff = np.round(star_param["teff"]+delta_teff, 2)
    print(f'Teff bad : ds/dx = {ep_par[0]:.6f} +/- {np.sqrt(ep_cov[0,0]):.6f}')
    print(f'Delta Teff = {delta_teff:.2f}')
    print(f'Old Teff = {star_param["teff"]:.2f} : New Teff = {new_teff:.2f}')
    if 'teff' in params_to_update:
        print('Teff updated.')
        star_param['teff'] = new_teff
    
# Logg
print()
if (fix_logg):
    print('Using Logg from Stefan-Boltzman')
    new_logg, logg_err = calculate_sb_logg(star, star_param['teff'],
                                           np.sqrt(ep_cov[0, 0]/0.0002))
    print(
        f'Old Logg = {star_param["logg"]:.2f} : New Logg = {new_logg:.2f} : Error = {logg_err:.2f}')
    if 'logg' in params_to_update:
        print('Logg updated')
        star_param['logg'] = new_logg

elif (ion_abund_diff < ion_abund_sigma) and (ion_abund_diff > -ion_abund_sigma):
    print(f'Logg good : delta eps(I - II) = {ion_abund_diff:.4f} +/- {ion_abund_sigma:.4f}')
else:
    delta_logg = np.sign(ion_abund_diff)*np.minimum(np.abs(ion_abund_diff) / 0.7, 1.0)
    new_logg = np.round(star_param["logg"] + delta_logg, 4)
    print(f'Logg bad : delta eps(I - II) = {ion_abund_diff:.4f} +/- {ion_abund_sigma:.4f}')
    print(f'Delta logg = {delta_logg:.2f}')
    print(f'Old logg = {star_param["logg"]:.2f} : New logg = {new_logg:.2f}')
    if 'logg' in params_to_update:
        print('logg updated.')
        star_param['logg'] = new_logg

# [Fe/H]
print()
new_feh = np.round(np.mean(fe1_abunds) - solar_abundances['Fe'], 3)
print(f'Old [Fe/H] = {star_param["feh"]:.2f} : New [Fe/H] = {new_feh:.3f} +/- {fe1_disp:.3f}')
if 'feh' in params_to_update:
    print('[Fe/H] updated.')
    star_param['feh'] = new_feh

# Vmicro
print()
if (fix_vmicro):
    print('Using Vmicro relation from Mashenina et al. 2017')
    new_vmicro = np.round(vmicro_mash(star_param['teff'], star_param['logg'], star_param['feh']), 2)
    print(f'Old Vmicro = {star_param["vmicro"]:.2f} : New Vmicro = {new_vmicro:.2f}')
    star_param['vmicro'] = new_vmicro
elif (eqw_par[0] < np.sqrt(eqw_cov[0,0])) & (eqw_par[0] > -np.sqrt(eqw_cov[0,0])):
    print(f'Vmicro good : ds/deqw = {eqw_par[0]:.6f} +/- {np.sqrt(eqw_cov[0,0]):.6f}')
else:
    delta_vmicro = eqw_par[0]/42.
    new_vmicro = np.round(star_param["vmicro"]+delta_vmicro, 2)
    print(f'Vmicro bad : ds/deqw = {eqw_par[0]:.6f} +/- {np.sqrt(eqw_cov[0,0]):.6f}')
    print(f'Delta Vmicro = {delta_vmicro:.2f}')
    print(f'Old Vmicro = {star_param["vmicro"]:.2f} : New Vmicro = {new_vmicro:.2f}')
    if 'vmicro' in params_to_update:
        print('Vmicro updated.')
        star_param['vmicro'] = new_vmicro
    
if update_param:
    write_json(param_filename, star_param)


############
# Plotting
############

fig = plt.figure(figsize = (8, 7))

ax = fig.add_subplot(211)

ax.plot(fe1_meas_ep, fe1_abunds, marker='.', color='k', ls='None', label=f'Fe I:  logeps = {np.mean(fe1_abunds):.3f}')
ax.plot(fe2_meas_ep, fe2_abunds, marker='.', color='r', ls='None', label=f'Fe II:  logeps = {np.mean(fe2_abunds):.3f}')

ax.plot(np.sort(fe1_meas_ep), line(np.sort(fe1_meas_ep), *ep_par), color = 'k', ls='dashed',
        label=f'ds/dx = {ep_par[0]:.4f} +/- {np.sqrt(ep_cov[0,0]):.4f}')

leg = ax.legend(loc=8, numpoints=1, framealpha=0.5, prop={'size':9}, ncol=3)
leg.get_frame().set_linewidth(0.0)

ax.set_xlabel(r'Excitation Potential (eV)', fontsize=14)
ax.set_ylabel('A(Fe)', fontsize=14)

ax.set_xlim([fe1_meas_ep.min()-0.2, fe1_meas_ep.max()+0.2])
ax.set_ylim([fe1_abunds.min()-0.5, fe1_abunds.max()+0.5])
ax.tick_params(axis='both', direction='in', labelsize=12)

ax.set_title(f'Initial Params:  Teff = {teff_i}, logg = {logg_i}, [Fe/H] = {feh_i}, Vmicro = {vmicro_i}', fontsize=10)


ax = fig.add_subplot(212)

ax.plot(fe1_norm_eqw, fe1_abunds, marker='.', color='k', ls='None', label=f'Fe I:  logeps = {np.mean(fe1_abunds):.3f}')
ax.plot(fe2_norm_eqw, fe2_abunds, marker='.', color='r', ls='None', label=f'Fe II:  logeps = {np.mean(fe2_abunds):.3f}')

ax.plot(np.sort(fe1_norm_eqw), line(np.sort(fe1_norm_eqw), *eqw_par), color = 'k', ls='dashed',
        label=f'ds/dx = {eqw_par[0]:.2f} +/- {eqw_cov[0,0]:.2f}')

leg = ax.legend(loc=8, numpoints=1, framealpha=0.5, prop={'size':9}, ncol=3)
leg.get_frame().set_linewidth(0.0)

ax.set_xlabel(r'EW / Wavelength', fontsize=14)
ax.set_ylabel('A(Fe)', fontsize=14)

ax.set_xlim([fe1_norm_eqw.min()-0.01, fe1_norm_eqw.max()+0.01])
ax.set_ylim([fe1_abunds.min()-0.5, fe1_abunds.max()+0.5])
ax.tick_params(axis='both', direction='in', labelsize=12)

fig.subplots_adjust(hspace=0.2)
fig.savefig(f'{star_name}/{star_name}_t{teff_i:.0f}_g{logg_i:2f}_m{feh_i}_v{vmicro_i}_param_check.pdf')
fig.tight_layout()
plt.close(fig)



# Need to reject any remaining lines with outlying abundance measurements and figure which version to use
# Measure gradients as a function of ep, differences between ionized and not, and as a funciton of eqw/wave or wave
# Evaluate if gradients are sufficiently small
# If not calculate how much to change the parameters
# update parameters
# Caveats
#     How do I do eqw, using the eqw from the file might include flux from blends and give a bit of a skewed view
#     does that even matter if I stick to the micro relation from lit?
    