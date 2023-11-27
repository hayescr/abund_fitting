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


star_list = ['cra2_gdr32672', 'cra2_gdr36912', 'hd122563']

pp = PdfPages(f'cra2_stars_spec_example.pdf')

# Make your figure and the axes that you are going to plot on
fig = plt.figure(figsize = (12,10))

plot_colors = ['k', 'k', 'gray']
plot_styles = ['solid', 'solid', 'solid']
plot_widths = [1.75, 1.75, 1.75]
plot_labels = [f'Cra II GDR3 2672', f'Cra II GDR3 6912', f'HD 122563']
offsets = [2, 1, 0]

xlims = [[4280, 4330], [5155, 5198], [6492, 6500]]

gs = GridSpec(2,3, figure=fig)

for i, xlim in enumerate(xlims):

    if i == 1:
        ax = fig.add_subplot(gs[0, 0:])
    elif i == 0:
        ax = fig.add_subplot(gs[1, 0:2])
    elif i == 2:
        ax = fig.add_subplot(gs[1, 2])

    for star_name, plot_color, plot_style, plot_width, plot_label, offset in zip(star_list, plot_colors, plot_styles, plot_widths, plot_labels, offsets):
        star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
        if i == 0:
            plot_width=1.25
            #ax.text(4281, offset+1.2, plot_label, ha='left', fontsize=16)

        elif i == 1:
            plot_width=1.5
            ax.text(5155.5, offset+1.2, plot_label, ha='left', fontsize=16)
        elif i == 2:
            plot_width=1.75
            #ax.text(6503.5, offset+1.2, plot_label, ha='right', fontsize=16)

        ax.plot(star_data['wavelength'], star_data['flux'] + offset, c=plot_color, ls=plot_style, lw=plot_width, label=f'{plot_label}')

    if i == 1:
        
        line_data = np.genfromtxt('hd122563_lineids_mgb.txt', dtype=None, names=True, encoding=None)
        line_data = line_data[line_data['rel_strength'] > 0.04]


        ax.scatter(line_data['wave'], np.ones(len(line_data))*3.275, marker = '|', s = 400
            , ec='k', linewidths=1.75, ls = 'None', c='k')
        for j, (linewave, lineelem, wave_offset) in enumerate(zip(line_data['wave'], line_data['element'], line_data['wave_offset'])):
            if (linewave > xlim[0]+2) & (linewave < xlim[1]-2):
                flux_offset = 0.0
                if (j!=0):
                    if np.abs(line_data['wave'][j-1]-linewave) < 0.5:
                        flux_offset = 0.0
                if j!=len(line_data)-1:
                    if np.abs(line_data['wave'][j+1]-linewave) < 0.5:
                        flux_offset = 0.2
        
                ax.text(linewave+wave_offset, 3.55+flux_offset, lineelem, fontsize = 14,
                        va='center', ha='center', color='k')

    if i == 2:
        
        line_data = np.genfromtxt('hd122563_lineids_ba6496.txt', dtype=None, names=True, encoding=None)
        line_data = line_data[line_data['rel_strength'] > 0.02]


        ax.scatter(line_data['wave'], np.ones(len(line_data))*3.275, marker = '|', s = 400
            , ec='k', linewidths=1.75, ls = 'None', c='k')
        for j, (linewave, lineelem, wave_offset) in enumerate(zip(line_data['wave'], line_data['element'], line_data['wave_offset'])):
            if (linewave > xlim[0]) & (linewave < xlim[1]):
                flux_offset = 0.0
                if (j!=0):
                    if np.abs(line_data['wave'][j-1]-linewave) < 0.5:
                        flux_offset = 0.0
                if j!=len(line_data)-1:
                    if np.abs(line_data['wave'][j+1]-linewave) < 0.5:
                        flux_offset = 0.2
        
                ax.text(linewave+wave_offset, 3.55+flux_offset, lineelem, fontsize = 14,
                        va='center', ha='center', color='k')

    if i == 0:

        ax.plot(np.array([4282., 4312.]), np.array([3.35, 3.35]), color='k', lw=1.75)
        ax.plot(np.array([4297., 4297.]), np.array([3.35, 3.43]), color='k', lw=1.75)
        ax.text(4297, 3.55, 'CH', fontsize = 14,
                va='center', ha='center', color='k')



        ax.plot(np.array([4323., 4324.5]), np.array([3.35, 3.35]), color='k', lw=1.75)
        ax.plot(np.array([4323.75, 4323.75]), np.array([3.35, 3.43]), color='k', lw=1.75)
        ax.text(4323.75, 3.55, 'CH', fontsize = 14,
                va='center', ha='center', color='k')
        #ax.text(5165, 1.4, r'C$_2$', fontsize = 16,
        #        va='center', ha='center', color='k')

    #ax.plot(np.array([5163., 5165., 5165.]), np.array([1.3, 1.05, 1.3]), color='k', lw=2.)
    #ax.text(5165, 1.4, r'C$_2$', fontsize = 16,
    #        va='center', ha='center', color='k')
        

    #leg = ax.legend(loc=1, numpoints=1, framealpha=0.5, prop={'size':14}, ncol=3)
    #leg.get_frame().set_linewidth(0.0)

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.tick_params(axis='both', direction='in', labelsize=16, which='both')

    ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    ax.yaxis.set_major_locator(MultipleLocator(1))


    ax.set_ylim([-0.2,4.2])

    ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)
    if i == 0:
        ax.set_ylabel('Normalized Flux + Offset', fontsize=16)
    elif i == 1:
        ax.set_ylabel('Normalized Flux + Offset', fontsize=16)
    else:
        ax.set_yticklabels([])

plt.tight_layout()

pp.savefig()
plt.close(fig)


pp.close()
