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


star_list = ['sgr2_gdr35584','sgr2_gdr32656','sgr2_gdr39936','aqu2_gdr33776','aqu2_gdr31472', 'hd122563']

pp = PdfPages(f'sgr2aqu2_stars_mgbtriplet.pdf')

# Make your figure and the axes that you are going to plot on
fig = plt.figure(figsize = (16,7))

ax = fig.add_subplot(111)

plot_colors = ['k','k','k','k', 'k', 'gray']
plot_styles = ['solid', 'solid', 'solid', 'solid', 'solid', 'solid']
plot_widths = [2., 2., 2., 2., 2., 2.]
plot_labels = ['sgr2_gdr35584','sgr2_gdr32656','sgr2_gdr39936','aqu2_gdr33776','aqu2_gdr31472', f'HD 122563']
offsets = [5, 4, 3, 2, 1, 0]

for star_name, plot_color, plot_style, plot_width, plot_label, offset in zip(star_list, plot_colors, plot_styles, plot_widths, plot_labels, offsets):
    star_data = np.genfromtxt(f'{star_name}_spec_normrv.txt', dtype=None, encoding=None, names=True)
    ax.plot(star_data['wavelength'], star_data['flux'] + offset, c=plot_color, ls=plot_style, lw=plot_width, label=f'{plot_label}')
    ax.text(5140.5, offset+1.2, plot_label, ha='left', fontsize=16)

line_data = np.genfromtxt('hd122563_lineids_mgb.txt', dtype=None, names=True, encoding=None)
line_data = line_data[line_data['rel_strength'] > 0.04]


ax.scatter(line_data['wave'], np.ones(len(line_data))*6.3, marker = '|', s = 480
    , ec='k', linewidths=2., ls = 'None', c='k')
for i, (linewave, lineelem, wave_offset) in enumerate(zip(line_data['wave'], line_data['element'], line_data['wave_offset'])):
    if (linewave > 5148) & (linewave < 5198):
        flux_offset = 0.0
        if (i!=0):
            if np.abs(line_data['wave'][i-1]-linewave) < 0.5:
                flux_offset = 0.0
        if i!=len(line_data)-1:
            if np.abs(line_data['wave'][i+1]-linewave) < 0.5:
                flux_offset = 0.2

        ax.text(linewave+wave_offset, 6.55+flux_offset, lineelem, fontsize = 16,
                va='center', ha='center', color='k')

#ax.plot(np.array([5163., 5165., 5165.]), np.array([1.3, 1.05, 1.3]), color='k', lw=2.)
#ax.text(5165, 1.4, r'C$_2$', fontsize = 16,
#        va='center', ha='center', color='k')

ax.set_xlabel(r'Wavelength $(\rm \AA)$', fontsize=16)
ax.set_ylabel('Normalized Flux + Offset', fontsize=16)

#leg = ax.legend(loc=1, numpoints=1, framealpha=0.5, prop={'size':14}, ncol=3)
#leg.get_frame().set_linewidth(0.0)

ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_xaxis().get_major_formatter().set_scientific(False)

ax.tick_params(axis='both', direction='in', labelsize=16, which='both')

ax.set_xlim([5140, 5198])
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.set_ylim([-0.2,7.2])

plt.tight_layout()

pp.savefig()
plt.close(fig)


pp.close()
