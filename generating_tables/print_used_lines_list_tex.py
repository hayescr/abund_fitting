import numpy as np
from abund_lines import *
import abund_utils as au

table = np.genfromtxt('cra2_abund_used_line_measurements.txt', encoding=None, dtype=None, names=True)

inspect = np.genfromtxt('cra2_abund_used_line_measurements_inspect_nlte.txt', encoding=None, dtype=None, names=True)
mpia = np.genfromtxt('cra2_abund_used_line_measurements_mpia_nlte.txt', encoding=None, dtype=None, names=True)
extra = np.genfromtxt('cra2_abund_used_line_measurements_extra_nlte.txt', encoding=None, dtype=None, names=True)

inspect = inspect[inspect['element']=='Na']
#mpia = mpia[mpia['element']!='Mg']
extra_k = extra[extra['element']=='K']
extra_mash = extra[extra['element']!='K']

outfile = open('line_by_line_measurements_table.txt', 'w')
outfile_tex = open('line_by_line_measurements_table_latex.txt', 'w')

for element, ions in abund_lines.items():
   for ion in ions:
       if (ion == 'convol') or (element == 'Tc'):
           continue
       if element not in ['C', 'N']:
           line_data = np.genfromtxt(f'line_data/{element}_{ion}_ts_lines.txt', encoding=None, dtype=None, usecols=(0, 1, 2), names=['wave', 'ep', 'loggf'])
       table_lines = table[(table['element']==element) & (np.in1d(table['wavelength'], ions[ion]))]
       for line in table_lines:
           if element not in ['C', 'N']:
               delta_wave = np.abs(line_data['wave']-line['wavelength'])
               ep = line_data['ep'][np.argmin(delta_wave)]
               if len(delta_wave[delta_wave < 0.1]) < 2:
                   loggf = f"{line_data['loggf'][np.argmin(delta_wave)]:.3f}"
               else:
                   loggf = ''
               for source, nltetable in zip([1,2,3,4],[inspect, mpia, extra_k, extra_mash]):
                   table_lines =  nltetable[(nltetable['star'] == line['star']) & (nltetable['element']==element) & (nltetable['wavelength'] == line['wavelength'])]
                   if len(table_lines) != 0:
                       nlte_corr = table_lines['nlte_correction'][0]
                       nlte_source = source
                       break
                   else:
                       nlte_corr = ''
                       nlte_source = ''
               print(f"{line['star']},{au.atomic_symbols[element]}.{au.ion_symbols[ion]},{line['wavelength']:.2f},{ep:.3f},{loggf},{line['logeps_chi2']:.3f},{nlte_corr},{nlte_source}", file=outfile)
               print(f"{line['star']} & {au.atomic_symbols[element]:3}.{au.ion_symbols[ion]} & {line['wavelength']:>8.2f} & {ep:5.3f} & {loggf} & {line['logeps_chi2']:5.3f} & {nlte_corr} & {nlte_source} \\\\", file=outfile_tex)

           else:
               if ion == 'C2':
                   species = '606.0'
               elif ion == 'CH':
                   species = '106.0'
               elif element == 'N':
                   species = '607.0'
               print(f"{line['star']},{species},{line['wavelength']:.2f},{''},{''},{line['logeps_chi2']:5.3f},{''},{''}", file=outfile)
               print(f"{line['star']} & {species:5} & {line['wavelength']:>8.2f} & {'     '} & {'      '} & {line['logeps_chi2']:5.3f} & & \\\\", file=outfile_tex)

outfile.close()
outfile_tex.close()