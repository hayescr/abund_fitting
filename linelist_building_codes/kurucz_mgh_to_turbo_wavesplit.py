import numpy as np
import abund_utils as au
from spec_tools import air_conversion
import os


wave_ranges = [
    [3000, 3500],
    [3500, 4000],
    [4000, 4500],
    [4500, 5000],
    [5000, 5500],
    [5500, 6000],
    [6000, 6500],
    [6500, 7000],
    [7000, 7500],
    [7500, 8000],
    [8000, 8500],
    [8500, 9000],
    [9000, 9500],
    [9500, 10000],
    [10000, 10500],
    [10500, 11000],
]


names = [
    'wavelength',
    'loggf',
    'J1',
    'e1',
    'J2',
    'e2',
    'label1',
    'label2',
    'isotope',
    'blank',
    'wave_num'
]


widths = [
    10, 7, 5, 10, 5, 11, 9, 8, 5, 4, 12
]  # For the mgh.asc.txt

# widths = [
#     9, 7, 5, 10, 5, 11, 9, 8, 5, 4, 12
# ]  # For the mgh.asc.txt

# with open('mgh.asc.txt') as file:
#     for line in file.readlines():
#         line_s = line.split()
#         if len(line_s) < 11:
#             print(line)
# kurucz = np.genfromtxt('gfallvac08oct17.dat.txt', encoding=None, dtype=None,
#                        names=names, delimiter=widths)

fulldata = np.genfromtxt('mgh.asc.txt', encoding=None,
                     dtype=None, names=names, delimiter=widths)

print(fulldata['label1'])
print(fulldata['wave_num'])


ev_in_cm1 = 8065.5410

# Convert wavelengths from nm to angstroms and convert to air
fulldata['wavelength'] *= 10.

# Convert energies from cm^-1 to eV
fulldata['e1'] /= ev_in_cm1
fulldata['e2'] /= ev_in_cm1

for isotope in [24, 25, 26]:
    if not os.path.isdir(f'{isotope}Mg1H_linelists'):
        os.mkdir(f'{isotope}Mg1H_linelists')


for wave_range in wave_ranges:
    data = fulldata[(fulldata['wavelength'] >= wave_range[0]) & 
                    (fulldata['wavelength'] < wave_range[1])]


    # elem_order = np.argsort(
    #     np.array([au.atomic_sym_to_num(element) for element in elements]))
    # elements = elements[elem_order]

    for isotope in [24, 25, 26]:
        isotab = data[data['isotope'] == isotope]

        with open(f'{isotope}Mg1H_linelists/{isotope}Mg1H_kurucz_all_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}.turbo.air', 'w') as outfile:
            elemnum = f'112.0010{isotope}'
            nlines = len(isotab)
            print(f"'{elemnum} '    1 {nlines:10.0f}",
                  file=outfile)
            print("'Kurucz Web '", file=outfile)
            for i in range(len(isotab)):
                wave = isotab['wavelength'][i]
                # Theoretically calculated energies are negative so,
                # we need to figure out which energy level is truly lower,
                # since the minimum one is in e1 and the max is in e2,
                # but we want the min and max of the abs for low and high

                # gu is the upper statistical weight used only for damping
                # and if raddmp is 0.  This the n degenerate states?
                if np.abs(isotab['e1'][i]) < np.abs(isotab['e2'][i]):
                    ep = np.abs(isotab['e1'][i])
                    gu = (2 * isotab['J2'][i]) + 1
                    labels = (f"{isotab['label1'][i].strip()} "
                              f"{isotab['label2'][i].strip()}")
                else:
                    ep = np.abs(isotab['e2'][i])
                    gu = (2 * isotab['J1'][i]) + 1
                    labels = (f"{isotab['label1'][i].strip()} "
                              f"{isotab['label2'][i].strip()}")

                loggf = isotab['loggf'][i]
                # Van der Waals damping
                vdw_damp = 0
                raddamp = 0.0
                gamma_stark = 0
                lower_lev = "'x'"
                upper_lev = "'x'"
                description = (f"'{labels}'")
                print(f"{wave:10.3f} {ep:6.3f} {loggf:7.3f} "
                      f"{vdw_damp:9.2f} {gu:6.1f} {raddamp:9.2e} "
                      f"{gamma_stark:6.2f} {lower_lev} {upper_lev} 0.0 "
                      f"1.0 {description}",
                      file=outfile)

