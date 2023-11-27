import numpy as np
import abund_utils as au
from spec_tools import air_conversion
import os

if not os.path.isdir('kurucz_atoms'):
    os.mkdir('kurucz_atoms')

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
    'element',
    'blank',
    'ion',
    'e1',
    'J1',
    'label1',
    'e2',
    'J2',
    'label2',
    'gamma_rad',
    'gamma_stark',
    'gamma_vdw',
    'ref',
    'nlte_index_1',
    'nlte_index_2',
    'isotope_hpf',
    'hpf_frac',
    'isotope_iso',
    'iso_frac',
    'hpf_eshift_1',
    'hpf_eshift_2',
    'f',
    'hpf_f_1',
    'hpf_note_1',
    'minus',
    'hpf_f_2',
    'hpf_note_2',
    'line_strength_class',
    'autoionizing',
    'lande_g_even',
    'lande_g_odd',
    'iso_wave_shift'
]

widths = [
    11, 7, 3, 1, 2, 12, 5, 11, 12, 5, 11,
    6, 6, 6, 4, 2, 2, 3, 6, 3, 6, 5, 5, 1, 1, 1, 1,
    1, 1, 1, 3, 5, 5, 6
]

# kurucz = np.genfromtxt('gfallvac08oct17.dat.txt', encoding=None, dtype=None,
#                        names=names, delimiter=widths)

fulldata = np.genfromtxt('gfallvac08oct17_300-1100.dat.txt', encoding=None,
                         dtype=None, names=names, delimiter=widths)


fulldata = fulldata[fulldata['ion'] < 2.5]

ev_in_cm1 = 8065.5410

# Convert wavelengths from nm to angstroms and convert to air
fulldata['wavelength'] *= 10.
fulldata['wavelength'] = air_conversion(fulldata['wavelength'])

# Convert energies from cm^-1 to eV
fulldata['e1'] /= ev_in_cm1
fulldata['e2'] /= ev_in_cm1


# gamma_ad (Radiative damping) is in log (as are the gamma_vdw (van der Waals)
# and gamma_stark) but turbospec needs it in not log so convert it and reset
# any missing values (0 in Kurucz) to 0 (changing the default to 1e5)
no_gamma_rad = fulldata['gamma_rad'] == 0.
fulldata['gamma_rad'] = 10**fulldata['gamma_rad']
# fulldata['gamma_rad'][no_gamma_rad] = 0.0
fulldata['gamma_rad'][no_gamma_rad] = 1.0e5

# change the missing gamma_vdw to 2.5 which seems to follow what is typical of TS
no_gamma_vdw = fulldata['gamma_vdw'] == 0.
fulldata['gamma_vdw'][no_gamma_vdw] = 2.5

for wave_range in wave_ranges:
    data = fulldata[(fulldata['wavelength'] >= wave_range[0]) & 
                    (fulldata['wavelength'] < wave_range[1])]

    elements = np.unique(data['element'])
    # elem_order = np.argsort(
    #     np.array([au.atomic_sym_to_num(element) for element in elements]))
    # elements = elements[elem_order]

    with open(f'kurucz_atoms/kurucz_gfallvac08oct17_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}.atoms.turbo.air', 'w') as outfile:
        for element in elements:
            if element in [1, 2]:
                continue
            elemtab = data[data['element'] == element]
            ions = np.unique(elemtab['ion'])
            for ion in ions:
                iontab = elemtab[elemtab['ion'] == ion]
                isotopes = np.unique(iontab['isotope_iso'])
                for isotope in isotopes:
                    isotab = iontab[iontab['isotope_iso'] == isotope]
                    # If we chose to make cuts, cut them from isotab here so
                    # nlines is accurate
                    elemnum = element + isotope / 1000.
                    nlines = len(isotab)
                    print(f"'{elemnum:7.4f}             '{ion+1:5.0f}{nlines:10.0f}",
                          file=outfile)
                    element_sym = au.atomic_num_to_sym(element)
                    eleminfo = f'{element_sym} {au.ion_num_to_sym(ion)}'
                    print(f"'{eleminfo: <7}'", file=outfile)
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
                            labels = (f"{isotab['label2'][i].strip()} "
                                      f"{isotab['label1'][i].strip()}")

                        loggf = (isotab['loggf'][i]
                                 + isotab['hpf_frac'][i] + isotab['iso_frac'][i])
                        # Van der Waals damping
                        vdw_damp = isotab['gamma_vdw'][i]
                        raddamp = isotab['gamma_rad'][i]
                        gamma_stark = isotab['gamma_stark'][i]
                        if np.isnan(vdw_damp):
                            vdw_damp = 2.5
                        if np.isnan(raddamp):
                            raddamp = 1.0e5
                        if np.isnan(gamma_stark):
                            gamma_stark = 0.0
                        lower_lev = "'x'"
                        upper_lev = "'x'"
                        description = (f"'{element_sym} {au.ion_num_to_sym(ion)} "
                                       f"{labels}'")
                        print(f"{wave:10.3f} {ep:6.3f} {loggf:7.3f} "
                              f"{vdw_damp:9.2f} {gu:6.1f} {raddamp:9.2e} "
                              f"{gamma_stark:6.2f} {lower_lev} {upper_lev} 0.0 "
                              f"1.0 {description}",
                              file=outfile)

