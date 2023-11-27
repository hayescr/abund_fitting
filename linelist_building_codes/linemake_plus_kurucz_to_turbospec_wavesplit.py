import numpy as np
import abund_utils as au
from spec_tools import air_conversion
import os


def turbo_fdamp(element, ion):
    # From Turbospectrum 2019 vald3line-BPz-freeformat.f
    
    if ion == 1:
        if element == 11:
            # Holweger 1971 A&A 10, 128
            fdamp=2.0
        elif element == 14:
            # Holweger 1973 A&A 26, 275
            fdamp=1.3
        elif element == 20:
            # O'Neill & Smith 1980 A&A 81, 100
            fdamp=1.8
        elif element == 26:
            # Simmons & Blackwell 1982 A&A 112, 209
            # Magain & Zhao 1996 A&A 305, 245
            fdamp=1.4
        else:
            # Mackle et al. 1975 A&A 38, 239
            fdamp=2.5
    else:
        if element == 20:
            # from fit of H&K in the HM model to the fts intensity spectrum
            fdamp=1.4
        elif element == 38:
            # from fit of Sr II 4077.724 in the HM model to the fts intensity spectrum
            fdamp=1.8
        elif element == 56:
            # Holweger & Muller 1974 Solar Physics 39, 19
            fdamp=3.0
        else:
            # Mackle et al. 1975 A&A 38, 239
            fdamp=2.5
    return fdamp


if not os.path.isdir('linemake_atoms'):
    os.mkdir('linemake_atoms')


wave_ranges2 = [
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

names2 = [
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

widths2 = [
    11, 7, 3, 1, 2, 12, 5, 11, 12, 5, 11,
    6, 6, 6, 4, 2, 2, 3, 6, 3, 6, 5, 5, 1, 1, 1, 1,
    1, 1, 1, 3, 5, 5, 6
]

# kurucz = np.genfromtxt('gfallvac08oct17.dat.txt', encoding=None, dtype=None,
#                        names=names, delimiter=widths)

fulldata2 = np.genfromtxt('gfallvac08oct17_300-1100.dat.txt', encoding=None,
                         dtype=None, names=names2, delimiter=widths2)


fulldata2 = fulldata2[fulldata2['ion'] < 2.5]

ev_in_cm1 = 8065.5410

# Convert wavelengths from nm to angstroms and convert to air
fulldata2['wavelength'] *= 10.
fulldata2['wavelength'] = air_conversion(fulldata2['wavelength'])

# Convert energies from cm^-1 to eV
fulldata2['e1'] /= ev_in_cm1
fulldata2['e2'] /= ev_in_cm1


# gamma_ad (Radiative damping) is in log (as are the gamma_vdw (van der Waals)
# and gamma_stark) but turbospec needs it in not log so convert it and reset
# any missing values (0 in Kurucz) to 0 (changing the default to 1e5)
no_gamma_rad = fulldata2['gamma_rad'] == 0.
fulldata2['gamma_rad'] = 10**fulldata2['gamma_rad']
# fulldata['gamma_rad'][no_gamma_rad] = 0.0
fulldata2['gamma_rad'][no_gamma_rad] = 1.0e5

# change the missing gamma_vdw to 2.5 which seems to follow what is typical of TS
#no_gamma_vdw = fulldata['gamma_vdw'] == 0.
#fulldata['gamma_vdw'][no_gamma_vdw] = 2.5

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


names = ['wavelength', 'mass', 'ep', 'loggf', 'vdw', 'dissociation', 'ew', 'source']
widths = [10, 10, 10, 10, 10, 10, 10, 10]

path = '/Users/hayesc/research/synthesis/linemake/linemake_ll_moog/'

for wave_range in wave_ranges:
    data2 = fulldata2[(fulldata2['wavelength'] >= wave_range[0]) & 
                    (fulldata2['wavelength'] < wave_range[1])]


    # filename = 'euII.mooghfs'
    filename = f'linemake_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}.atoms.moog'
    new_filename = f'linemake_atoms/linemake_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}.atoms.turbo.air'

    data = np.genfromtxt(path+filename, dtype=None, names=names, encoding=None, 
                         skip_header=1, delimiter=widths)

    data = data[data['wavelength'] > 0.]

    elem_info = [f'{eli:0.4f}' for eli in data['mass']]

    elem_arr = np.array([int(f'{eli.split(".")[0]}') for eli in elem_info])
    mass_arr = np.array([int(f'{eli.split(".")[1][1:]}') for eli in elem_info])
    ion_arr = np.array([int(f'{eli.split(".")[1][0]}') for eli in elem_info])
    elements = np.unique(elem_arr)

    with open(new_filename, 'w') as outfile:
        for element in elements:
            if element == 1:
                continue
            elem_tab = data[elem_arr==element]
            mass_elem_tab = mass_arr[elem_arr==element]
            ion_elem_tab = ion_arr[elem_arr==element]
            ions = np.unique(ion_elem_tab)

            elemtab2 = data2[data2['element'] == element]

            for ion in ions:
                iontab2 = elemtab2[elemtab2['ion'] == ion]

                ion_tab = elem_tab[ion_elem_tab == ion]
                mass_ion_tab = mass_elem_tab[ion_elem_tab == ion]
                isotopes = np.unique(mass_ion_tab)
                fdamp = turbo_fdamp(element, ion+1)
                for isotope in isotopes:
                    isotab2 = iontab2[iontab2['isotope_iso'] == isotope]

                    kurucz_ep = np.minimum(isotab2['e1'], isotab2['e2'])
                    kurucz_j = isotab2['J2']
                    states_flipped = np.abs(isotab2['e1']) > np.abs(isotab2['e2'])
                    kurucz_j[states_flipped] = isotab2['J1'][states_flipped]
                    kurucz_gu = (2 * kurucz_j) + 1

                    isotab = ion_tab[mass_ion_tab == isotope]
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
                        #if np.abs(isotab['e1'][i]) < np.abs(isotab['e2'][i]):
                        #    ep = np.abs(isotab['e1'][i])
                        #    gu = (2 * isotab['J2'][i]) + 1
                        #    labels = (f"{isotab['label1'][i].strip()} "
                        #              f"{isotab['label2'][i].strip()}")

                        loggf = isotab['loggf'][i]
                        ep = isotab['ep'][i]
                        # Van der Waals damping

                        kurucz_matches = (np.abs(isotab2['wavelength']-wave) < 0.03)
                        kurucz_matches &= (np.abs(kurucz_ep - ep) < 0.03)

                        ktab = isotab2[kurucz_matches]
                        ktab_ep = kurucz_ep[kurucz_matches]
                        ktab_gu = kurucz_gu[kurucz_matches]

                        diff = np.abs(ktab['wavelength']-wave) + np.abs(ktab_ep - ep)

                        if len(ktab) > 1:
                            best_match = np.argmin(diff)
                            gu = ktab_gu[best_match]
                            vdw_damp = ktab['gamma_vdw'][best_match]
                            raddamp = ktab['gamma_rad'][best_match]
                            gamma_stark = ktab['gamma_stark'][best_match]
                        elif len(ktab) == 1:
                            gu = ktab_gu[0]
                            vdw_damp = ktab['gamma_vdw'][0]
                            raddamp = ktab['gamma_rad'][0]
                            gamma_stark = ktab['gamma_stark'][0]
                        else:
                            gu = 0.0
                            vdw_damp = fdamp
                            raddamp = 1.0e5
                            gamma_stark = 0.0

                        labels = isotab['source'][i]
                        if np.isnan(gu):
                            gu = 0.0
                        if np.isnan(vdw_damp) | (vdw_damp == 0.0):
                            vdw_damp = fdamp
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

