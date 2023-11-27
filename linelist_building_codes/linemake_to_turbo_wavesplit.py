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
            for ion in ions:
                ion_tab = elem_tab[ion_elem_tab == ion]
                mass_ion_tab = mass_elem_tab[ion_elem_tab == ion]
                isotopes = np.unique(mass_ion_tab)
                fdamp = turbo_fdamp(element, ion+1)
                for isotope in isotopes:
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

                        gu = 0.0
                        labels = isotab['source'][i]

                        loggf = isotab['loggf'][i]
                        ep = isotab['ep'][i]
                        # Van der Waals damping
                        vdw_damp = fdamp
                        raddamp = 1.0e5
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

