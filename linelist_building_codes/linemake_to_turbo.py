import numpy as np
import abund_utils as au
from spec_tools import air_conversion

path = '/Users/hayesc/research/synthesis/linemake/mooglists/'

# filename = 'euII.mooghfs'
filename = 'dyII.moog'

ion = 'II'

data = np.genfromtxt(path+filename, dtype=None, names=['wavelength', 'mass', 'ep', 'loggf', 'source'], encoding=None)

data = data[data['wavelength'] > 0.]

elem_info = [f'{eli:0.4f}' for eli in data['mass']]

elem_arr = np.array([int(f'{eli.split(".")[0]}') for eli in elem_info])
mass_arr = np.array([int(f'{eli.split(".")[1][1:]}') for eli in elem_info])
ion_arr = np.array([int(f'{eli.split(".")[1][0]}') for eli in elem_info])
elements = np.unique(elem_arr)

with open(f'linemake.{filename}.atoms.turbo.air', 'w') as outfile:
    for element in elements:
        elem_tab = data[elem_arr==element]
        mass_elem_tab = mass_arr[elem_arr==element]
        ion_elem_tab = ion_arr[elem_arr==element]
        ions = np.unique(ion_elem_tab)
        for ion in ions:
            ion_tab = elem_tab[ion_elem_tab == ion]
            mass_ion_tab = mass_elem_tab[ion_elem_tab == ion]
            isotopes = np.unique(mass_ion_tab)
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
                    vdw_damp = 2.5
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
