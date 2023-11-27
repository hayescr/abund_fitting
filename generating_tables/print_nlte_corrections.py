import numpy as np
from abund_lines import *

def nlte_corrections():

    mpia_nlte = np.genfromtxt('cra2_abund_used_line_measurements_mpia_nlte.txt', dtype=None, encoding=None, names=True)
    inspect_nlte = np.genfromtxt('cra2_abund_used_line_measurements_inspect_nlte.txt', dtype=None, encoding=None, names=True)
    extra_nlte = np.genfromtxt('cra2_abund_used_line_measurements_extra_nlte.txt', dtype=None, encoding=None, names=True)

    abund_data = np.genfromtxt('cra2_abund_error_table_ions.txt', dtype=None, encoding=None, names=True)

    mpia_elements = np.unique(mpia_nlte['element'])
    inspect_elements = np.unique(inspect_nlte['element'])
    extra_elements = np.unique(extra_nlte['element'])

    inspect_use_elements = ['Na']
    extra_use_elements = ['K', 'Ba']

    mpia_nlte_corrections = {}
    inspect_nlte_corrections = {}
    extra_nlte_corrections = {}

    star_list = ['cra2_gdr32672', 'cra2_gdr36912']

    print('MPIA NLTE Corrections')
    for star in star_list:
        print(star)
        star_nltes = mpia_nlte[mpia_nlte['star']==star]
        mpia_nlte_corrections[star]={}
        for element in mpia_elements:
            element_nltes = star_nltes[star_nltes['element']==element]
            element_nlte_dict = {}
            for ion in abund_lines[element]:
                if ion == 'convol':
                    continue
                ion_nltes = element_nltes[np.in1d(element_nltes['wavelength'], abund_lines[element][ion])]
                if len(ion_nltes) > 0:
                    element_nlte_dict[ion] = np.mean(ion_nltes['nlte_correction'])
                    print(element, ion, len(ion_nltes['nlte_correction']), f'{element_nlte_dict[ion]:.2f}')
            mpia_nlte_corrections[star][element] = element_nlte_dict
        print()

    print()
    print('INPSECT NLTE Corrections')
    for star in star_list:
        print(star)
        star_nltes = inspect_nlte[inspect_nlte['star']==star]
        inspect_nlte_corrections[star]={}
        for element in inspect_elements:
            element_nltes = star_nltes[star_nltes['element']==element]
            element_nlte_dict = {}
            for ion in abund_lines[element]:
                if ion == 'convol':
                    continue
                ion_nltes = element_nltes[np.in1d(element_nltes['wavelength'], abund_lines[element][ion])]
                if len(ion_nltes) > 0:
                    element_nlte_dict[ion] = np.mean(ion_nltes['nlte_correction'])
                    print(element, ion, len(ion_nltes['nlte_correction']), f'{element_nlte_dict[ion]:.2f}')
            inspect_nlte_corrections[star][element] = element_nlte_dict
        print()

    print()
    print('Extra NLTE Corrections')
    for star in star_list:
        print(star)
        star_nltes = extra_nlte[extra_nlte['star']==star]
        extra_nlte_corrections[star]={}
        for element in extra_elements:
            element_nltes = star_nltes[star_nltes['element']==element]
            element_nlte_dict = {}
            for ion in abund_lines[element]:
                if ion == 'convol':
                    continue
                ion_nltes = element_nltes[np.in1d(element_nltes['wavelength'], abund_lines[element][ion])]
                if len(ion_nltes) > 0:
                    element_nlte_dict[ion] = np.mean(ion_nltes['nlte_correction'])
                    print(element, ion, len(ion_nltes['nlte_correction']), f'{element_nlte_dict[ion]:.2f}')
            extra_nlte_corrections[star][element] = element_nlte_dict
        print()


    final_nlte_corrections = dict(mpia_nlte_corrections)

    for star in star_list:
        for element in inspect_use_elements:
            if element in inspect_nlte_corrections[star]:
                final_nlte_corrections[star][element] = inspect_nlte_corrections[star][element]
    for star in star_list:
        for element in extra_use_elements:
            if element in extra_nlte_corrections[star]:
                final_nlte_corrections[star][element] = extra_nlte_corrections[star][element]

    return final_nlte_corrections
    
if __name__ == "__main__":
    nlte_corrections()