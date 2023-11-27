import numpy as np
from spec_tools import air_conversion


def calc_loggf(wavenumber, jlow, Ajj):
    # enter the wavenumber in 1/cm and Ajj in 1/s
    mass_e = 9.1093837015e-31  # kg
    charge_e = 1.602176634e-19  # C
    vac_electric_perm = 8.8541878128e-12  # F/m
    speed_of_light = 299792458.  # m/s
    pi_num = 3.14159265
    constant_factor = (mass_e * vac_electric_perm
                       * speed_of_light) / (2 * pi_num * charge_e**2.)
    wavenumber_mks = 100 * wavenumber
    gf = constant_factor * ((2 * jlow) + 1) * Ajj / (wavenumber_mks**2.)
    # print(gf)
    return np.log10(gf)


def exomol_reader(states_filename, trans_filename, wave_range=None):

    state_names = ['id', 'e', 'gtot', 'j']
    states = np.genfromtxt(states_filename, names=state_names,
                           dtype=None, encoding=None, usecols=(0, 1, 2, 3))

    trans_names = ['id_u', 'id_l', 'Ajj', 'wave_num']
    trans = np.genfromtxt(trans_filename, names=trans_names,
                          dtype=None, encoding=None, usecols=(0, 1, 2, 3))

    ajj = []
    wave_num = []
    el = []
    eu = []
    jl = []
    ju = []
    labels = []

    wavelengths = 1e8 / trans['wave_num']
    wavelengths = air_conversion(wavelengths)

    if wave_range is not None:
        cuts = (wavelengths < wave_range[1]) & (wavelengths >= wave_range[0])
        wavelengths = wavelengths[cuts]
        trans = trans[cuts]

    for i in range(len(trans)):
        ajj += [trans['Ajj'][i]]
        wave_num += [trans['wave_num'][i]]
        state_l = states[states['id'] == trans['id_l'][i]]
        state_u = states[states['id'] == trans['id_u'][i]]
        el += [state_l['e'][0]]
        eu += [state_u['e'][0]]
        jl += [state_l['j'][0]]
        ju += [state_u['j'][0]]

    ajj = np.array(ajj)
    wave_num = np.array(wave_num)
    el = np.array(el)
    eu = np.array(eu)
    jl = np.array(jl)
    ju = np.array(ju)
    labels = np.array(labels)

    ev_in_cm1 = 8065.5410

    eps = el / ev_in_cm1
    loggfs = calc_loggf(wave_num, jl, ajj)
    gus = (2 * ju) + 1

    wave_sort = np.argsort(wavelengths)

    wavelengths = wavelengths[wave_sort]
    eps = eps[wave_sort]
    loggfs = loggfs[wave_sort]
    gus = gus[wave_sort]

    return wavelengths, eps, loggfs, gus
