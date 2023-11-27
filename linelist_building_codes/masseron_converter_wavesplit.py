import numpy as np
import abund_utils as au
from astropy.io import fits
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


mols = ['12CH', '13CH']
molnames = ['12C1H', '13C1H']
elemnums = [106.001012, 106.001013]

for molname in molnames:
    if not os.path.isdir(f'{molname}_linelists'):
        os.mkdir(f'{molname}_linelists')

with fits.open('masseron2014_ch_linelist.fits') as file:
    fulldata = file[1].data
    
for wave_range in wave_ranges:

    data = fulldata[(fulldata['lam_Air'] >= wave_range[0]) & (fulldata['lam_Air'] < wave_range[1])]

    data = data[np.argsort(data['lam_Air'])]

    ev_in_cm1 = 8065.5410

    data['Elow'] /= ev_in_cm1
    data['Eup'] /= ev_in_cm1

    for mol, molname, elemnum in zip(mols, molnames, elemnums):
        isotab = data[data['mol'] == mol]
        nlines = len(isotab)
        filename = f'{molname}_linelists/{molname}_mas14_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}.turbo.air'
        with open(filename, 'w') as outfile:
            print(f"'{elemnum:11.6f} '    1 {nlines:10.0f}",
                  file=outfile)
            print("'Masseron 2014 CH line list '", file=outfile)

            for i in range(len(isotab)):
                wave = isotab['lam_Air'][i]

                # gu is the upper statistical weight used only for damping
                # and if raddmp is 0.  This the n degenerate states?

                ep = isotab['Elow'][i]
                gu = (2 * isotab['Ju'][i]) + 1
                labels = (f"{isotab['trans'][i]} {isotab['branch'][i]}"
                          f"{isotab['vu'][i]} {isotab['Ju'][i]} {isotab['symu'][i]} - "
                          f"{isotab['vl'][i]} {isotab['J'][i]} {isotab['syml'][i]}"
                          f"{isotab['Ref'][i]}")

                loggf = np.log10(isotab['gf'][i])

                # Van der Waals damping
                vdw_damp = 0.0
                raddamp = isotab['gamrad'][i]
                gamma_stark = 0.0

                lower_lev = "'x'"
                upper_lev = "'x'"
                description = (f"'{labels}'")
                print(f"{wave:10.3f} {ep:6.3f} {loggf:7.3f} "
                      f"{vdw_damp:9.2f} {gu:6.1f} {raddamp:9.2e} "
                      f"{gamma_stark:6.2f} {lower_lev} {upper_lev} 0.0 "
                      f"1.0 {description}",
                      file=outfile)

