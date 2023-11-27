import numpy as np

hydride_list = ['mooghyd03000', 'mooghyd04000', 'mooghyd06000', 'mooghyd10000']

hydride_list = ['/Users/hayesc/research/synthesis/linemake/mooglists/'
                + file for file in hydride_list]

wavelengths = []
isos = []
eps = []
loggfs = []
labels = []

for hydride_file in hydride_list:
    with open(hydride_file, 'r') as file:
        for line in file.readlines():
            if line[11:14] == '112':
                wavelengths += [float(line[0:10].strip())]
                isos += [line[18:20]]
                eps += [float(line[21:30].strip())]
                loggfs += [float(line[31:40].strip())]
                labels += [line[61:].strip()]


wavelengths = np.array(wavelengths)
isos = np.array(isos)
eps = np.array(eps)
loggfs = np.array(loggfs)
labels = np.array(labels)

for isotope in ['24', '25', '26']:
    wave_iso = wavelengths[isos == isotope]
    eps_iso = eps[isos == isotope]
    loggfs_iso = loggfs[isos == isotope]
    labels_iso = labels[isos == isotope]

    with open(f'{isotope}MgH_hinkle.turbo.air', 'w') as outfile:
        elemnum = f'112.0010{isotope}'
        nlines = len(wave_iso)
        print(f"'{elemnum} '    1 {nlines:10.0f}",
              file=outfile)
        print("'Placco Linemake Hinkle 2013 '", file=outfile)
        for wave, ep, loggf, label in zip(wave_iso, eps_iso, loggfs_iso, labels_iso):
            vdw_damp = 0.0
            raddamp = 1.0
            gamma_stark = 0.0
            gu = 0.0

            lower_lev = "'x'"
            upper_lev = "'x'"
            description = (f"'{label}'")
            print(f"{wave:10.3f} {ep:6.3f} {loggf:7.3f} "
                  f"{vdw_damp:9.2f} {gu:6.1f} {raddamp:9.2e} "
                  f"{gamma_stark:6.2f} {lower_lev} {upper_lev} 0.0 "
                  f"1.0 {description}",
                  file=outfile)
