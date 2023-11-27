from exomol_converter import exomol_reader
import os

mol_filename_list_base = ['16O-1H__MoLLIST', '12C-14N__Trihybrid',
                     '12C-15N__Trihybrid', '13C-14N__Trihybrid',
                     '13C-15N__Trihybrid', '14N-1H__MoLLIST', '12C2__8states',
                     '12C-13C__8states', '13C2__8states', '28Si-1H__SiGHTLY',
                     '29Si-1H__SiGHTLY', '30Si-1H__SiGHTLY',
                     '56Fe-1H__MoLLIST', '40Ca-1H__MoLLIST']
elemnums = ['0108.001016', '0607.012014', '0607.012015', '0607.013014',
            '0607.013015', '0107.001014', '0606.012012', '0606.012013',
            '0606.013013', '0114.001028', '0114.001029', '0114.001030',
            '0126.001056', '0120.001040']
descriptions = ['Exomol Brooke 2016 Yousefi 2018',
                'Exomol Brooke 2014 Syme 2021',
                'Exomol Brooke 2014 Syme 2021',
                'Exomol Brooke 2014 Syme 2021',
                'Exomol Brooke 2014 Syme 2021',
                'Exomol Brooke 2014 2018 Fernando 2018',
                'Exomol Yurchenko 2018 McKemmish 2020',
                'Exomol Yurchenko 2018 McKemmish 2020',
                'Exomol Yurchenko 2018 McKemmish 2020',
                'Exomol Yurchenko 2018',
                'Exomol Yurchenko 2018',
                'Exomol Yurchenko 2018',
                'Exomol Dulick 2003',
                'Exomol Li 2012 Shayesteh 2004 Alavi 2017',
                ]
outfile_titles_base = ['16O1H', '12C14N', '12C15N', '13C14N', '13C15N', '14N1H',
                  '12C12C', '12C13C', '13C13C', '28Si1H', '29Si1H', '30Si1H',
                  '56Fe1H', '40Ca1H']

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

for name in outfile_titles_base:
    if not os.path.isdir(f'{name}_linelists'):
        os.mkdir(f'{name}_linelists')



for wave_range in wave_ranges:

    outfile_titles = [f'{name}_linelists/{name}_exomol_{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}' for name in outfile_titles_base]
    mol_filename_list = [f'exomol_files/{name}' for name in mol_filename_list_base]

    '''
    Things I need:
    elemnum
    outfile_title
    source_description
    '''

    for name, elemnum, source_description, outfile_title in zip(mol_filename_list,
                                                                elemnums,
                                                                descriptions,
                                                                outfile_titles):
        states_name = f'{name}.states'
        trans_name = f'{name}.trans'
        wavelengths, eps, loggfs, gus = exomol_reader(states_name, trans_name,
                                                      wave_range=wave_range)

        with open(f'{outfile_title}.turbo.air', 'w') as outfile:
            nlines = len(wavelengths)
            print(f"'{elemnum} '    1 {nlines:10.0f}",
                  file=outfile)
            print(f"'{source_description} '", file=outfile)
            for wave, ep, loggf, gu in zip(wavelengths, eps, loggfs, gus):
                # Van der Waals damping
                vdw_damp = 0
                raddamp = 0.0
                gamma_stark = 0
                lower_lev = "'x'"
                upper_lev = "'x'"
                description = ("")
                print(f"{wave:10.3f} {ep:6.3f} {loggf:7.3f} "
                      f"{vdw_damp:9.2f} {gu:6.1f} {raddamp:9.2e} "
                      f"{gamma_stark:6.2f} {lower_lev} {upper_lev} 0.0 "
                      f"1.0 {description}",
                      file=outfile)

        print(f'{outfile_title} done')

