import os

class LinelistManager:
    available_molecules = {'16O1H' : 'exomol', 
                           '12C14N' : 'exomol', 
                           '12C15N' : 'exomol', 
                           '13C14N' : 'exomol', 
                           '13C15N' : 'exomol', 
                           '14N1H' : 'exomol', 
                           '12C12C' : 'exomol', 
                           '12C13C' : 'exomol', 
                           '13C13C' : 'exomol', 
                           '28Si1H' : 'exomol', 
                           '29Si1H' : 'exomol', 
                           '30Si1H' : 'exomol', 
                           '56Fe1H' : 'exomol', 
                           '40Ca1H' : 'exomol', 
                           '24Mg1H' : 'kurucz_all', 
                           '25Mg1H' : 'kurucz_all', 
                           '26Mg1H' : 'kurucz_all',
                           '12C1H' : 'mas14',
                           '13C1H' : 'mas14',
                          }
    
    wave_ranges = [[3000, 3500],
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
    
    atoms_suffix = '.atoms.turbo.air'
    mols_suffix = '.turbo.air'
    h_line_data = 'DATA/Hlinedata'
    
    'COM-v19.1/'

    def __init__(self, hlines=True, atoms=True, molecules=None, source='kurucz',
                 turbospec_path='.', ts_linelist_path='COM-v19.1/new_linelists'):
        self.atoms = atoms
        self.hlines = hlines
        molecule_list = self.check_molecules(molecules)
        self.molecules = molecule_list
        self.ts_linelist_path = ts_linelist_path
        self.turbospec_path = turbospec_path
        if source == 'kurucz':
            self.atoms_prefix = 'kurucz_gfallvac08oct17_'
            self.atoms_folder = 'kurucz_atoms'
        elif source == 'linemake':
            self.atoms_prefix = 'linemake_'
            self.atoms_folder = 'linemake_atoms'

        
        
    def check_molecules(self, molecules):
        if molecules is None:
            return []
        else:
            for molecule in molecules:
                if molecule not in self.available_molecules:
                    print('Sorry there is no linelist for {molecule}')
            molecule_list = [molecule for molecule in molecules if molecule in self.available_molecules]
            return molecule_list
        
    def build_library(self, synth_range=None):
        
        ll_path = f'{self.turbospec_path}/{self.ts_linelist_path}'
        linelists = []

        if self.hlines:
            linelists += [f'{self.turbospec_path}/{self.h_line_data}']

        
        if synth_range is None:
            if self.atoms:
                linelists += [f'{ll_path}/{self.atoms_prefix}300-1100{self.atoms_suffix}']
            for molecule in self.molecules:
                linelists += [f'{ll_path}/{molecule}_{self.available_molecules[molecule]}{self.mols_suffix}']
        else:
            library_ranges = [f'{wave_range[0]/10:.0f}-{wave_range[1]/10:.0f}' for wave_range in self.wave_ranges if self.overlap(wave_range, synth_range)]
            if self.atoms:
                linelists += [f'{ll_path}/{self.atoms_folder}/{self.atoms_prefix}{library_range}{self.atoms_suffix}' for library_range in library_ranges]
            for molecule in self.molecules:
                linelists += [f'{ll_path}/{molecule}_linelists/{molecule}_{self.available_molecules[molecule]}_{library_range}{self.mols_suffix}' for library_range in library_ranges]

        linelists = [linelist for linelist in linelists if os.path.exists(linelist)]
                
        return linelists
                
    def overlap(self, wave_range, synth_range):
        return (min(wave_range[1], synth_range[1]) - max(wave_range[0], synth_range[0])) > 0