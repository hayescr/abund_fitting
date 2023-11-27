element = 'Fe'
ion = 'I'

line_lists=['linemake_300-350.atoms.turbo.air',    'linemake_350-400.atoms.turbo.air',    'linemake_400-450.atoms.turbo.air',    'linemake_450-500.atoms.turbo.air',    'linemake_500-550.atoms.turbo.air',    'linemake_550-600.atoms.turbo.air',    'linemake_600-650.atoms.turbo.air',    'linemake_650-700.atoms.turbo.air',    'linemake_700-750.atoms.turbo.air',    'linemake_750-800.atoms.turbo.air',    'linemake_800-850.atoms.turbo.air',    'linemake_850-900.atoms.turbo.air',    'linemake_900-950.atoms.turbo.air',    'linemake_950-1000.atoms.turbo.air',    'linemake_1000-1050.atoms.turbo.air',    'linemake_1050-1100.atoms.turbo.air']

outfile = open(f'{element}_{ion}_ts_lines.txt', 'w')

for linelist in line_lists:
    filename = f'/Users/hayesc/research/Turbospectrum2019/COM-v19.1/new_linelists/linemake_atoms/{linelist}'

    with open(filename, 'r') as file_in:
        printout = False
        for line in file_in.readlines():
            if printout:
                if line.startswith('  '):
                    print(line, end='', file=outfile)
                else:
                    if element_found:
                        break
            
            if line.startswith(f"'{element} {ion}"):
                printout = True
                element_found = True

outfile.close()