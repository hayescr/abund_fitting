element = 'Fe'
ion = 'I'

line_lists=['linemake_300-350.atoms.turbo.air',

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