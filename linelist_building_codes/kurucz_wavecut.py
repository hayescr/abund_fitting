waverange = [300., 1100.]

outfile_name = f'gfallvac08oct17_{waverange[0]:.0f}-{waverange[1]:.0f}.dat.txt'

with open('gfallvac08oct17.dat.txt', 'r') as infile:
    with open(outfile_name, 'w') as outfile:
        for line in infile.readlines():
            wave = float(line[0:11].strip())
            if (wave > waverange[0]) & (wave < waverange[1]):
                print(line, end='', file=outfile)
