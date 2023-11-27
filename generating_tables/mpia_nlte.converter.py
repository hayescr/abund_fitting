import io

longstring = '''lines(A)   4890.765   4891.502   4903.318   4994.130   5049.828   5079.228   5166.290   5191.462   5194.949   5216.280   5225.533   5232.950   5250.216   5281.798   5339.937   5383.369   5397.135   5405.786   5415.199   5434.524   5569.618   5572.853   5586.766   5615.656   5624.551   6065.492   6136.627   6137.702   6191.569   6230.735   6252.555   6265.134   6335.331   6393.601   6400.011   6430.846   6677.995   4515.338   4583.838   4923.932   5197.575
star1      0.182      0.134      0.120      0.235      0.155      0.154      0.223      0.123      0.154      0.211      0.000      0.127      0.056      0.108      0.076      0.102      0.195      0.170      0.098      0.231      0.073      0.073      0.077      0.074      0.068      0.097      0.061      0.049      0.132      0.114      0.120      0.094      0.115      0.134      0.113      0.000      0.094     -0.011     -0.012      0.030     -0.002
star2      0.262      0.227      0.210      0.336      0.252      0.238      0.336      0.216      0.244      0.341      0.000      0.225      0.093      0.198      0.132      0.162      0.348      0.275      0.161      0.401      0.127      0.126      0.134      0.132      0.123      0.171      0.108      0.085      0.238      0.201      0.216      0.169      0.208      0.246      0.186      0.000      0.167     -0.015     -0.014      0.018     -0.005
'''

outfilename='nlte.tmp.txt'

filelines = io.StringIO(longstring).readlines()

#with open(filename, 'r') as file:
#    filelines = file.readlines()

with open(outfilename, 'w') as outfile:
    lines = filelines[0].split()[1:]
    print('# wavelength', file=outfile, end=' ')
    for fileline in filelines[1:]:
        print(fileline.split()[0], end=' ', file=outfile)
    print('', file=outfile)
    for i, line in enumerate(lines):
        print(line, file=outfile, end=' ')
        for fileline in filelines[1:]:
            print(f'{float(fileline.split()[i+1]):6.3f}', file=outfile, end=' ')
        print('', file=outfile)


        