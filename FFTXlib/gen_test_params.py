import re, sys

info = """
This script collects the relevant simulation parameters and generates a 
set of input options for `fft_test.x` that replicate the execution of
vloc_psi done in the production run.
This simplifies the identification of the optimal fft tasking paramete.

Usage: python run_test.py pw_out_file
"""

match_alat = re.compile(r'lattice parameter \(alat\)\s+=\s+([+-]?([0-9]*[.])?[0-9]+)')
match_nbnd = re.compile(r'number of Kohn-Sham states=\s+([+-]?([0-9]*[.])?[0-9]+)')
match_ecutwfc = re.compile(r'kinetic-energy cutoff\s+=\s+([+-]?([0-9]*[.])?[0-9]+)')
match_ecutrho = re.compile(r'charge density cutoff\s+=\s+([+-]?([0-9]*[.])?[0-9]+)')
match_k = re.compile(r'number of k points=\s+(\d+)')


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print(info)
    else:
        with open(sys.argv[1],'r') as f:
            data = f.read(30000)
            gamma = False
            maxk = ''
            alat = match_alat.findall(data)
            nbnd = match_nbnd.findall(data)
            ewfc = match_ecutwfc.findall(data)
            erho = match_ecutrho.findall(data)
            if len(alat[0]) == 0 or len(nbnd[0]) == 0 or len(ewfc[0]) == 0 or len(erho[0]) == 0:
                print("Could not parse file. Sorry.")
            alat = alat[0][0]; nbnd=nbnd[0][0]; ewfc=ewfc[0][0];erho=erho[0][0]
            a1 = []
            a2 = []
            a3 = []
            lines = data.splitlines()
            for i, line in enumerate(lines):
                if 'gamma-point specific algorithms are used' in line:
                    gamma = True
                if 'crystal axes' in line:
                    a1 = [float(x)*float(alat) for x in lines[i+1].split()[3:6]]
                    a2 = [float(x)*float(alat) for x in lines[i+2].split()[3:6]]
                    a3 = [float(x)*float(alat) for x in lines[i+3].split()[3:6]]
                if 'number of k points' in line:
                    nk = int(match_k.findall(line)[0])
                    if '2pi/alat' in lines[i+1]:
                        nrm2 = 0
                        for k in range(nk):
                            kn, v = re.findall(r"\(([-\d\s\.]*)\)", lines[i+k+2])
                            v2 = [float(x)**2 for x in v.split()]
                            if sum(v2) > nrm2:
                                maxk=v
                    
            print ("To analize performances run with:")
            buf = ("mpirun -np X ./fft_test.x -ntg Y -ecutwfc {ewfc} -ecutrho {erho} " + \
                  "-av1 {av1} -av2 {av2} -av3 {av3} -nbnd {nbnd} -gamma {gamma}").format(\
                       ewfc=ewfc, erho=erho, \
                       av1=' '.join([str(x) for x in a1]), \
                       av2=' '.join([str(x) for x in a2]), \
                       av3=' '.join([str(x) for x in a3]), \
                       nbnd=nbnd, gamma=('.true.' if gamma else '.false.'))
            if maxk:
                buf += " -kmax "+maxk
            print(buf)

