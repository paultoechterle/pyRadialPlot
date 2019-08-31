
import csv

def read_radialplotter_file(filename):
    """ Parser P. Vermeesh RadialPlotter csv file

        returns: Ns, Ni, dpars as python lists"""

    with open(filename, "r") as f:
        file = csv.reader(f)
        name = next(file)
        line = next(file)
        zeta, zeta_err = float(line[0]), float(line[1])
        line = next(file)
        rhod, rhod_err = float(line[0]), int(float(line[1]))
        
        Ns = []
        Ni = []
        dpars = []
        
        for line in file:
            Ns.append(int(line[0]))
            Ni.append(int(line[1]))
            dpars.append(float(line[2]))

    return Ns, Ni, dpars
