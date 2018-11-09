import numpy as np
import sys
from Bio.PDB import PDBParser
from scipy.stats.mstats import mquantiles

jsk = sys.argv[1]


waters = ['OW', 'HW1', 'HW2']
heads = ['N', 'C13', 'H13A', 'H13B', 'H13C', 'C14', 'H14A', 'H14B', 'H14C',
       'C15', 'H15A', 'H15B', 'H15C', 'C12', 'H12A', 'H12B', 'C11', 'H11A',
       'H11B', 'P', 'O13', 'O14', 'O11', 'O12', 'C1', 'HA', 'HB', 'C2', 'HS',
       'O21', 'C21', 'O22', 'C3', 'HX', 'HY', 'O31', 'C31', 'O32']


pdbfile = '../data/simulation/{}/surf_pres{}/frame{}.pdb'.format('slipids', '_30', jsk)

parser = PDBParser()
structure = parser.get_structure('model', pdbfile)
dimensions = np.zeros(3)
f = open(pdbfile, 'r')
for line in f:
    if line.startswith('CRYST1') and not np.any(dimensions):
        line_list = line.split()
        dimensions[0] = float(line_list[1])
        dimensions[1] = float(line_list[2])
        dimensions[2] = float(line_list[3])
        break
f.close()

bin_width = 1.

k = 0
for model in structure:
    k+=1

tb = np.zeros((int(np.ceil(dimensions[2] / bin_width)), k))
hb = np.zeros((int(np.ceil(dimensions[2] / bin_width)), k))
wb = np.zeros((int(np.ceil(dimensions[2] / bin_width)), k))

for i, model in enumerate(structure):
    for chain in model:
        for residue in chain:
            for atoms in residue:
                if residue.resname == 'DSP':
                    if atoms.name in heads:
                        hb[int(atoms.coord[2] / bin_width), i] += 1
                    else:
                        tb[int(atoms.coord[2] / bin_width), i] += 1
                else:
                    wb[int(atoms.coord[2] / bin_width), i] += 1
wb = wb[::-1, :]
hb = hb[::-1, :]
tb = tb[::-1, :]

wab = wb / len(waters)
hab = hb / len(heads)
tab = tb / (142-len(heads))

water_bin = wab / (dimensions[0] * dimensions[1] * bin_width)
head_bin = hab / (dimensions[0] * dimensions[1] * bin_width)
tail_bin = tab / (dimensions[0] * dimensions[1] * bin_width)

f_out = open('../output/simulation/slipids_nb{}.txt'.format(jsk), 'w')

for i in range(int(np.ceil(dimensions[2] / bin_width))):
    out = '{} {} {} {} \n'.format(i, water_bin.mean(axis=1)[i], head_bin.mean(axis=1)[i], tail_bin.mean(axis=1)[i])
    f_out.write(out)

f_out.close()

print('{} done'.format(jsk))
