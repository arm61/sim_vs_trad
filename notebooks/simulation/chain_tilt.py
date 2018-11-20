
# coding: utf-8

# In[1]:


import numpy as np
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set(palette='colorblind')
import corner
mpl.rcParams['axes.labelsize']=28
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18


# In[2]:


import sys
from scipy.stats.mstats import mquantiles

ff = sys.argv[1]

sp = sys.argv[2]

if ff == 'berger':
    chain1 = ['C15', 'C33']
    chain2 = ['C36', 'C54']
elif ff == 'martini':
    chain1 = ['GL1', 'C4A']
    chain2 = ['GL2', 'C4B']
elif ff == 'slipids':
    chain1 = ['C21', '8C21']
    chain2 = ['C31', '8C31']
if ff == 'berger':
    heads = ['C1', 'C2', 'C3', 'N4', 'C5', 'C6', 'O7', 'P8', 'O9', 'O10', 'O11', 
             'C12', 'C13', 'O14', 'C15', 'O16', 'C34', 'O35', 'C36', 'O37']
    waters = ['OW', 'HW1', 'HW2']
if ff == 'martini':
    heads = ['NC3', 'PO4', 'GL1', 'GL2']
    waters = ['W', 'WP', 'WM']
if ff == 'slipids':
    waters = ['OW', 'HW1', 'HW2']
    heads = ['N', 'C13', 'H13A', 'H13B', 'H13C', 'C14', 'H14A', 'H14B', 'H14C',
           'C15', 'H15A', 'H15B', 'H15C', 'C12', 'H12A', 'H12B', 'C11', 'H11A',
           'H11B', 'P', 'O13', 'O14', 'O11', 'O12', 'C1', 'HA', 'HB', 'C2', 'HS',
           'O21', 'C21', 'O22', 'C3', 'HX', 'HY', 'O31', 'C31', 'O32']
at = np.array([])
al = np.array([])
dh_tot = np.array([])
wph_tot = np.array([])
for i in range(1, 2):
    print(ff, sp, i)
    pdbfile = '../../data/simulation/{}/surf_pres{}/frame{}.pdb'.format(ff, sp, i)
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
    all_theta = np.array([])
    all_length = np.array([])
    labels = [chain1, chain2]
    for model in structure:
        for chain in model:
            for residue in chain:
                for i in labels:
                    if residue.resname != 'SOL' and residue.resname != ' PW':
                        names = ['a', 'a']
                        tail_atoms = np.zeros((2,3))
                        k = 0
                        for atoms in residue:
                            if atoms.name in i:
                                tail_atoms[k, :] = atoms.coord
                                names[k] = atoms.name
                                k += 1
                        d = np.abs(tail_atoms[1, :] - tail_atoms[0, :])
                        for j in range(0, 3):
                            if d[j] > dimensions[j] / 2.:
                                d[j] = np.abs(d[j] - dimensions[j])
                        xy = np.sqrt(np.square(d[0]) + np.square(d[1]))
                        theta = np.arctan(xy / d[2])
                        dist = np.sqrt(np.square(d[0]) + np.square(d[1]) + np.square(d[2]))
                        all_length = np.append(all_length, dist)
                        all_theta = np.append(all_theta, np.rad2deg(theta))
    at = np.append(at, all_theta)
    al = np.append(al, all_length)
    bin_width = 0.01

    k = 0
    for model in structure:
        k+=1 

    hb = np.zeros((int(np.ceil(dimensions[2] / bin_width)), k))
    wb = np.zeros((int(np.ceil(dimensions[2] / bin_width)), k))


    for i, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                if residue.resname != 'SOL' and residue.resname != ' PW':
                    for atoms in residue:
                        if atoms.name in heads:
                            hb[int(atoms.coord[2] / bin_width), i] += 1
                else:
                    for atoms in residue:
                        if atoms.name in waters:
                            wb[int(atoms.coord[2] / bin_width), i] += 1

    wb = wb[::-1, :]
    hb = hb[::-1, :]

    wab = wb / len(waters)
    hab = hb / len(heads)

    water_bin = wab / (dimensions[0] * dimensions[1] * bin_width)
    head_bin = hab / (dimensions[0] * dimensions[1] * bin_width)

    total_head = np.sum(head_bin, axis=0)

    dh = np.zeros(water_bin.shape[1])
    e = np.zeros(water_bin.shape[1])
    s = np.zeros(water_bin.shape[1])

    for j in range(0, water_bin.shape[1]):
        summing = 0
        start = 0
        end = 0
        for i in range(0, water_bin.shape[0]):
            if summing > total_head[j] * 0.2:
                start = i * bin_width
                break
            else:
                summing += head_bin[i, j]

        summing = 0
        for i in range(0, water_bin.shape[0]):
            if summing > total_head[j] * 0.8:
                end = i * bin_width
                break
            else:
                summing += head_bin[i, j]

        dh[j] = end-start
        e[j] = end
        s[j] = start
    dh_tot = np.append(dh_tot, dh)

    a = mquantiles(dh, prob=[0.025, 0.5, 0.975])
    a

    wph = np.zeros(water_bin.shape[1])
    for j in range(0, water_bin.shape[1]): 
        w = 0
        h = 0
        for i in np.arange(0, int(np.ceil(dimensions[2] / bin_width))):
            if i * bin_width > s[j] and i * bin_width <= e[j] + bin_width:
                w += wab[i, j]
                h += hab[i, j]
                #print(wab[j, 0], hab[j, 0])
        wph[j] = w/h
    wph_tot = np.append(wph_tot, wph)

weightal = np.ones_like(wph_tot)/float(len(wph_tot))
plt.hist(wph_tot, bins=15, histtype='stepfilled', weights=weightal)
plt.xlabel(r'wph')
plt.ylabel(r'PDF(wph)')
plt.tight_layout()
plt.savefig('../../reports/figures/{}{}_wph.pdf'.format(ff, sp))
plt.close()
f_out  = open('../../output/simulation/{}{}_wph.txt'.format(ff, sp), 'w')
a = mquantiles(wph_tot, prob=[0.025, 0.5, 0.975])
k = [a[1], a[1] - a[0], a[2] - a[1]]
q = '{:.2f}'.format(k[0])
e = '{:.2f}'.format(k[1])
w = '{:.2f}'.format(k[2])
f_out.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
f_out.close()
weightal = np.ones_like(dh_tot)/float(len(dh_tot))
plt.hist(dh_tot, bins=15, histtype='stepfilled', weights=weightal)
plt.xlabel(r'$d_h$/Å')
plt.ylabel(r'PDF($d_h$)')
plt.tight_layout()
plt.savefig('../../reports/figures/{}{}_dh.pdf'.format(ff, sp))
plt.close()
f_out  = open('../../output/simulation/{}{}_dh.txt'.format(ff, sp), 'w')
a = mquantiles(dh_tot, prob=[0.025, 0.5, 0.975])
k = [a[1], a[1] - a[0], a[2] - a[1]]
q = '{:.2f}'.format(k[0])
e = '{:.2f}'.format(k[1])
w = '{:.2f}'.format(k[2])
f_out.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
f_out.close()

weightal = np.ones_like(al)/float(len(al))
plt.hist(al, bins=25, histtype='stepfilled', weights=weightal)
plt.xlabel(r'$t_t$/Å')
plt.ylabel(r'PDF($t_t$)')
plt.tight_layout()
plt.savefig('../../reports/figures/{}{}_tt.pdf'.format(ff, sp))
plt.close()
f_out  = open('../../output/simulation/{}{}_tt.txt'.format(ff, sp), 'w')
a = mquantiles(al, prob=[0.025, 0.5, 0.975])
k = [a[1], a[1] - a[0], a[2] - a[1]]
q = '{:.2f}'.format(k[0])
e = '{:.2f}'.format(k[1])
w = '{:.2f}'.format(k[2])
f_out.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
f_out.close()

weightat = np.ones_like(at)/float(len(at))
plt.hist(at, bins=25, histtype='stepfilled', weights=weightat)
plt.xlabel(r'$\theta_t$/$^\circ$')
plt.ylabel(r'PDF($\theta_t$)')
plt.tight_layout()
plt.savefig('../../reports/figures/{}{}_angle.pdf'.format(ff, sp))
plt.close()
f_out = open('../../output/simulation/{}{}_angle.txt'.format(ff, sp), 'w')
a = mquantiles(at, prob=[0.025, 0.5, 0.975])
k = [a[1], a[1] - a[0], a[2] - a[1]]
q = '{:.2f}'.format(k[0])
e = '{:.2f}'.format(k[1])
w = '{:.2f}'.format(k[2])
f_out.write('$' + str(q) + '^{+' + str(w) + '}_{-' + str(e) + '}$')
f_out.close()

