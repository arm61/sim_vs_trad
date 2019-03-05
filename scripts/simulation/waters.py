import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(palette="colorblind")
import matplotlib as mpl

def xydist(a, b, cell):
    dx = b[0]-a[0]
    dx = dx % cell[0]
    dy = b[1]-a[1]
    dy = dy % cell[1]
    dist = np.sqrt(np.square(dx) + np.square(dy))
    return dist

import sys

ff = sys.argv[1]
sp = sys.argv[2]

aks = 1
sda = np.array([])
for f in range(1, 11):
    print('frame{}'.format(f))
    u = mda.Universe('../../data/simulation/{}/surf_pres_{}/frame{}.pdb'.format(ff, sp, f))
    p_pos = np.zeros((len(u.trajectory), 100, 3))
    for k, ts in enumerate(u.trajectory):
        print(ts)
        count = 0
        count_w = 0
        for i in u.atoms:
            if i.name == 'P':
                p_pos[k, count, :] = i.position
                count += 1
            if i.name == 'OW':
                count_w += 1

    nearest = np.zeros((count_w, len(u.trajectory)), dtype=int)
    for k, ts in enumerate(u.trajectory):
        print(ts)
        count = 0
        for i in u.atoms:
            if i.name == 'OW':
                best = 0
                dist = ts.dimensions[0]
                for j in range(0, p_pos.shape[1]):
                    new_d = xydist(i.position, p_pos[k, j], ts.dimensions)
                    if new_d < dist:
                        dist = new_d
                        best = j
                nearest[count, k] = best
                count += 1
        count = 0
        for i in u.atoms:
            if i.name == 'OW':
                a = p_pos[k, nearest[count, k]][2] - i.position[2]
                b = i.position[2]
                c = p_pos[k, nearest[count, k]][2]
                d = a - b + c
                sda = np.append(sda, d)
                count += 1
        sda = np.reshape(sda, (aks, count_w))
        aks += 1

jj = len(np.arange(int(np.floor(np.min(sda))), int(np.ceil(np.max(sda))), 1))
ks = np.zeros((sda.shape[0], jj))
for dd, sd in enumerate(sda):
    z = np.arange(int(np.floor(np.min(sda))), np.ceil(np.max(sda)), 1)
    jk = int(np.floor(np.min(sda))) / 1
    for i in sd:
        ks[dd, int(i-jk)] += 1
water = ks / (0.5 * ts.dimensions[0] * ts.dimensions[1])

outarray = np.array([z+0.5, np.average(water, axis=0), np.std(water, axis=0)])
np.savetxt('../../output/simulation/waters_{}_{}.txt'.format(ff, sp), outarray)
