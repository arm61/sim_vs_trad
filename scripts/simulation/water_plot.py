import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(palette="colorblind")
import matplotlib as mpl

mpl.rcParams["xtick.labelsize"] = 8
mpl.rcParams["ytick.labelsize"] = 8
mpl.rcParams["axes.facecolor"] = "w"
mpl.rcParams["lines.linewidth"] = 2
mpl.rcParams["xtick.top"] = False
mpl.rcParams["xtick.bottom"] = True
mpl.rcParams["ytick.left"] = True
mpl.rcParams["grid.linestyle"] = "--"
mpl.rcParams["legend.fontsize"] = 8
mpl.rcParams["legend.facecolor"] = [1, 1, 1]
mpl.rcParams["legend.framealpha"] = 0.75
mpl.rcParams["axes.labelsize"] = 8
mpl.rcParams["axes.linewidth"] = 1
mpl.rcParams["axes.edgecolor"] = "k"
mpl.rcParams["axes.titlesize"] = 8

import sys

sp = sys.argv[1]

data = np.loadtxt('../../output/simulation/waters_slipids_{}.txt'.format(sp))
def get_value(file):
    f = open(file, "r")
    for line in f:
        k = line
    if "^" in k:
        l = k.split("$")[1].split("^")[0]
    else:
        l = k.split("$")[1].split("\\pm")[0]
    return float(l)

dh = get_value('../../output/traditional/dspc_{}-d_h_{}.tex'.format(sp, sp))
phih = get_value('../../output/traditional/dspc_{}-phih_{}.tex'.format(sp, sp))
from scipy.optimize import curve_fit

def rough(x, a):
    l1 = np.array([a, 0, 0, 0])
    l2 = np.array([dh, 0.03283*phih/100, 0, 0])
    l3 = np.array([0, 0.03283, 0, 0])
    layers = np.array([l1, l2, l3])
    from scipy.stats import norm
    # work out the step in SLD at an interface
    delta_rho = layers[1:, 1] - layers[:-1, 1]
    sld = np.ones_like(x, dtype=float) * layers[0, 1]
    # the roughness of each step
    sigma = np.clip(layers[1:, 3], 1e-3, None)
    zed = np.asfarray(x)
    dist = np.cumsum(layers[:-1, 0])
    # accumulate the SLD of each step.
    for i in range(np.size(layers, 0) - 2 + 1):
        sld += delta_rho[i] * norm.cdf(zed, scale=sigma[i], loc=dist[i])
    return sld


dd, ss = curve_fit(rough, data[0], data[1], bounds=((-25), (15)))

fig = plt.figure(figsize=(5, 25/8))
plt.errorbar(data[0], data[1], yerr=data[2], marker='o')
x = np.linspace(data[0][0], data[0][-1], 10000)
plt.plot(x, rough(x, dd[0]))
plt.xlim([-20, 20])
plt.xlabel('$z$/Å')
plt.ylabel('Intrinsic density of water/Å$^{-3}$')
plt.ylim([0, 0.04])
plt.tight_layout()
plt.savefig('../../reports/figures/water_{}.pdf'.format(sp))
