import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns; sns.set('paper', palette='colorblind')
import pandas as pd

data = pd.read_csv('../data/surf_iso.csv')

data_sorted = data.sort_values(['apm'])

mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.facecolor'] = 'w'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['xtick.top'] = False
mpl.rcParams['xtick.bottom'] = True
mpl.rcParams['ytick.left'] = True
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.facecolor'] = [1,1,1]
mpl.rcParams['legend.framealpha'] = 0.75
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.edgecolor'] = 'k'

plt.figure(figsize=(5, 25/11))
plt.plot(data_sorted['apm'], data_sorted['sp'], c='#0173B2')
plt.ylim([0, 80])
#plt.xticks([40, 50, 60, 70, 80])
plt.ylabel('$\pi$/mNm$^{-1}$')
plt.xlabel('Area per molecule/Å$^{2}$')
plt.tight_layout()
plt.savefig('../reports/figures/apm.pdf')
#plt.show()
