import numpy as np
import os.path

import refnx, scipy

# the ReflectDataset object will contain the data
from refnx.dataset import ReflectDataset

# the reflect module contains functionality relevant to reflectometry
from refnx.reflect import ReflectModel

# the analysis module contains the curvefitting engine
from refnx.analysis import Objective, Transform, CurveFitter
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set(palette="colorblind")


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
import sim_lengths as sl
import md_simulation as md


forcefield = sys.argv[1]
surface_pressure = sys.argv[2]
lt = float(sys.argv[3])
rough = float(sys.argv[4])
traj_dir = '../../data/simulation/{}/surf_pres_{}/'.format(forcefield, surface_pressure)
anal_dir = "../../output/simulation/".format(surface_pressure)
print('{} {}'.format(forcefield, surface_pressure))

head = ['D', 'D', 'H', 'H', 'H', 'D', 'D']
tail = ['H', 'H', 'H', 'D', 'D', 'D', 'D']
sol = ['acmw', 'D', 'D', 'acmw', 'D', 'acmw', 'D']
contrasts = ['d13acmw', 'd13d2o', 'hd2o', 'd70acmw', 'd70d2o', 'd83acmw', 'd83d2o']

fig = plt.figure(figsize=(5, 25/6))
gs = mpl.gridspec.GridSpec(1, 3)
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2])
all_chi = np.array([])
abc = {'trad': '(a)', 'slipids': '(b)', 'berger': '(c)', 'martini': '(d)'}
for ci, contrast in enumerate(contrasts):
    for k in range(0, len(contrasts)):
        if contrasts[k] == contrast:
            break

    models = []
    datasets = []
    structures = []

    lgts = sl.get_lgts(head[k], tail[k], sol[k], forcefield)
    l = np.array([])
    timesteps = 0
    for i in range(1, 2):
        print('frame{}'.format(i))
        try:
            del sim
        except:
            pass
        if forcefield == 'martini':
            co = 30
        else:
            co = 15
        sim = md.MDSimulation(traj_dir + 'frame{}.pdb'.format(i), flip=True,
                              verbose=True, layer_thickness=lt, roughness=rough)

        sim.assign_scattering_lengths('neutron', atom_types=lgts[0], scattering_lengths=lgts[1])

        sim.run()
        layers_to_cut = int(co / lt) + 1
        timesteps += sim.layers.shape[0]
        l = np.append(l, sim.layers[:, :-layers_to_cut, :])

    n = l.reshape(timesteps, sim.layers.shape[1]-layers_to_cut, sim.layers.shape[2])

    data_dir = '../../data/experimental/surf_pres_{}/'.format(surface_pressure)
    dataset = ReflectDataset(os.path.join(data_dir, '{}{}.dat'.format(contrast, surface_pressure)))

    refy = np.zeros((n.shape[0], dataset.x.size))
    sldy = []
    chi = np.zeros((n.shape[0]))
    print(n.shape[0])
    for i in range(n.shape[0]):
        sim.av_layers = n[i, :, :]
        model = ReflectModel(sim)
        model.scale.setp(1, vary=True, bounds=(0.00000001, np.inf))
        model.bkg.setp(dataset.y[-1], vary=True, bounds=(dataset.y[-1]*0.9, dataset.y[-1]*1.1))
        objective = Objective(model, dataset, transform=Transform('YX4'))
        fitter = CurveFitter(objective)
        res = fitter.fit()
        refy[i] = model(dataset.x, x_err=dataset.x_err)*(dataset.x)**4
        sldy.append(sim.sld_profile()[1])
        chi[i] = objective.chisqr()
        all_chi = np.append(all_chi, objective.chisqr())
        if i == 0:
            ax1.errorbar(dataset.x,
                         dataset.y*(dataset.x)**4 * 10**(ci-1),
                         yerr=dataset.y_err*(
                             dataset.x)**4 * 10**(ci-1),
                         linestyle='', marker='o',
                         color=sns.color_palette()[ci])
        ax1.plot(dataset.x,
                 model(dataset.x,
                           x_err=dataset.x_err)*(
                     dataset.x)**4 * 10**(ci-1),
                 color=sns.color_palette()[ci], alpha=0.1)
        zs, sld = sim.sld_profile()
        if zs.min() > -20:
            x2 = np.linspace(-20, zs.min(), 100)
            zs = np.append(x2, zs)
            y2 = np.zeros_like(x2)
            sld = np.append(y2, sld)
        if zs.max() < 80:
            x3 = np.linspace(zs.max(), 81, 100)
            y3 = np.ones_like(x3) * sld[-1]
            zs = np.append(zs, x3)
            sld = np.append(sld, y3)
        ax2.plot(zs, sld + ci*10, color=sns.color_palette()[ci],
                 alpha=0.1)
        ax2.set_xlim([-20, 80])
    ax1.plot(dataset.x,
             np.average(refy, axis=0) * 10**(ci-1),
             color='k', zorder=10)
    file_open = open('{}dspc_{}_{}_{}_chi_short.txt'.format(anal_dir, forcefield, surface_pressure, contrast), 'w')
    file_open.write('{:.2f}'.format(np.average(chi)))
    file_open.close()
    print(contrast)
file_open = open('{}dspc_{}_{}_all_chi_short.txt'.format(anal_dir, forcefield, surface_pressure), 'w')
file_open.write('${:.2f}\\pm{:.2f}$'.format(np.average(all_chi), np.std(all_chi)))
file_open.close()
ax1.set_ylabel(r'$Rq^4$/Å$^{-4}$')
ax1.set_yscale('log')
ax1.set_xlabel(r'$q$/Å$^{-1}$')
ax2.set_xlabel(r'$z$/Å')
ax2.set_ylabel(r'SLD/$10^{-6}$Å$^{-2}$')
plt.tight_layout()
fig_dir = "../../reports/figures/"
plt.savefig('{}dspc_{}_{}_ref_sld_short.pdf'.format(fig_dir, forcefield, surface_pressure))
plt.close()
