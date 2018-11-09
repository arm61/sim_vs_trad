
# coding: utf-8

# In[ ]:


import numpy as np
import os.path

import refnx, scipy

# the ReflectDataset object will contain the data
from refnx.dataset import ReflectDataset

# the reflect module contains functionality relevant to reflectometry
from refnx.reflect import ReflectModel, MDSimulation

# the analysis module contains the curvefitting engine
from refnx.analysis import Objective, Transform, CurveFitter

import sys
sys.path.insert(0, '../../bin')
sys.path.insert(0, '../../models')
import sim_lengths as sl


# In[ ]:


# version numbers used in this analysis
refnx.version.version, scipy.version.version


# In[ ]:


forcefield = sys.argv[1]
surface_pressure = sys.argv[2]
contrast = sys.argv[3]
traj_dir = '../../data/simulation/{}/surf_pres_{}/'.format(forcefield, surface_pressure)
print('{} {} {}'.format(forcefield, surface_pressure, contrast))

head = ['D', 'D', 'H', 'H', 'H', 'D', 'D']
tail = ['H', 'H', 'H', 'D', 'D', 'D', 'D']
sol = ['acmw', 'D', 'D', 'acmw', 'D', 'acmw', 'D']
contrasts = ['d13acmw', 'd13d2o', 'hd2o', 'd70acmw', 'd70d2o', 'd83acmw', 'd83d2o']


# In[ ]:


for k in range(0, len(contrasts)):
    if contrasts[k] == contrast:
        break


# In[ ]:


models = []
datasets = []
structures = []


# In[ ]:


lgts = sl.get_lgts(head[k], tail[k], sol[k], forcefield)
l = np.array([])
timesteps = 0
for i in range(1, 11):
    print('frame{}'.format(i))
    try:
        del sim
    except:
        pass
    if forcefield == 'martini':
        lt = 4
        rough = 0.2
        co = 30
    else:
        lt = 1
        rough = 0.
        co = 15
    sim = MDSimulation(os.path.join(traj_dir, 'frame{}.pdb'.format(i)), flip=True, 
                       verbose=True, layer_thickness=lt, roughness=rough)

    sim.assign_scattering_lengths('neutron', atom_types=lgts[0], scattering_lengths=lgts[1])

    sim.run()
    layers_to_cut = int(co / lt) + 1
    timesteps += sim.layers.shape[0]
    l = np.append(l, sim.layers[:, :-layers_to_cut, :])


# In[ ]:


n = l.reshape(timesteps, sim.layers.shape[1]-layers_to_cut, sim.layers.shape[2])


# In[ ]:


data_dir = '../../data/experimental/surf_pres_{}/'.format(surface_pressure)
dataset = ReflectDataset(os.path.join(data_dir, '{}{}.dat'.format(contrast, surface_pressure)))


# In[ ]:


refy = np.zeros((n.shape[0], dataset.x.size))
sldy = []
chi = np.zeros((n.shape[0]))


# In[ ]:


for i in range(n.shape[0]):
    sim.av_layers = n[i, :, :]
    model = ReflectModel(sim)
    model.scale.setp(100, vary=True, bounds=(0.00000001, np.inf))
    model.bkg.setp(dataset.y[-1], vary=False)
    objective = Objective(model, dataset, transform=Transform('YX4'))
    fitter = CurveFitter(objective)
    res = fitter.fit()
    refy[i] = model(dataset.x, x_err=dataset.x_err)*(dataset.x)**4
    sldy.append(sim.sld_profile()[1])
    chi[i] = objective.chisqr()


# In[ ]:


import scipy
from scipy.stats.mstats import mquantiles
chi_arr = mquantiles(chi, prob=[0.025, 0.5, 0.975])
chia = '{:.0f}'.format(chi_arr[1])

chi_out = open('../../output/simulation/{}_{}_{}_chisq.txt'.format(contrast, forcefield, surface_pressure), 'w')
chi_out.write('$' + str(chia) + '$')
chi_out.close()


# In[ ]:


ref_out = open('../../output/simulation/{}_{}_{}_ref.txt'.format(contrast, forcefield, surface_pressure), 'w')
for i in range(dataset.x.size):
    ref_out.write('{} {} {} '.format(dataset.x[i], dataset.y[i]*(dataset.x[i])**4, dataset.y_err[i]*(dataset.x[i])**4))
    for j in range(n.shape[0]):
        ref_out.write('{} '.format(refy[j, i]))
    ref_out.write('\n')
ref_out.close()


# In[ ]:


sld_out = open('../../output/simulation/{}_{}_{}_sld.txt'.format(contrast, forcefield, surface_pressure), 'w')
sldya = np.asarray(sldy)
for i in range(sim.sld_profile()[0].shape[0]):
    sld_out.write('{} '.format(sim.sld_profile()[0][i]))
    for j in range(sldya.shape[0]):
        sld_out.write('{} '.format(sldya[j, i]))
    sld_out.write('\n')
sld_out.close()

