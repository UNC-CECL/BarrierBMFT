"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 6 January 2022
"""

import time
import math
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from datetime import datetime

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools import plot as B3Dfunc

# ==================================================================================================================================================================================
# Define batch parameters

Num = 1  # Number of runs at each combinations of parameter values
SimDur = 125  # [Yr] Duration of each simulation

# Parameter values
# rslr = [2, 6, 10, 14, 18]
# co = [10, 30, 50, 70, 90]
# slope = [0.001, 0.003, 0.005, 0.007, 0.009]
rslr = [2, 10, 18]
co = [10, 50, 90]
slope = [0.003]

SimNum = len(rslr) * len(co) * len(slope)

# ==================================================================================================================================================================================
# Load data
filename = '/Users/reevesi/PycharmProjects/BarrierBMFT/Output/Batch_2022_0107_18_39/'

BarrierWidth = np.load(filename + 'Widths_Barrier.npy')
BarrierWidth = BarrierWidth[0, :, :, 0]

BBMarshWidth = np.load(filename + 'Widths_BBMarsh.npy')
BBMarshWidth = BBMarshWidth[0, :, :, 0]

BayWidth = np.load(filename + 'Widths_Bay.npy')
BayWidth = BayWidth[0, :, :, 0]

MLMarshWidth = np.load(filename + 'Widths_MLMarsh.npy')
MLMarshWidth = MLMarshWidth[0, :, :, 0]

ForestWidth = np.load(filename + 'Widths_Forest.npy')
ForestWidth = ForestWidth[0, :, :, 0]

ShorelineChange = np.load(filename + 'ShorelineChange.npy')
ShorelineChange = ShorelineChange[0, :, :, 0]


# ==================================================================================================================================================================================
# Plot

xtic = ['', '2', '10', '18']
ytic = ['', '10', '50', '90']

xlab = 'RSLR [mm/yr]'
ylab = 'Susp. Sed. Concentration [mg/L]'

all = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth))
vmin = int(np.min(all))
vmax = int(np.max(all))

Fig = plt.figure(figsize=(13.5, 3))
plt.rcParams.update({'font.size': 10})

ax = Fig.add_subplot(151)
ax.matshow(BarrierWidth, origin='lower', cmap='Greys', aspect='auto')  # , vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab)
plt.title('Barrier')

ax = Fig.add_subplot(152)
ax.matshow(BBMarshWidth, origin='lower', cmap='Greys', aspect='auto')  # , vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Back-Barrier Marsh')

ax = Fig.add_subplot(153)
ax.matshow(BayWidth, origin='lower', cmap='Greys', aspect='auto')  # , vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.title('Bay')

ax = Fig.add_subplot(154)
ax.matshow(MLMarshWidth, origin='lower', cmap='Greys', aspect='auto')  # , vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Mainland Marsh')

ax = Fig.add_subplot(155)
cax = ax.matshow(ForestWidth, origin='lower', cmap='Greys', aspect='auto')  # , vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
cbar = Fig.colorbar(cax)
cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Forest')

# sumplot = BarrierWidth + BBMarshWidth + BayWidth + MLMarshWidth + ForestWidth
# ax = Fig.add_subplot(166)
# cax = ax.matshow(sumplot, origin='lower', cmap='Greys', aspect='auto')
# ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.title('')

plt.tight_layout()
print()
plt.show()

