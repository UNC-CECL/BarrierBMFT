"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 19 January 2022
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
SimDur = 100  # [Yr] Duration of each simulation

# Parameter values
rslr = [1, 3, 5, 7, 9]
co = [15, 30, 45, 60, 75]
slope = [0.005]

SimNum = len(rslr) * len(co) * len(slope)

# ==================================================================================================================================================================================
# Load data
filename = '/Users/ianreeves/PycharmProjects/BarrierBMFT/Output/Batch_2022_0118_17_15/'

BarrierWidth = np.load(filename + 'Widths_Barrier.npy')
BarrierWidth = BarrierWidth[0, :, :, 0] * -1
BarrierWidth = np.rot90(BarrierWidth)
BarrierWidth = np.flipud(BarrierWidth)

BBMarshWidth = np.load(filename + 'Widths_BBMarsh.npy')
BBMarshWidth = BBMarshWidth[0, :, :, 0] * -1
BBMarshWidth = np.rot90(BBMarshWidth)
BBMarshWidth = np.flipud(BBMarshWidth)

BayWidth = np.load(filename + 'Widths_Bay.npy')
BayWidth = BayWidth[0, :, :, 0] * -1
BayWidthOrig = BayWidth
BayWidth = np.rot90(BayWidth)
BayWidth = np.flipud(BayWidth)

MLMarshWidth = np.load(filename + 'Widths_MLMarsh.npy')
MLMarshWidth = MLMarshWidth[0, :, :, 0] * -1
MLMarshWidth = np.rot90(MLMarshWidth)
MLMarshWidth = np.flipud(MLMarshWidth)

ForestWidth = np.load(filename + 'Widths_Forest.npy')
ForestWidth = ForestWidth[0, :, :, 0] * -1
ForestWidthOrig = ForestWidth
ForestWidth = np.rot90(ForestWidth)
ForestWidth = np.flipud(ForestWidth)

ShorelineChange = np.load(filename + 'ShorelineChange.npy')
ShorelineChange = ShorelineChange[0, :, :, 0]
ShorelineChange = np.rot90(ShorelineChange)
ShorelineChange = np.flipud(ShorelineChange)


# SimEl = np.load(filename + 'Sim25_elevation.npy', allow_pickle=True)

# ==================================================================================================================================================================================
# Plot

xtic = ['', '1', '3', '5', '7', '9']
ytic = ['', '15', '30', '45', '60', '75']

xlab = 'RSLR [mm/yr]'
ylab = 'Back-Barrier SSC [mg/L]'

all = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth))
vmin = int(np.min(all))
vmax = int(np.max(all))
maximum = max(abs(vmin), abs(vmax))
vmax = maximum
vmin = -maximum


cmap = 'RdBu'
Fig = plt.figure(figsize=(28, 6))
plt.rcParams.update({'font.size': 14, 'font.family': 'Arial'})

ax = Fig.add_subplot(151)
cax = ax.matshow(BarrierWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab)
plt.title('Barrier')


ax = Fig.add_subplot(152)
cax = ax.matshow(BBMarshWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Back-Barrier Marsh')

ax = Fig.add_subplot(153)
cax = ax.matshow(BayWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.title('Bay')

ax = Fig.add_subplot(154)
cax = ax.matshow(MLMarshWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Mainland Marsh')

ax = Fig.add_subplot(155)
cax = ax.matshow(ForestWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
cbar = Fig.colorbar(cax)
cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Forest')
plt.tight_layout()

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


# plt.figure(figsize=(12, 6))
# plt.xlabel("Distance Cross-Shore [m]")
# plt.ylabel("Elevation [m]")
# plt.title("Sim 25")
# for t in range(0, SimDur, 10):
#     elev = SimEl[t]
#     plt.plot(elev)


plt.show()

