"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 7 February 2022
"""

import numpy as np
import matplotlib.pyplot as plt

# ==================================================================================================================================================================================
# Define batch parameters

Num = 5  # Number of runs at each combinations of parameter values
SimDur = 250  # [Yr] Duration of each simulation

# Parameter values
rslr = [3, 6, 9, 12, 15]
co = [20, 30, 40, 50, 60]
slope = [0.005]

SimNum = len(rslr) * len(co) * len(slope)

# ==================================================================================================================================================================================
# Load data
# filename = '/BarrierBMFT/Output/Batch_2022_0207_13_37/'
filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0301_21_51/'
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0302_23_57/'
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0303_11_59/'

BarrierWidth = np.load(filename + 'Widths_Barrier.npy')
BarrierWidth = np.mean(BarrierWidth[:, :, :, 0], axis=0)
BarrierWidth = np.rot90(BarrierWidth)
BarrierWidth = np.flipud(BarrierWidth)

BBMarshWidth = np.load(filename + 'Widths_BBMarsh.npy')
BBMarshWidth = np.mean(BBMarshWidth[:, :, :, 0], axis=0)
BBMarshWidth = np.rot90(BBMarshWidth)
BBMarshWidth = np.flipud(BBMarshWidth)

BayWidth = np.load(filename + 'Widths_Bay.npy')
BayWidth = np.mean(BayWidth[:, :, :, 0], axis=0)
BayWidthOrig = BayWidth
BayWidth = np.rot90(BayWidth)
BayWidth = np.flipud(BayWidth)

MLMarshWidth = np.load(filename + 'Widths_MLMarsh.npy')
MLMarshWidth = np.mean(MLMarshWidth[:, :, :, 0], axis=0)
MLMarshWidth = np.rot90(MLMarshWidth)
MLMarshWidth = np.flipud(MLMarshWidth)

ForestWidth = np.load(filename + 'Widths_Forest.npy')
ForestWidth = np.mean(ForestWidth[:, :, :, 0], axis=0)
ForestWidthOrig = ForestWidth
ForestWidth = np.rot90(ForestWidth)
ForestWidth = np.flipud(ForestWidth)

ShorelineChange = np.load(filename + 'ShorelineChange.npy')
ShorelineChange = np.mean(ShorelineChange[:, :, :, 0], axis=0)
ShorelineChange = ShorelineChange * 10  # Convert to m
ShorelineChange = np.rot90(ShorelineChange)
ShorelineChange = np.flipud(ShorelineChange)

# SimEl = np.load(filename + 'Sim25_elevation.npy', allow_pickle=True)

# ==================================================================================================================================================================================
# Plot

xtic = ['', '3', '6', '9', '12', '15']
ytic = ['', '20', '30', '40', '50', '60']

xlab = 'RSLR [mm/yr]'
ylab = 'Back-Barrier SSC [mg/L]'

sumwidth = BarrierWidth + BBMarshWidth + BayWidth + MLMarshWidth + ForestWidth
all_widths = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth))
maximum = max(abs(int(np.min(all_widths))), abs(int(np.max(all_widths))))
vmax = maximum
vmin = -maximum


cmap = 'RdBu'
Fig = plt.figure(figsize=(16, 3.4))
plt.rcParams.update({'font.size': 10, 'font.family': 'Arial'})

ax = Fig.add_subplot(161)
cax = ax.matshow(BarrierWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab)
plt.title('Barrier')

ax = Fig.add_subplot(162)
cax = ax.matshow(BBMarshWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Back-Barrier Marsh')

ax = Fig.add_subplot(163)
cax = ax.matshow(BayWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.title('Bay')

ax = Fig.add_subplot(164)
cax = ax.matshow(MLMarshWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Mainland Marsh')

ax = Fig.add_subplot(165)
cax = ax.matshow(ForestWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Forest')
plt.tight_layout()

ax = Fig.add_subplot(166)
cax = ax.matshow(sumwidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
# cax = ax.matshow(ShorelineChange, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title('Total Landscape Width')
# plt.title('Ocean Shoreline Change')
plt.tight_layout()

# --------------------------

# plt.figure(figsize=(12, 6))
# plt.xlabel("Distance Cross-Shore [m]")
# plt.ylabel("Elevation [m]")
# plt.title("Sim 25")
# for t in range(0, SimDur, 10):
#     elev = SimEl[t]
#     plt.plot(elev)

# --------------------------
# Subtract forest width change from ML marsh change
# Fig = plt.figure(figsize=(8, 8))
# ax = Fig.add_subplot(111)
# plt.xlabel("RSLR [mm/yr]")
# plt.ylabel("Back-Barrier SSC [mg/L]")
# plt.title("Mainland Marsh (Minus Forest Change)")
# cax = ax.matshow(MLMarshWidth + ForestWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
# ax.xaxis.set_ticks_position('bottom')
# # cbar = Fig.colorbar(cax)
# # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.tight_layout()

# Subtract total landscape width change from bay change
Fig = plt.figure(figsize=(8, 8))
ax = Fig.add_subplot(111)
plt.xlabel("RSLR [mm/yr]")
plt.ylabel("Back-Barrier SSC [mg/L]")
plt.title("Bay (Minus Total Landscape Change)")
cax = ax.matshow(BayWidth + sumwidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
ax.xaxis.set_ticks_position('bottom')
# cbar = Fig.colorbar(cax)
# cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.tight_layout()


plt.show()

