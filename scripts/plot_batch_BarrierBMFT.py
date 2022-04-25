"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 25 April 2022
"""

import numpy as np
import matplotlib.pyplot as plt

# ==================================================================================================================================================================================
# Define batch parameters

SimDur = 400  # [Yr] Duration of each simulation

# Parameter values
rslr = [3, 6, 9, 12, 15]
co = [40, 50, 60, 70, 80]
slope = [0.001]

SimNum = len(rslr) * len(co) * len(slope)

Sim_plot = 5  # Simulation number for plotting
Data_Sim_plot = 1  # Data file number for plotting

plot_param_space = True
plot_elev = True
plot_minus = False
plot_dune_width = True

# ==================================================================================================================================================================================
# Specify data

# ~ File bank ~
# filename = '/BarrierBMFT/Output/Batch_2022_0207_13_37/'
# filename = '/Users/reevesi/PycharmProjects/BarrierBMFT/Output/Batch_2022_0408_13_55/'
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0408_19_32/'
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0423_15_12/'  # 38
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0423_15_19/'  # 44
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_28/'  # 49
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_30/'  # 51
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_31/'  # 52
# filename = '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_32/'  # 53

data_files = [
    '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_30/',  # 51
    '/Users/reevesi/DesktopBackup/BarrierBMFT/Data/Batch_2022_0425_00_31/',  # 52
]

file_num = len(data_files)

# Initialize data arrays
BarrierWidth = np.zeros([file_num, len(rslr), len(co)])
BBMarshWidth = np.zeros([file_num, len(rslr), len(co)])
BayWidth = np.zeros([file_num, len(rslr), len(co)])
MLMarshWidth = np.zeros([file_num, len(rslr), len(co)])
ForestWidth = np.zeros([file_num, len(rslr), len(co)])
BBPondWidth = np.zeros([file_num, len(rslr), len(co)])
MLPondWidth = np.zeros([file_num, len(rslr), len(co)])
ShorelineChange = np.zeros([file_num, len(rslr), len(co)])

# ==================================================================================================================================================================================
# Load data

for n in range(file_num):
    filename = data_files[n]

    nBarrierWidth = np.load(filename + 'Widths_Barrier.npy')
    nBarrierWidth = np.mean(nBarrierWidth[:, :, :, 0], axis=0)
    nBarrierWidth = np.rot90(nBarrierWidth)
    BarrierWidth[n, :, :] = np.flipud(nBarrierWidth)

    nBBMarshWidth = np.load(filename + 'Widths_BBMarsh.npy')
    nBBMarshWidth = np.mean(nBBMarshWidth[:, :, :, 0], axis=0)
    nBBMarshWidth = np.rot90(nBBMarshWidth)
    BBMarshWidth[n, :, :] = np.flipud(nBBMarshWidth)

    nBayWidth = np.load(filename + 'Widths_Bay.npy')
    nBayWidth = np.mean(nBayWidth[:, :, :, 0], axis=0)
    nBayWidth = np.rot90(nBayWidth)
    BayWidth[n, :, :] = np.flipud(nBayWidth)

    nMLMarshWidth = np.load(filename + 'Widths_MLMarsh.npy')
    nMLMarshWidth = np.mean(nMLMarshWidth[:, :, :, 0], axis=0)
    nMLMarshWidth = np.rot90(nMLMarshWidth)
    MLMarshWidth[n, :, :] = np.flipud(nMLMarshWidth)

    nForestWidth = np.load(filename + 'Widths_Forest.npy')
    nForestWidth = np.mean(nForestWidth[:, :, :, 0], axis=0)
    nForestWidth = np.rot90(nForestWidth)
    ForestWidth[n, :, :] = np.flipud(nForestWidth)

    nBBPondWidth = np.load(filename + 'Widths_BBMarshPond.npy')
    nBBPondWidth = np.mean(nBBPondWidth[:, :, :, 0], axis=0)
    nBBPondWidth = np.rot90(nBBPondWidth)
    BBPondWidth[n, :, :] = np.flipud(nBBPondWidth)

    nMLPondWidth = np.load(filename + 'Widths_MLMarshPond.npy')
    nMLPondWidth = np.mean(nMLPondWidth[:, :, :, 0], axis=0)
    nMLPondWidth = np.rot90(nMLPondWidth)
    MLPondWidth[n, :, :] = np.flipud(nMLPondWidth)

    nShorelineChange = np.load(filename + 'ShorelineChange.npy')
    nShorelineChange = np.mean(nShorelineChange[:, :, :, 0], axis=0)
    nShorelineChange = nShorelineChange * 10  # Convert to m
    nShorelineChange = np.rot90(nShorelineChange)
    ShorelineChange[n, :, :] = np.flipud(nShorelineChange)

# Average across data files
BarrierWidth = np.mean(BarrierWidth, axis=0)
BBMarshWidth = np.mean(BBMarshWidth, axis=0)
BayWidth = np.mean(BayWidth, axis=0)
MLMarshWidth = np.mean(MLMarshWidth, axis=0)
ForestWidth = np.mean(ForestWidth, axis=0)
BBPondWidth = np.mean(BBPondWidth, axis=0)
MLPondWidth = np.mean(MLPondWidth, axis=0)
ShorelineChange = np.mean(ShorelineChange, axis=0)


# ==================================================================================================================================================================================
# Plot

if plot_param_space:

    xtic = ['', '3', '6', '9', '12', '15']
    ytic = ['', '40', '50', '60', '70', '80']

    xlab = 'RSLR [mm/yr]'
    ylab = 'Back-Barrier SSC [mg/L]'

    sumwidth = BarrierWidth + BBMarshWidth + BayWidth + MLMarshWidth + ForestWidth
    all_widths = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth))
    # all_widths = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth, ShorelineChange))
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
    # cax = ax.matshow(sumwidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    cax = ax.matshow(ShorelineChange, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)
    # plt.title('Total Landscape Width')
    plt.title('Ocean Shoreline Change')
    plt.tight_layout()

    # --------------------------
    # Marsh + Pond
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
    cax = ax.matshow(BBMarshWidth + BBPondWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
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
    cax = ax.matshow(MLMarshWidth + MLPondWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
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
    # cax = ax.matshow(sumwidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    cax = ax.matshow(ShorelineChange, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)
    # plt.title('Total Landscape Width')
    plt.title('Ocean Shoreline Change')
    plt.tight_layout()

# --------------------------
if plot_elev:
    SimEl = np.load(data_files[Data_Sim_plot - 1] + 'Sim' + str(Sim_plot) + '_elevation.npy', allow_pickle=True)
    plt.figure(figsize=(12, 6))
    plt.xlabel("Distance Cross-Shore [m]")
    plt.ylabel("Elevation [m]")
    plt.title("Sim" + str(Sim_plot))
    for t in range(0, len(SimEl), 50):
        elev = SimEl[t]
        plt.plot(elev)

# --------------------------
if plot_minus:

    # Subtract forest width change from ML marsh change
    Fig = plt.figure(figsize=(8, 8))
    ax = Fig.add_subplot(111)
    plt.xlabel("RSLR [mm/yr]")
    plt.ylabel("Back-Barrier SSC [mg/L]")
    plt.title("Mainland Marsh (Minus Forest Change)")
    cax = ax.matshow(MLMarshWidth + ForestWidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)
    plt.tight_layout()

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

# --------------------------
if plot_dune_width:

    SimWidths = np.load(data_files[Data_Sim_plot - 1] + 'Sim' + str(Sim_plot) + '_widthsTS.npy', allow_pickle=True)
    SimStats = np.load(data_files[Data_Sim_plot - 1] + 'Sim' + str(Sim_plot) + '_statsTS.npy', allow_pickle=True)

    aHd = SimStats[:, 4]

    widths = SimWidths
    total_width = np.sum(widths, axis=1)
    barrier = widths[:, 0]
    BBmarsh = widths[:, 1]
    bay = widths[:, 2]
    MLmarsh = widths[:, 3]
    forest = widths[:, 4]
    BBpond = widths[:, 5]
    MLpond = widths[:, 6]
    barrier_marsh = barrier + BBmarsh

    # ------
    fig, ax1 = plt.subplots()
    fig.set_size_inches(12, 9)

    ax1.set_xlabel('Year')
    ax1.set_ylabel('BB Marsh Width [m]', color='red')
    ax1.plot(BBmarsh, c="red")
    ax1.tick_params(axis='y', labelcolor='red')

    ax2 = ax1.twinx()
    ax2.set_ylabel("Dune Height [m]", color='black')
    ax2.plot(aHd, c="black")
    ax2.tick_params(axis='y', labelcolor='black')


# --------------------------
plt.show()

