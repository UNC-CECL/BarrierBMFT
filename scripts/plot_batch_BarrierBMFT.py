"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 25 April 2022
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ==================================================================================================================================================================================
# Define batch parameters

SimDur = 400  # [Yr] Duration of each simulation

# Parameter values
rslr = [3, 6, 9, 12, 15]
co = [40, 50, 60, 70, 80]
slope = [0.005]

SimNum = len(rslr) * len(co) * len(slope)

Sim_plot = 2  # Simulation number for plotting
Data_Sim_plot = 2  # Data file number for plotting

plot_param_space = True
plot_param_space_pond = False
plot_elev = False
plot_minus = False
plot_dune_width = False

title = 'Slow Dune Growth, Moderate Slope'

# ==================================================================================================================================================================================
# Specify data

data_files = [
    # Group 4: Fast, shallow
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_28/',  # 49
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_34/',  # 71
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_00_33/',  # 77
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_17/',  # 84
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_31/',  # 102
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0505_18_51/',  # 114

    # Group 5: Fast, moderate
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_30/',  # 51
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_31/',  # 52
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_29/',  # 57
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_35/',  # 72
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_36/',  # 73
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_00_36/',  # 79
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_00_37/',  # 80
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_20/',  # 86
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0511_00_09/',  # 116
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_11/',  # 119
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_14/',  # 120

    # Group 6: Fast, steep
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_32/',  # 53
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_30/',  # 58
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_31/',  # 59
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_38/',  # 74
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_39/',  # 75
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_00_39/',  # 82
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_22/',  # 87

    # Group 7: Slow, shallow
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_21/',  # 61
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_00_08/',  # 69
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_37/',  # 70
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_18_27/',  # 92

    # Group 8: Slow, moderate
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_22/',  # 62
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_23/',  # 63
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_01/',  # 66
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_16_14/',  # 100
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0511_18_35/',  # 118
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_16/',  # 121
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_17/',  # 122

    # Group 9: Slow, steep
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_24/',  # 64
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_25/',  # 65
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0428_18_31/',  # 76
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_16_37/',  # 101

    # Group 10: Fast, moderate, Rinr=1
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_28/',  # 89
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_11_36/',  # 94
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_12_09/',  # 98
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_13_52/',  # 99
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0502_17_28/',  # 110

    # Group 11: Fast, moderate, Rinr=2
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_29/',  # 90
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_12_30/',  # 91
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_34/',  # 104
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_35/',  # 105
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_36/',  # 106
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0505_00_20/',  # 112

    # Group 12: Moderate, moderate
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_11_42/',  # 95
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_11_43/',  # 96
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_11_46/',  # 97
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_38/',  # 107
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_23_39/',  # 108
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0505_00_22/',  # 113

    # Group 13: Fast, moderate, Cbb=0.5/0.65
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_23_40/',  # 123
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0516_12_30/',  # 124
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0518_13_05/',  # 125
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

    xlab = 'RSLR (mm/yr)'
    ylab = 'External SSC (mg/L)'

    sumwidth = BarrierWidth + BBMarshWidth + BayWidth + MLMarshWidth + ForestWidth
    all_widths = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth))
    # all_widths = np.concatenate((BarrierWidth, BBMarshWidth, BayWidth, MLMarshWidth, ForestWidth, ShorelineChange))
    maximum = max(abs(int(np.min(all_widths))), abs(int(np.max(all_widths))))
    vmax = maximum
    vmin = -maximum
    # vmax = 1300
    # vmin = -1300
    print('min:', vmin, ', max:', vmax)

    # colors = ['maroon', 'red', 'white', 'blue', 'midnightblue']
    # nodes = [0.0, 0.4, 0.5, 0.6, 1.0]
    # cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

    cmap = 'seismic_r'
    Fig = plt.figure(figsize=(16, 3.55))
    plt.rcParams.update({'font.size': 12.5, 'font.family': 'Arial'})

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

    ax = Fig.add_subplot(166)
    # cax = ax.matshow(sumwidth, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    cax = ax.matshow(ShorelineChange, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width (m)', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)
    # plt.title('Total Landscape Width')
    plt.title('Ocean Shoreline Change')

    Fig.suptitle(title, fontsize=18)
    plt.tight_layout()

    # plt.savefig("/Users/ianreeves/Desktop/Figure2_Cbar.svg")
    # plt.savefig("/Users/ianreeves/Desktop/Figure2_Cbar.png")

# --------------------------
if plot_param_space_pond:
    # Marsh + Pond
    Fig = plt.figure(figsize=(16, 3.55))
    # plt.rcParams.update({'font.size': 10, 'font.family': 'Arial'})

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

    Fig.suptitle(title, fontsize=18)
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
    Fig = plt.figure(figsize=(10, 5))
    plt.rcParams.update({'font.size': 14, 'font.family': 'Arial'})
    ax = Fig.add_subplot(122)
    plt.xlabel("RSLR (mm/yr)")
    plt.ylabel("External SSC (mg/L)")
    plt.title("Bay Width Change via Mainland Edge")
    cax = ax.matshow((MLMarshWidth + MLPondWidth + ForestWidth) * -1, origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)

    # Bay Width Change From Back-Barrier Marsh
    ax = Fig.add_subplot(121)
    plt.xlabel("RSLR (mm/yr)")
    plt.ylabel("External SSC (mg/L)")
    plt.title("Bay Width Change via Back-Barrier Edge")
    cax = ax.matshow((BayWidth + MLMarshWidth + MLPondWidth + ForestWidth), origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    ax.xaxis.set_ticks_position('bottom')
    # cbar = Fig.colorbar(cax)
    # cbar.set_label('Change in Width [m]', rotation=270, labelpad=20)
    ax.set_xticklabels(xtic)
    ax.set_yticklabels(ytic)
    plt.tight_layout()

    # plt.savefig("/Users/ianreeves/Desktop/Figure4.svg")
    # plt.savefig("/Users/ianreeves/Desktop/Figure4.png")

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

