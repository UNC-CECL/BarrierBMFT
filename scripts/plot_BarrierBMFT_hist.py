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


Sim = 5  # Parameter space location
Sim_plot = 55  # Simulation number for plotting  # mod:55-8, shallow:55-5, steep:55-7
Data_Sim_plot = 8  # Data file number for plotting
title = 'Moderate Slope'
leg1 = 'Fast dune growth (r=0.75)'  # Legend
leg2 = 'Slow dune growth (r=0.45)'  # Legend
leg3 = 'VCR Observations'

# ==================================================================================================================================================================================
# Specify data

data_files_a = [
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
]

data_files_b = [
    # Group 7: Slow, shallow
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_21/',  # 61
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_00_08/',  # 69
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_37/',  # 70
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0429_18_27/',  # 92

    # Group 8: Slow, moderate
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_22/',  # 62
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_23/',  # 63
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_01/',  # 66
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0430_16_14/',  # 100
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0511_18_35/',  # 118
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_16/',  # 121
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0512_01_17/',  # 122

    # Group 9: Slow, steep
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_24/',  # 64
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_25/',  # 65
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0428_18_31/',  # 76
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0501_16_37/',  # 101
]

file_num_a = len(data_files_a)
file_num_b = len(data_files_b)

# Initialize data arrays
BBMarshWidth_a = []
MLMarshWidth_a = []
BBPondWidth_a = []
MLPondWidth_a = []
DuneHeight_a = []

BBMarshWidth_b = []
MLMarshWidth_b = []
BBPondWidth_b = []
MLPondWidth_b = []
DuneHeight_b = []


# ==================================================================================================================================================================================
# Load data

for n in range(file_num_a):
    filename = data_files_a[n]
    widthsTS = np.load(filename + 'Sim' + str(Sim) + '_widthsTS.npy', allow_pickle=True)
    statsTS = np.load(filename + 'Sim' + str(Sim) + '_statsTS.npy', allow_pickle=True)

    BBmarsh_w = widthsTS[:, 1]
    MLmarsh_w = widthsTS[:, 3]
    BBpond_w = widthsTS[:, 5]
    MLpond_w = widthsTS[:, 6]
    aHd = statsTS[:, 4]  # [m] Average dune height
    x_f_BB = statsTS[:, 1]  # [m] Average dune height

    BBMarshWidth_a = np.append(BBMarshWidth_a, BBmarsh_w)
    MLMarshWidth_a = np.append(MLMarshWidth_a, MLmarsh_w)
    BBPondWidth_a = np.append(BBPondWidth_a, BBpond_w)
    MLPondWidth_a = np.append(MLPondWidth_a, MLpond_w)
    DuneHeight_a = np.append(DuneHeight_a, aHd)

for n in range(file_num_b):
    filename = data_files_b[n]
    widthsTS = np.load(filename + 'Sim' + str(Sim) + '_widthsTS.npy', allow_pickle=True)
    statsTS = np.load(filename + 'Sim' + str(Sim) + '_statsTS.npy', allow_pickle=True)

    BBmarsh_w = widthsTS[:, 1]
    MLmarsh_w = widthsTS[:, 3]
    BBpond_w = widthsTS[:, 5]
    MLpond_w = widthsTS[:, 6]
    aHd = statsTS[:, 4]  # [m] Average dune height

    BBMarshWidth_b = np.append(BBMarshWidth_b, BBmarsh_w)
    MLMarshWidth_b = np.append(MLMarshWidth_b, MLmarsh_w)
    BBPondWidth_b = np.append(BBPondWidth_b, BBpond_w)
    MLPondWidth_b = np.append(MLPondWidth_b, MLpond_w)
    DuneHeight_b = np.append(DuneHeight_b, aHd)


# ==================================================================================================================================================================================
# Plot

# Dune and BB width histograms
DuneBinStart = 0.1
DuneBinStop = 1.6
DuneBinWidth = 0.1
DuneBinNum = int(((DuneBinStop - DuneBinStart) / DuneBinWidth) + 1)
DuneBin = np.linspace(DuneBinStart, DuneBinStop, DuneBinNum)

MarshBinStart = 0
MarshBinStop = 1000
MarshBinWidth = 50
MarshBinNum = int(((MarshBinStop - MarshBinStart) / MarshBinWidth) + 1)
MarshBin = np.linspace(MarshBinStart, MarshBinStop, MarshBinNum)

Fig = plt.figure(figsize=(15, 12))
plt.rcParams.update({'font.size':16})

ax = Fig.add_subplot(2, 2, 3)
ax.hist(DuneHeight_a, bins=DuneBin, density=True, color='red', alpha=0.5)
ax.hist(DuneHeight_b, bins=DuneBin, density=True, color='black', alpha=0.5)
ax.xaxis.set_ticks_position('bottom')
plt.xlim([0, DuneBinStop])
plt.legend([leg1, leg2])
plt.ylabel('PDF (Model)', labelpad=15)
plt.xlabel('Dune Height [m]', labelpad=15)

ax = Fig.add_subplot(2, 2, 4)
ax.hist(BBMarshWidth_a, bins=MarshBin, density=True, color='red', alpha=0.5)
ax.hist(BBMarshWidth_b, bins=MarshBin, density=True, color='black', alpha=0.5)
ax.xaxis.set_ticks_position('bottom')
plt.ylabel('PDF (Model)', labelpad=15)
plt.xlabel('Back-Barrier Marsh Width [m]', labelpad=15)

# ADD OBSERVATIONAL DATA - Walters et al. (2014)
BBWidthObs = np.loadtxt('/Users/ianreeves/Desktop/WaltersBBMarshWidthObs-2010-2.csv', delimiter=',')
BBWidthObs = BBWidthObs[:, 0]

ax2 = ax.twinx()
# ax2.set_ylim([0, 0.005])
ax2.hist(BBWidthObs, bins=MarshBin, density=True, facecolor='none', edgecolor='black', linewidth=2, alpha=0.8)
ax2.set_ylabel("PDF (Observations)", labelpad=15)
plt.legend([leg3])


# ------
SimWidths = np.load(data_files_a[Data_Sim_plot - 1] + 'Sim' + str(Sim_plot) + '_widthsTS.npy', allow_pickle=True)
SimStats = np.load(data_files_a[Data_Sim_plot - 1] + 'Sim' + str(Sim_plot) + '_statsTS.npy', allow_pickle=True)

aHd = SimStats[:, 4]
x_f_BB = SimStats[:, 1]
x_f_BB = [(x - x_f_BB[0])for x in x_f_BB]
x_s = SimStats[:, 7]
x_s = [(x - x_s[0]) for x in x_s]

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
ax1 = Fig.add_subplot(2, 2, 1)

ax1.set_xlabel('Year')
ax1.set_ylabel('BB Marsh Width [m]', color='red')
ax1.plot(BBmarsh, c="red", linewidth=2)
ax1.tick_params(axis='y', labelcolor='red')

ax2 = ax1.twinx()
ax2.set_ylabel("Dune Height [m]", color='black')
ax2.plot(aHd, c="black", linewidth=2)
ax2.tick_params(axis='y', labelcolor='black')


# ------
# Back-barrier Shoreline Change
# ax1 = Fig.add_subplot(2, 2, 4)
#
# ax1.set_xlabel('Year')
# ax1.set_ylabel('BB Marsh Width [m]', color='red')
# ax1.plot(BBmarsh, c="red")
# ax1.tick_params(axis='y', labelcolor='red')
#
# ax2 = ax1.twinx()
# ax2.set_ylabel("Back-Barrier Shoreline Position [m]", color='black')
# ax2.plot(x_f_BB, c="black")
# ax2.tick_params(axis='y', labelcolor='black')

# Ocean Shoreline Change
ax1 = Fig.add_subplot(2, 2, 2)

ax1.set_xlabel('Year')
ax1.set_ylabel('BB Marsh Width [m]', color='red')
ax1.plot(BBmarsh, c="red", linewidth=2)
ax1.tick_params(axis='y', labelcolor='red')

ax2 = ax1.twinx()
ax2.set_ylabel("Ocean Shoreline Position [m]", color='black')
ax2.plot(x_s, c="black", linewidth=2)
ax2.tick_params(axis='y', labelcolor='black')


# --------------------------
Fig.suptitle(title, fontsize=18)
plt.tight_layout()

plt.savefig("/Users/ianreeves/Desktop/Figure3.svg")
plt.show()
























