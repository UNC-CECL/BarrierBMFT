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


Sim = 5
title = 'Low Slope'


# ==================================================================================================================================================================================
# Specify data

data_files_a = [
    # Group 4: Bistable, shallow
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_28/',  # 49
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_18_15/',  # 55  <-- 1 drowned

    # Group 5: Bistable, moderate
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_30/',  # 51
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_31/',  # 52
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_29/',  # 57

    # Group 6: Bistable, steep
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_32/',  # 53
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_30/',  # 58
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_31/',  # 59
]

data_files_b = [
    # Group 7: Low, shallow
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_21/',  # 61
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_00_08/',  # 69  <-- running
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0427_12_37/',  # 70  <-- running

    # Group 8: Low, moderate
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_22/',  # 62
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_23/',  # 63
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_01/',  # 66

    # Group 9: Low, steep
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_24/',  # 64
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0426_14_25/',  # 65
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
DuneBinStart = 0
DuneBinStop = 1.6
DuneBinWidth = 0.1
DuneBinNum = int(((DuneBinStop - DuneBinStart) / DuneBinWidth) + 2)
DuneBin = np.linspace(DuneBinStart, DuneBinStop, DuneBinNum)

MarshBinStart = 0
MarshBinStop = 900
MarshBinWidth = 50
MarshBinNum = int(((MarshBinStop - MarshBinStart) / MarshBinWidth) + 2)
MarshBin = np.linspace(MarshBinStart, MarshBinStop, MarshBinNum)

Fig = plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size':12})

ax = Fig.add_subplot(1, 2, 1)
ax.hist(DuneHeight_a, bins=DuneBin, density=True, color='red', alpha=0.5)
ax.hist(DuneHeight_b, bins=DuneBin, density=True, color='black', alpha=0.5)
ax.xaxis.set_ticks_position('bottom')
plt.ylabel('PDF', labelpad=15)
plt.xlabel('Dune Height [m]', labelpad=15)

ax = Fig.add_subplot(1, 2, 2)
ax.hist(BBMarshWidth_a, bins=MarshBin, density=True, color='red', alpha=0.5)
ax.hist(BBMarshWidth_b, bins=MarshBin, density=True, color='black', alpha=0.5)
ax.xaxis.set_ticks_position('bottom')
plt.legend(['Bistable (r=0.75)', 'Low (r=0.45)'])
plt.ylabel('PDF', labelpad=15)
plt.xlabel('Back-Barrier Marsh Width [m]', labelpad=15)

Fig.suptitle(title, fontsize=18)
plt.tight_layout()
plt.show()





























