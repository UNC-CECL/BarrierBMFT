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


Sim = 1

# ==================================================================================================================================================================================
# Specify data

# ~ File bank ~
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_28/'  # 49
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_29/'  # 50
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_30/'  # 51
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_31/'  # 52
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_32/'  # 53
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_33/'  # 54
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_18_15/'  # 55
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_27/'  # 56
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_29/'  # 57
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_30/'  # 58
# '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_31/'  # 59

data_files = [
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_30/',  # 51
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_31/',  # 52
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_27/',  # 56  <-- drowning
    # '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_16_29/',  # 57

    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_00_28/',  # 49
    '/Volumes/LACIE SHARE/Reeves/BarrierBMFT/Data-Results/Batch_2022_0425_18_15/',  # 55
]

file_num = len(data_files)

# Initialize data arrays
BBMarshWidth = []
MLMarshWidth = []
BBPondWidth = []
MLPondWidth = []
DuneHeight = []


# ==================================================================================================================================================================================
# Load data

for n in range(file_num):
    filename = data_files[n]
    widthsTS = np.load(filename + 'Sim' + str(Sim) + '_widthsTS.npy', allow_pickle=True)
    statsTS = np.load(filename + 'Sim' + str(Sim) + '_statsTS.npy', allow_pickle=True)

    BBmarsh_w = widthsTS[:, 1]
    MLmarsh_w = widthsTS[:, 3]
    BBpond_w = widthsTS[:, 5]
    MLpond_w = widthsTS[:, 6]
    aHd = statsTS[:, 4]  # [m] Average dune height

    BBMarshWidth = np.append(BBMarshWidth, BBmarsh_w)
    MLMarshWidth = np.append(MLMarshWidth, MLmarsh_w)
    BBPondWidth = np.append(BBPondWidth, BBpond_w)
    MLPondWidth = np.append(MLPondWidth, MLpond_w)
    DuneHeight = np.append(DuneHeight, aHd)


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
MarshBinWidth = 25
MarshBinNum = int(((MarshBinStop - MarshBinStart) / MarshBinWidth) + 2)
MarshBin = np.linspace(MarshBinStart, MarshBinStop, MarshBinNum)

Fig = plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size':12})

ax = Fig.add_subplot(1, 2, 1)
ax.hist(DuneHeight, bins=DuneBin, density=True)
ax.xaxis.set_ticks_position('bottom')
plt.ylabel('PDF', labelpad=15)
plt.xlabel('Dune Height [m]', labelpad=15)

ax = Fig.add_subplot(1, 2, 2)
ax.hist(BBMarshWidth, bins=MarshBin, density=True)
ax.xaxis.set_ticks_position('bottom')
plt.ylabel('PDF', labelpad=15)
plt.xlabel('Back-Barrier Marsh Width [m]', labelpad=15)

plt.show()





























