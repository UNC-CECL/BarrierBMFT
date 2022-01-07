"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 17 August 2021
"""

import time
import math
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools import plot as B3Dfunc

# ==================================================================================================================================================================================
# Create an instance of the BMI class
barrierbmft = BarrierBMFT()


# ==================================================================================================================================================================================
# Run the BarrierBMFT model

# Record start time
Time = time.time()

# Loop through time
for time_step in range(int(barrierbmft.bmftc.dur)):

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

    # Run time step
    barrierbmft.update(time_step)

    # Check for breaks
    if barrierbmft.BMFTC_Break or barrierbmft.Barrier3D_Break:
        break

# Print elapsed time of simulation
print()
SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")

# Show figure(s)
plt.show()


# ==================================================================================================================================================================================
# Plot

# ===========
# plt.figure()
# plt.plot(barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.endyear - 1, :])
# plt.xlabel("Distance")
# plt.ylabel("Back-Barrier Elevation [m MSL]")
#
# plt.figure()
# plt.plot(barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.endyear - 1, :])
# plt.xlabel("Distance")
# plt.ylabel("Mainland Elevation [m MSL]")

# ===========
plt.figure()
plt.plot(barrierbmft.bmftc_ML.fetch[barrierbmft.bmftc_ML.startyear: barrierbmft.bmftc_ML.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Fetch [m]")
plt.rcParams.update({"font.size": 11})

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(7, 16)
plt.subplot(3, 1, 1)
marsh_width_TS = barrierbmft.bmftc_BB.Forest_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear] - barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear]
plt.plot(marsh_width_TS)
plt.xlabel("Time [yr]")
plt.ylabel("Back-Barrier Marsh Width [m]")
plt.subplot(3, 1, 2)
plt.plot(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("Back-Barrier Marsh Edge Location [m]")
plt.subplot(3, 1, 3)
plt.plot(barrierbmft.bmftc_BB.Forest_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("PyBMFT-C Back-Barrier Forest Edge Location [m]")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(7, 16)
plt.subplot(3, 1, 1)
bbscts = [(x - barrierbmft.barrier3d.model.x_b_TS[0]) * 10 for x in barrierbmft.barrier3d.model.x_b_TS][1:]
plt.plot(bbscts)
plt.ylabel("Barrier3D Forest Edge [m]")
plt.subplot(3, 1, 2)
bbbmftc = (barrierbmft.bmftc_BB.Forest_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear] - barrierbmft.bmftc_BB.Forest_edge[barrierbmft.bmftc_BB.startyear]) * (-1)
plt.plot(bbbmftc)
plt.ylabel("PyBMFT-C Forest Edge [m]")
plt.subplot(3, 1, 3)
plt.plot(bbscts - bbbmftc)
plt.xlabel("Time [yr]")
plt.ylabel("Difference [m]")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(14, 18)
# plt.rcParams.update({"font.size": 12})

# Interior Width
plt.subplot(6, 1, 1)
plt.plot(barrierbmft.bmftc_ML.fetch[barrierbmft.bmftc_ML.startyear: barrierbmft.bmftc_ML.endyear])
plt.ylabel("Bay Fetch [m]")

# Back-Barrier Shoreline Change
plt.subplot(6, 1, 2)
bbscts = [(x - barrierbmft.barrier3d.model.x_b_TS[0]) * 10 for x in barrierbmft.barrier3d.model.x_b_TS]
plt.plot(bbscts)
plt.ylabel("BB Shoreline Change (m)")  # Avergae Interior Width

# Interior Width
plt.subplot(6, 1, 3)
aiw = [a * 10 for a in barrierbmft.barrier3d.model.InteriorWidth_AvgTS]
plt.plot(aiw)
plt.ylabel("Avg. Width (m)")  # Avergae Interior Width

# Shoreline Change
scts = [(x - barrierbmft.barrier3d.model.x_s_TS[0]) * 10 for x in barrierbmft.barrier3d.model.x_s_TS]
plt.subplot(6, 1, 4)
plt.plot(scts)
plt.ylabel("Shoreline Position (m)")

# Overwash Flux
plt.subplot(6, 1, 5)
plt.plot(barrierbmft.barrier3d.model.QowTS)
plt.ylabel("Qow (m^3/m)")

# Shoreface Flux
plt.subplot(6, 1, 6)
plt.plot(barrierbmft.barrier3d.model.QsfTS)
plt.ylabel("Qsf (m^3/m)")
plt.xlabel("Time [yr]")

# ===========
# BB_transect = np.flip(barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.dur, int(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_BB.dur]):])
# if barrierbmft.x_b_TS_ML[barrierbmft.bmftc_ML.dur] < 0:
#     ML_transect = np.append(np.ones([abs(int(barrierbmft.x_b_TS_ML[barrierbmft.bmftc_ML.dur]))]) * barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.dur, 1], barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.dur, :])
# elif barrierbmft.x_b_TS_ML[barrierbmft.bmftc_ML.dur] > 0:
#     ML_transect = barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.dur, int(barrierbmft.x_b_TS_ML[barrierbmft.bmftc_ML.dur]):]
# whole_transect = np.append(BB_transect, ML_transect)
# plt.figure()
# plt.plot(whole_transect)
# plt.xlabel("Distance")
# plt.ylabel("Elevation [m MSL]")

# ===========
# plt.figure()
# plt.plot((barrierbmft.barrier3d.model.x_b_TS - barrierbmft.barrier3d.model.x_b_TS[0]) * 10)
# plt.xlabel("Time [yr]")
# plt.ylabel("Barrier3D Back-Barrier Shoreline Position [m]")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(7, 15)
plt.subplot(2, 1, 1)
plt.plot(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_BB.startyear: barrierbmft.bmftc_BB.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("BMFTC BB Marsh Edge [m]")
plt.subplot(2, 1, 2)
plt.plot(barrierbmft.bmftc_ML.Marsh_edge[barrierbmft.bmftc_ML.startyear: barrierbmft.bmftc_ML.endyear + 1])
plt.xlabel("Time [yr]")
plt.ylabel("BMFTC ML Marsh Edge [m]")

# ===========
# plt.figure()
# plt.plot(barrierbmft.cumul_len_change)
# plt.xlabel("Time [yr]")
# plt.ylabel("Cumulative Length Gain/Loss From Interpolation Rounding [m]")

# ===========
# plt.figure()
# plt.plot(barrierbmft.delta_fetch_BB_TS, c="blue")
# plt.plot(barrierbmft.delta_fetch_ML_TS, c="red")
# plt.xlabel("Time [yr]")
# plt.ylabel("Change in Marsh Edge Location from Bay Processes [m]")
# plt.legend(["Back-Barrier", "Mainland"])

# ===========
# Barrier Elevation (end)
# B3Dfunc.plot_ElevTMAX(
#     barrierbmft.bmftc.dur,
#     barrierbmft.barrier3d.model._DuneDomain,
#     barrierbmft.barrier3d.model._DomainTS,
#     barrierbmft.barrier3d.model._BermEl,
#     barrierbmft.barrier3d.model._Shrub_ON,
#     barrierbmft.barrier3d.model._PercentCoverTS,
#     barrierbmft.barrier3d.model._DeadPercentCoverTS,
#     barrierbmft.barrier3d.model._DuneWidth,
# )


# ===========
# Width of coastal landscape

widths = barrierbmft.LandscapeTypeWidth_TS
total_width = np.sum(widths, axis=1)
barrier = widths[:, 0] / total_width
BBmarsh = widths[:, 1] / total_width
bay = widths[:, 2] / total_width
MLmarsh = widths[:, 3] / total_width
forest = widths[:, 4] / total_width
barrier_marsh = barrier + BBmarsh

plt.figure(figsize=(15, 7))
plt.plot(barrier, c="black")
plt.plot(BBmarsh, c="red")
plt.plot(bay, c="blue")
plt.plot(MLmarsh, c="yellow")
plt.plot(forest, c="green")
plt.plot(barrier_marsh, c="purple")
plt.xlabel("Time [yr]")
plt.ylabel("Proportion of Entire Landscape Width")
plt.legend(["Barrier", "Back-Barrier Marsh", "Bay", "Mainland Marsh", "Forest", "Barrier + Back-Barrier Marsh"])



# ===========
# Barrier Animation
# BeachWidth = 6
# OriginY = 10
# AniDomainWidth = int(
#     max(barrierbmft.barrier3d.model.InteriorWidth_AvgTS) + BeachWidth + abs(barrierbmft.barrier3d.model.ShorelineChange) + OriginY + 15 + (marsh_width_TS[-1] / 10))  # was +15
#
# for t in range(barrierbmft.barrier3d.model.TMAX):
#     # Build beach elevation domain
#     BeachDomain = np.zeros([BeachWidth, barrierbmft.barrier3d.model.BarrierLength])
#     berm = math.ceil(BeachWidth * 0.65)
#     BeachDomain[berm: BeachWidth + 1, :] = barrierbmft.barrier3d.model.BermEl
#     add = (barrierbmft.barrier3d.model.BermEl - barrierbmft.barrier3d.model.SL) / berm
#     for i in range(berm):
#         BeachDomain[i, :] = barrierbmft.barrier3d.model.SL + add * i
#
#     # Make animation frame domain
#     Domain = barrierbmft.barrier3d.model.DomainTS[t] * 10
#     Dunes = (barrierbmft.barrier3d.model.DuneDomain[t, :, :] + barrierbmft.barrier3d.model.BermEl) * 10
#     Dunes = np.rot90(Dunes)
#     Dunes = np.flipud(Dunes)
#     Beach = BeachDomain * 10
#     Domain = np.vstack([Beach, Dunes, Domain])
#     Domain[Domain < -1] = -1
#     AnimateDomain = np.ones([AniDomainWidth + 1, barrierbmft.barrier3d.model.BarrierLength]) * -1
#     widthTS = len(Domain)
#     scts = [(x - barrierbmft.barrier3d.model.x_s_TS[0]) for x in barrierbmft.barrier3d.model.x_s_TS]
#     if scts[t] >= 0:
#         OriginTstart = OriginY + math.floor(scts[t])
#     else:
#         OriginTstart = OriginY + math.ceil(scts[t])
#     OriginTstop = OriginTstart + widthTS
#     AnimateDomain[OriginTstart:OriginTstop, 0: barrierbmft.barrier3d.model.BarrierLength] = Domain
#
#     # Plot and save
#     elevFig1 = plt.figure(figsize=(7, 12))
#     ax = elevFig1.add_subplot(111)
#     cax = ax.matshow(AnimateDomain, origin="lower", cmap="terrain", vmin=-2, vmax=4.0)  # , interpolation='gaussian') # analysis:ignore
#
#     ax.xaxis.set_ticks_position("bottom")
#     # cbar = elevFig1.colorbar(cax)
#     plt.xlabel("Alongshore Distance (dam)")
#     plt.ylabel("Cross-Shore Diatance (dam)")
#     plt.title("Interior Elevation")
#     plt.tight_layout()
#     timestr = "Time = " + str(t) + " yrs"
#     newpath = "Output/SimFrames/"
#     if not os.path.exists(newpath):
#         os.makedirs(newpath)
#     plt.text(1, 1, timestr)
#     name = "Output/SimFrames/elev_" + str(t)
#     elevFig1.savefig(name)  # dpi=200
#     plt.close(elevFig1)
#
# frames = []
# for filenum in range(1, barrierbmft.barrier3d.model.TMAX):
#     filename = "Output/SimFrames/elev_" + str(filenum) + ".png"
#     frames.append(imageio.imread(filename))
# imageio.mimsave("Output/SimFrames/elev.gif", frames, "GIF-FI")
# print()
# print("[ * GIF successfully generated * ]")

# ===========
# Transect Animation

for t in range(int(barrierbmft.bmftc.dur)):

    # Combine transects
    BB_transect = np.flip(barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.startyear + t - 1, int(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_ML.startyear + t]):])
    if barrierbmft.x_b_TS_ML[t] < 0:
        ML_transect = np.append(np.ones([abs(int(barrierbmft.x_b_TS_ML[t]))]) * barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, 1], barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, :])
    elif barrierbmft.x_b_TS_ML[t] > 0:
        ML_transect = barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, int(barrierbmft.x_b_TS_ML[t]):]

    whole_transect = np.append(BB_transect, ML_transect)

    # Extract marsh & forest locations to plot
    x_forest = [barrierbmft.bmftc_BB.B - int(barrierbmft.bmftc_BB.Forest_edge[barrierbmft.bmftc_BB.startyear + t]), len(BB_transect) - int(barrierbmft.x_b_TS_ML[t]) + int(barrierbmft.bmftc_ML.Forest_edge[barrierbmft.bmftc_ML.startyear + t])]
    x_marsh = [len(BB_transect) - 1, len(BB_transect) - int(barrierbmft.x_b_TS_ML[t]) + int(barrierbmft.bmftc_ML.Marsh_edge[barrierbmft.bmftc_ML.startyear + t])]
    y_forest = whole_transect[x_forest]
    y_marsh = whole_transect[x_marsh]

    # Plot and save
    transectFig = plt.figure(figsize=(15, 7))
    plt.plot(whole_transect, c="black")
    plt.scatter(x_forest, y_forest, c="green")
    plt.scatter(x_marsh, y_marsh, c="brown")
    plt.ylim(-1.5, 4)
    plt.xlabel("Cross-shore Distance (m)")
    plt.ylabel("Elevation (m)")
    plt.tight_layout()
    timestr = "Time = " + str(t) + " yrs"
    newpath = "Output/SimFrames/"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    plt.text(0, 0, timestr)
    name = "Output/SimFrames/transect_" + str(t)
    transectFig.savefig(name)  # dpi=200
    plt.close(transectFig)

frames = []
for filenum in range(int(barrierbmft.bmftc.dur)):
    filename = "Output/SimFrames/transect_" + str(filenum) + ".png"
    frames.append(imageio.imread(filename))
imageio.mimsave("Output/SimFrames/transect.gif", frames, "GIF-FI")
print()
print("[ * GIF successfully generated * ]")

# ===========
plt.show()
