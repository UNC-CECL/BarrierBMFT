"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 17 August 2021
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools import plot as B3Dfunc

# ==================================================================================================================================================================================
# Create an instance of the BMI class
barrierbmft = BarrierBMFT(
    coupling_on=True
)


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
plt.figure()
plt.plot(barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.endyear - 1, :])
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")


# ===========
plt.figure()
plt.plot(barrierbmft.bmftc_ML.fetch[barrierbmft.bmftc_ML.startyear: barrierbmft.bmftc_ML.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Fetch [m]")

# ===========
plt.figure()

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
plt.ylabel("Back-Barrier Forest Edge Location [m]")


# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(14, 18)
plt.rcParams.update({"font.size": 12})

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
barrier_transect = np.mean(barrierbmft.barrier3d.model.InteriorDomain, axis=1) * 10
x = np.linspace(1, len(barrier_transect) * 10, num=len(barrier_transect) * 10)
xp = np.linspace(1, len(barrier_transect), num=len(barrier_transect)) * 10
barrier_transect = np.interp(x, xp, barrier_transect)  # Interpolate from dam to m
subaqueous = np.where(barrier_transect <= 0)[0]
bmf_transect = barrierbmft.bmftc.elevation[barrierbmft.bmftc.endyear - 1, :] - barrierbmft.bmftc.msl[-1] - barrierbmft.bmftc.amp
bmf_transect = bmf_transect[len(subaqueous):]
whole_transect = np.append(barrier_transect, bmf_transect)
plt.figure()
plt.plot(whole_transect)
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")


# ===========
plt.figure()
plt.plot((barrierbmft.barrier3d.model.x_b_TS - barrierbmft.barrier3d.model.x_b_TS[0]) * 10)
plt.xlabel("Time [yr]")
plt.ylabel("Back-Barrier Shoreline Position [m]")


# ===========
# Barrier Elevation (end)
B3Dfunc.plot_ElevTMAX(
    barrierbmft.bmftc.dur,
    barrierbmft.barrier3d.model._DuneDomain,
    barrierbmft.barrier3d.model._DomainTS,
    barrierbmft.barrier3d.model._BermEl,
    barrierbmft.barrier3d.model._Shrub_ON,
    barrierbmft.barrier3d.model._PercentCoverTS,
    barrierbmft.barrier3d.model._DeadPercentCoverTS,
    barrierbmft.barrier3d.model._DuneWidth,
)


# ===========
plt.show()

