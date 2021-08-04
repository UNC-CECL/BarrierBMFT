"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 4 August 2021
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools import plot as B3Dfunc

# ==================================================================================================================================================================================
# Create an instance of the BMI class
barrierbmft = BarrierBMFT()


# # specify data directory with initial conditions
# datadir = "Barrier3D/"
# input_file = "barrier3d-parameters.yaml"
# barrierbmft.barrier3d.initialize(datadir + input_file)


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
plt.plot(barrierbmft.bmftc.elevation[barrierbmft.bmftc.endyear - 1, :])
plt.xlabel("Distance")
plt.ylabel("Elevation [m MSL]")

# ===========
plt.figure()
plt.plot(barrierbmft.bmftc.fetch[barrierbmft.bmftc.startyear: barrierbmft.bmftc.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Fetch [m]")

# ===========
plt.figure()
fig = plt.gcf()
fig.set_size_inches(14, 18)
plt.rcParams.update({"font.size": 17})

# Interior Width
plt.subplot(5, 1, 1)
plt.plot(barrierbmft.bmftc.fetch[barrierbmft.bmftc.startyear: barrierbmft.bmftc.endyear])
plt.xlabel("Time [yr]")
plt.ylabel("Bay Fetch [m]")

# Interior Width
plt.subplot(5, 1, 2)
aiw = [a * 10 for a in barrierbmft.barrier3d._model._InteriorWidth_AvgTS]
plt.plot(aiw)
plt.ylabel("Avg. Width (m)")  # Avergae Interior Width

# Shoreline Change
scts = [(x - barrierbmft.barrier3d._model._x_s_TS[0]) * 10 for x in barrierbmft.barrier3d._model._x_s_TS]
plt.subplot(5, 1, 3)
plt.plot(scts)
plt.ylabel("Shoreline Position (m)")

# Overwash Flux
plt.subplot(5, 1, 4)
plt.plot(barrierbmft.barrier3d._model._QowTS)
plt.ylabel("Qow (m^3/m)")

# Shoreface Flux
plt.subplot(5, 1, 5)
plt.plot(barrierbmft.barrier3d._model._QsfTS)
plt.ylabel("Qsf (m^3/m)")

# ===========
plt.show()
