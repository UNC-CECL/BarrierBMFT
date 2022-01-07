"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 6 January 2022
"""

import time
import math
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
import ternary
from datetime import datetime

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools import plot as B3Dfunc

# ==================================================================================================================================================================================
# Define batch parameters

Num = 1  # Number of runs at each combinations of parameter values
SimDur = 125  # [Yr] Duration of each simulation

# Parameter values
rslr = [2, 6, 10, 14, 18]
co = [10, 30, 50, 70, 90]
slope = [0.001, 0.003, 0.005, 0.007, 0.009]


SimNum = len(rslr) * len(co) * len(slope)

# ==================================================================================================================================================================================
# Make new directory and data arrays

name = datetime.today().strftime('%Y_%m%d_%H_%M')
directory = 'Output/Batch_' + name
print('Batch_' + name)
os.makedirs(directory)

BarrierWidth = np.zeros([Num, len(rslr), len(co), len(slope)])
BBMarshWidth = np.zeros([Num, len(rslr), len(co), len(slope)])
BayWidth = np.zeros([Num, len(rslr), len(co), len(slope)])
MLMarshWidth = np.zeros([Num, len(rslr), len(co), len(slope)])
ForestWidth = np.zeros([Num, len(rslr), len(co), len(slope)])


# ==================================================================================================================================================================================
# Run the BarrierBMFT batch

# Record start time
Time = time.time()

for r in range(len(rslr)):
    for c in range(len(co)):
        for s in range(len(slope)):

            # Create an instance of the BMI class and set input parameter values
            barrierbmft = BarrierBMFT(
                time_step_count=SimDur,
                relative_sea_level_rise=rslr[r],
                reference_concentration=co[c],
                slope_upland=slope[s],
            )

            # Run simulation: Loop through time
            for time_step in range(int(barrierbmft.bmftc.dur)):

                # Print time step to screen
                print("\r", "Time Step: ", time_step, end="")

                # Run time step
                barrierbmft.update(time_step)

                # Check for breaks
                if barrierbmft.BMFTC_Break or barrierbmft.Barrier3D_Break:
                    break

            # Calculate widths at end of simulation
            widths = barrierbmft.LandscapeTypeWidth_TS
            barrier_w = widths[-1, 0]
            BBmarsh_w = widths[-1, 1]
            bay_w = widths[-1, 2]
            MLmarsh_w = widths[-1, 3]
            forest_w = widths[-1, 4]

            # Store width in arrays
            BarrierWidth[0, r, c, s] = barrier_w
            BBMarshWidth[0, r, c, s] = BBmarsh_w
            BayWidth[0, r, c, s] = bay_w
            MLMarshWidth[0, r, c, s] = MLmarsh_w
            ForestWidth[0, r, c, s] = forest_w

# Save batch data arrays
np.save(directory + '/Widths_Barrier.npy', BarrierWidth)
np.save(directory + '/Widths_BBMarsh.npy', BBMarshWidth)
np.save(directory + '/Widths_Bay.npy', BayWidth)
np.save(directory + '/Widths_MLMarsh.npy', MLMarshWidth)
np.save(directory + '/Widths_Forest.npy', ForestWidth)


# Print elapsed time of batch run
print()
BatchDuration = time.time() - Time
print()
print("Elapsed Time: ", BatchDuration, "sec")





# ==================================================================================================================================================================================
# Plot

# # Barrier
# BarrierWidth = BarrierWidth[0, :, :, :]  # Note: in future need to average across all duplicates
#
# # Reformat data
# Barrier_Loc = []
# Barrier_Val = []
# for r in range(len(rslr)):
#     for c in range(len(co)):
#         for s in range(len(slope)):
#             Barrier_Loc.append((rslr[r], co[c], slope[s] * 10))
#             Barrier_Val.append(BarrierWidth[r, c, s])
#
#
# # Plot ternary diagram
# scale = 50
# figure, tax = ternary.figure(scale=scale)
# figure.set_size_inches(10, 10)
#
# # Set labels
# tax.set_title("Barrier", fontsize=20)
# tax.right_axis_label("Sediment Concentration")
# tax.left_axis_label("Upland Slope")
# tax.bottom_axis_label("RSLR")
#
# # Set custom axis limits
# tax.set_axis_limits({'b': [min(rslr), max(rslr)], 'l': [min(slope) * 10, max(slope) * 10], 'r': [min(co), max(co)]})
# tax.get_ticks_from_axis_limits()
# tax.set_custom_ticks(fontsize=10, offset=0.02)
#
# tax.scatter(Barrier_Loc, marker='s', color='red')
#
# tax.ax.set_aspect('equal', adjustable='box')
# tax._redraw_labels()
#
# tax.boundary(linewidth=2.0)
# tax.gridlines(multiple=5, color='blue')
# # tax.ticks(axis='lbr', linewidth=1, multiple=5)
# tax.clear_matplotlib_ticks()
# tax.get_axes().axis('off')
#
# tax.show()

# ####################
# filename = '/Users/ianreeves/PycharmProjects/BarrierBMFT/Output/Batch_2022_0106_23_42/Widths_Barrier.npy'
# BarrierWidth = np.load(filename)
# BarrierWidth = BarrierWidth[0, :, :, :]







