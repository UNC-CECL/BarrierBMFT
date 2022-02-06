"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 20 January 2022
"""

import time
import os
import numpy as np
from datetime import datetime

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools.input_files import yearly_storms

# ==================================================================================================================================================================================
# Define batch parameters

Num = 2  # Number of runs at each combinations of parameter values
SimDur = 3  # [Yr] Duration of each simulation

# Parameter values
rslr = [3, 4, 9, 12, 15]
co = [20, 30, 40, 50, 60]
slope = [0.003]

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
ShorelineChange = np.zeros([Num, len(rslr), len(co), len(slope)])


# ==================================================================================================================================================================================
# Run the BarrierBMFT batch

# Record start time
Time = time.time()
Sim = 0

for n in range(Num):

    StormSeries = yearly_storms(
        datadir='Input/Barrier3D',
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        model_years=10,
        bPlot=False,
    )

    for r in range(len(rslr)):
        for c in range(len(co)):
            for s in range(len(slope)):

                Sim += 1

                # Create an instance of the BMI class and set input parameter values
                barrierbmft = BarrierBMFT(
                    time_step_count=SimDur,
                    relative_sea_level_rise=rslr[r],
                    reference_concentration=co[c],
                    slope_upland=slope[s],
                )

                # Update storm series
                barrierbmft.barrier3d.StormSeries = StormSeries

                # Run simulation: Loop through time
                for time_step in range(int(barrierbmft.bmftc.dur)):

                    # Print time step to screen
                    print("\r", "Time Step: ", time_step, "     Sim", Sim, "/", SimNum, "     (n", n, ", r", r, ", c", c, ")", end="")

                    # Run time step
                    barrierbmft.update(time_step)

                    # Check for breaks
                    if barrierbmft.BMFTC_Break or barrierbmft.Barrier3D_Break:
                        break

                # Calculate change in widths
                widths = barrierbmft.LandscapeTypeWidth_TS
                barrier_w = widths[-1, 0] - widths[0, 0]
                BBmarsh_w = widths[-1, 1] - widths[0, 1]
                bay_w = widths[-1, 2] - widths[0, 2]
                MLmarsh_w = widths[-1, 3] - widths[0, 3]
                forest_w = widths[-1, 4] - widths[0, 4]

                # Calculate shoreline change
                sc = int(barrierbmft.barrier3d.model.ShorelineChange)

                # Store in arrays
                BarrierWidth[n, r, c, s] = barrier_w
                BBMarshWidth[n, r, c, s] = BBmarsh_w
                BayWidth[n, r, c, s] = bay_w
                MLMarshWidth[n, r, c, s] = MLmarsh_w
                ForestWidth[n, r, c, s] = forest_w
                ShorelineChange[n, r, c, s] = sc

                # Save elevation array
                whole_transect = []

                for t in range(int(barrierbmft.bmftc.dur)):
                    BB_transect = np.flip(barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.startyear + t - 1, int(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_ML.startyear + t]):])
                    if barrierbmft.x_b_TS_ML[t] < 0:
                        ML_transect = np.append(np.ones([abs(int(barrierbmft.x_b_TS_ML[t]))]) * barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, 1], barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, :])
                    elif barrierbmft.x_b_TS_ML[t] > 0:
                        ML_transect = barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, int(barrierbmft.x_b_TS_ML[t]):]

                    whole_transect.append(np.append(BB_transect, ML_transect))  # Combine transects

                np.save(directory + '/Sim' + str(Sim) + '_elevation.npy', whole_transect)

# Save batch data arrays
np.save(directory + '/Widths_Barrier.npy', BarrierWidth)
np.save(directory + '/Widths_BBMarsh.npy', BBMarshWidth)
np.save(directory + '/Widths_Bay.npy', BayWidth)
np.save(directory + '/Widths_MLMarsh.npy', MLMarshWidth)
np.save(directory + '/Widths_Forest.npy', ForestWidth)
np.save(directory + '/ShorelineChange.npy', ShorelineChange)


# Print elapsed time of batch run
print()
BatchDuration = time.time() - Time
print()
print("Elapsed Time: ", BatchDuration, "sec")





# ==================================================================================================================================================================================
# Plot









