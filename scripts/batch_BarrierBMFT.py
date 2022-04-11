"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the PyBMFT-C model

Copyright Ian RB Reeves
Last updated: 14 March 2022
"""

import time
import os
import numpy as np
from datetime import datetime
import multiprocessing
import shutil
import random
from joblib import Parallel, delayed

from barrierbmft.barrierbmft import BarrierBMFT
from barrier3d.tools.input_files import yearly_storms

# ==================================================================================================================================================================================
# Define batch parameters

Num = 3  # Number of runs at each combinations of parameter values
SimDur = 400  # [Yr] Duration of each simulation

# Parameter values
rslr = [3, 6, 9, 12, 15]
co = [20, 30, 40, 50, 60]
slope = [0.005]

add = 0  # Pad input labels, optional


# ==================================================================================================================================================================================
# Batch code

# Makes storm series
def makeStormSeries(
        datadir,
        storm_list_filename,
        StormStart,
        TMAX,
        mean_storm,
        SD_storm,
        MHW,
        output_filename,
):

    # Load input file
    StormList = np.load(datadir + storm_list_filename)

    # Time series
    StormSeries = np.zeros([StormStart, 5])
    for t in range(StormStart, TMAX):
        # Calculate number of storms in year
        numstorm = round(np.random.normal(mean_storm, SD_storm))

        if numstorm < 0:
            numstorm = 0
        stormTS = np.zeros([numstorm, 5])

        # Select storms for year
        for n in range(numstorm):
            storm = random.randint(1, len(StormList) - 1)

            dur = StormList[storm, 1]  # Duration
            Rhigh = StormList[storm, 2]  # TWL
            period = StormList[storm, 4]  # Tp
            Rlow = StormList[storm, 6]  # Rlow

            stormTS[n, 0] = t
            stormTS[n, 1] = Rhigh / 10 - MHW
            stormTS[n, 2] = Rlow / 10 - MHW
            stormTS[n, 3] = period
            stormTS[n, 4] = round(dur / 2)  # Divided by two assuming TWL only for only half of storm

        # Save
        StormSeries = np.vstack([StormSeries, stormTS])
        np.save(datadir + output_filename, StormSeries)


def RunBatch(n):
    SimNum = len(rslr) * len(co) * len(slope)
    WidthData = np.zeros([8, len(rslr), len(co), len(slope)])
    Sim = 0 + (SimNum * n)

    # Make new parameters file for this set of simulations
    new_param_file = "Batch" + str(n + add) + "-barrier3d-parameters.yaml"
    shutil.copy("Input/Barrier3D/barrier3d-parameters.yaml", "Input/Barrier3D/" + new_param_file)

    # Make new parameters file for this set of simulations
    storm_file = "StormSeries_VCR_StormList_" + str(n + add) + ".npy"
    makeStormSeries(
        datadir='Input/Barrier3D/',
        storm_list_filename="StormList.npy",
        StormStart=2,
        TMAX=SimDur + 2,
        mean_storm=8,
        SD_storm=5.9,
        MHW=0.046,  # [dam]
        output_filename=storm_file,
    )

    for r in range(len(rslr)):
        for c in range(len(co)):
            for s in range(len(slope)):

                Sim += 1

                # Create an instance of the BMI class and set input parameter values
                barrierbmft = BarrierBMFT(
                    name="Batch",
                    time_step_count=SimDur,
                    relative_sea_level_rise=rslr[r],
                    reference_concentration=co[c],
                    slope_upland=slope[s],
                    storm_file=storm_file,  # Update new storm series
                    parameter_file=new_param_file
                )

                # Run simulation: Loop through time
                for time_step in range(int(barrierbmft.bmftc.dur)):

                    # Print time step to screen
                    # print("\r", "Time Step: ", time_step, "     Sim", Sim, "/", (SimNum * (n + 1)), "     (n", n, ", r", r, ", c", c, ")", end="")

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
                BBpond_w = widths[-1, 5] - widths[0, 5]
                MLpond_w = widths[-1, 6] - widths[0, 6]

                # Calculate shoreline change
                sc = int(barrierbmft.barrier3d.model.ShorelineChange)

                # Store in arrays
                WidthData[0, r, c, s] = barrier_w
                WidthData[1, r, c, s] = BBmarsh_w
                WidthData[2, r, c, s] = bay_w
                WidthData[3, r, c, s] = MLmarsh_w
                WidthData[4, r, c, s] = forest_w
                WidthData[5, r, c, s] = BBpond_w
                WidthData[6, r, c, s] = MLpond_w
                WidthData[7, r, c, s] = sc

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

    return WidthData


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
BBMarshPondWidth = np.zeros([Num, len(rslr), len(co), len(slope)])
MLMarshPondWidth = np.zeros([Num, len(rslr), len(co), len(slope)])

# ==================================================================================================================================================================================
# Run the BarrierBMFT batch

# Record start time
Time = time.time()

# num_cores = multiprocessing.cpu_count()
num_cores = 3
print('Cores: ' + str(num_cores))
print('Running...')

results = Parallel(n_jobs=num_cores)(delayed(RunBatch)(n) for n in range(Num))

# Reorganize
for N in range(Num):
    PS = results[N]
    BarrierWidth[N, :, :, :] = PS[0, :, :, :]
    BBMarshWidth[N, :, :, :] = PS[1, :, :, :]
    BayWidth[N, :, :, :] = PS[2, :, :, :]
    MLMarshWidth[N, :, :, :] = PS[3, :, :, :]
    ForestWidth[N, :, :, :] = PS[4, :, :, :]
    BBMarshPondWidth[N, :, :, :] = PS[5, :, :, :]
    MLMarshPondWidth[N, :, :, :] = PS[6, :, :, :]
    ShorelineChange[N, :, :, :] = PS[7, :, :, :]

# Save batch data arrays
np.save(directory + '/Widths_Barrier.npy', BarrierWidth)
np.save(directory + '/Widths_BBMarsh.npy', BBMarshWidth)
np.save(directory + '/Widths_Bay.npy', BayWidth)
np.save(directory + '/Widths_MLMarsh.npy', MLMarshWidth)
np.save(directory + '/Widths_Forest.npy', ForestWidth)
np.save(directory + '/Widths_BBMarshPond.npy', BBMarshPondWidth)
np.save(directory + '/Widths_MLMarshPond.npy', MLMarshPondWidth)
np.save(directory + '/ShorelineChange.npy', ShorelineChange)

# Print elapsed time of batch run
print()
BatchDuration = time.time() - Time
print()
print("Elapsed Time: ", BatchDuration, "sec")
