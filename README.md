# BarrierBMFT
Coupled Barrier-Bay-Marsh-Forest Transect Model

## About
*BarrierBMFT* is an coupled model framework for exploring morphodynamic interactions across components of the entire
coastal barrier system, from the ocean shoreface to the mainland forest. The model frameqork couples *Barrier3D* (Reeves 
et al., 2021, a spatially explicit model of barrier evolution, with the Python version of the Coastal Landscape Transect 
model (*CoLT*; Valentine et al., 2023), known as *PyBMFT-C* (Bay-Marsh-Forest Transect Model with Carbon). In the coupled 
model framework, BarrierBMFT, two *PyBMFT-C* simulations drive evolution of back-barrier marsh, bay, mainland marsh, and 
forest ecosystems, and a *Barrier3D* simulation drives evolution of barrier and back-barrier marsh ecosystems. As these model
components simultaneously advance, they dynamically evolve together by sharing information annually to capture the effects 
of key cross-landscape couplings. *BarrierBMFT* contains no new governing equations or parameterizations itself, but rather is 
a framework for trading information between *Barrier3D* and *PyBMFT-C*. Detailed desrciptions of *BarrierBMFT* and the
coupled models involved can be found in publications listed below under *Further Reading*.

_Copyright (C) 2021 Ian R.B. Reeves (principal developer) licensed under the GNU General Public License v3.0_

## Requirements
BarrierBMFT requires Python 3, and the libraries listed in the project's `requirements.txt` file.

## Installation

First, download the source code for *BarrierBMFT*, *PyBMFT-C*, and *Barrier3D* into separate subdirectories within the 
same project directory. To get the source code you can either clone the repositories with git:

    git clone git@github.com/UNC-CECL/barrierbmft
    git clone git@github.com/UNC-CECL/pybmft-c
    git clone git@github.com/UNC-CECL/barrier3d

or download zip files:

    https://github.com/UNC-CECL/BarrierBMFT/archive/refs/heads/main.zip
    https://github.com/UNC-CECL/Barrier3D/archive/refs/heads/master.zip
    https://github.com/UNC-CECL/PyBMFT-C/archive/refs/heads/main.zip

You should now have directories organized as:

    Your_Project_Directory __________ BarrierBMFT
                                  |
                                  |__ Barrier3D
                                  |
                                  |__ PyBMFT-C

Lastly, for each of the of *BarrierBMFT*, *PyBMFT-C*, and *Barrier3D* directories, run the following from their top-level folders 
(the ones that contain `setup.py`) to install each of them into the current environment:

    pip install -e .

## Input Files & Parameters

#### BarrierBMFT:
A main set of commonly-manipulated parameters can be adjusted in the initializarion of the *BarrierBMFT* class in
the run script (see example below).

#### Barrier3D:
Input files are located in the `BarrierBMFT/Input/Barrier3D` directory. *Barrier3D* parameter values can be
adjusted in the `barrier3d-parameters.yaml` file; a decription of each variable is available in 
`Barrier3D/barrier3d/configuration.py`. In the initialization of *BarrierBMFT*, the simulation duration,
relative sea-level rise (RSLR) rate, and bay depth are set according the values of the same parameters in *PyBMFT-C*,
therefore the value of these parameters within the `barrier3d-parameters.yaml` file has no effect on the simulation.

#### PyBMFT-C:
Input files are located in the `BarrierBMFT/Input/PyBMFT-C` directory. *PyBMFT-C* parameter values can be set
in the initialization of the *PyBMFT-C* classes within `barrierbmft.py`. Note that because *BarrierBMFT* runs two separate
instances ofthe *PyBMFT-C* model during a simulation, both instances can be set with different parameters, though a warning
will be printed to screen if specific parameter values do not match between both instances in the initialization.

## Example Simulation

The following describes the approach for running a basic *BarrierBMFT* simulation. For a complete example run 
script, see `BarrierBMFT/scripts/run_BarrierBMFT.py`. 

To run *BarrierBMFT*, first set the working direcory to the main `BarrierBMFT` directory.

Then, import *BarrierBMFT* and dependencies:
    
    from barrierbmft.barrierbmft import BarrierBMFT
    import numpy as np
    import matplotlib.pyplot as plt

Next, initialize an instance of the *BarrierBMFT* class. Here, you can set certain parameter values such as the RSLR 
rate (mm/yr) or the external suspended sediment supply (i.e., reference concentration; mg/L):

    # Create an instance of the class
    barrierbmft = BarrierBMFT(relative_sea_level_rise=12, reference_concentration=60)
    print(barrierbmft.name)

Next, run the simulation by progressing through time:

    # Loop through time
    for time_step in range(int(barrierbmft.bmftc.dur)):

        # Run time step
        barrierbmft.update(time_step)
    
        # Check for breaks (e.g., barrier drowning)
        if barrierbmft.BMFTC_Break or barrierbmft.Barrier3D_Break:
            break

Once the simulation finishes, we can plot some results. For example, the full barrier-to-mainland-forest transect every
at simulation end:

    # Time index at simulatione end
    t = barrierbmft.bmftc_BB.dur - 1
    
    # Back-barrier transect
    BB_transect = barrierbmft.bmftc_BB.elevation[barrierbmft.bmftc_BB.startyear + t - 1, int(barrierbmft.bmftc_BB.Marsh_edge[barrierbmft.bmftc_BB.startyear + t]):]  # Trim off bay of ML transect (use identical bay in BB transect)
    BB_transect = np.flip(BB_transect)
    
    # Mainland transect
    if barrierbmft.x_b_TS_ML[t] < 0:  # If far end of bay extends beyond end of transect, append additional bay cells for visualization
        ML_transect = np.append(np.ones([abs(int(barrierbmft.x_b_TS_ML[t]))]) * barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, 1], barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, :])
    elif barrierbmft.x_b_TS_ML[t] > 0:  # Far end of bay to end of mainland forest
        ML_transect = barrierbmft.bmftc_ML.elevation[barrierbmft.bmftc_ML.startyear + t - 1, int(barrierbmft.x_b_TS_ML[t]):]

    # Combine transects and plot
    whole_transect = np.append(BB_transect, ML_transect)
    plt.figure()
    plt.plot(whole_transect)
    plt.xlabel("Distance")
    plt.ylabel("Elevation [m MSL]")

Or, plot the change in extent over time for each ecosystem:

    # Landscape extent
    widths = barrierbmft.LandscapeTypeWidth_TS
    barrier = widths[:, 0]
    BBmarsh = widths[:, 1]
    bay = widths[:, 2]
    MLmarsh = widths[:, 3]
    forest = widths[:, 4]
    BBpond = widths[:, 5]
    MLpond = widths[:, 6]
    barrier_marsh = barrier + BBmarsh

    plt.figure()
    
    # Barrier
    plt.subplot(2, 3, 1)
    plt.plot(barrier, c="red")
    plt.ylabel("Barrier [m]")
    
    # Back-Barrier Marsh
    plt.subplot(2, 3, 2)
    plt.plot(BBmarsh, c="red")
    plt.plot(BBpond, c="red", linestyle='dashed')
    plt.ylabel("Back-Barrier Marsh [m]")
    
    # Barrier + Marsh
    plt.subplot(2, 3, 3)
    plt.plot(barrier_marsh, c="red")
    plt.ylabel("Barrier + Marsh [m]")
    
    # Bay
    plt.subplot(2, 3, 4)
    plt.plot(bay, c="red")
    plt.ylabel("Bay [m]")
    
    # Mainland marsh
    plt.subplot(2, 3, 5)
    plt.plot(MLmarsh, c="red")
    plt.plot(MLpond, c="red", linestyle='dashed')
    plt.ylabel("Mainland Marsh [m]")
    plt.xlabel("Time [yr]")
    
    # Forest
    plt.subplot(2, 3, 6)
    plt.plot(forest, c="red")
    plt.ylabel("Forest [m]")
    plt.tight_layout()

## Further Reading

#### BarrierBMFT
    Reeves, I.R.B., Moore, L.J., Valentine, K., Fagherazzi, S., & Kirwan, M.L. (in review). Cross-landscape couplings influence 
    ecosystem extent in coastal barrier systems.
#### Barrier3D
    Reeves, I.R.B., Moore, L.J., Murray, A.B., Anarde, K.A., & Goldstein, E.B. (2021). Dune dynamics drive discontinuous 
    barrier retreat. Geophysical Research Letters, 48(13), e2021GL092958. https://doi.org/10.1029/2021GL092958.
#### PyBMFT-C / CoLT
    Valentine, K., Herbert, E. R., Walters, D. C., Chen, Y., Smith, A. J., & Kirwan, M. L. (2023). Climate-driven tradeoffs between 
    landscape connectivity and the maintenance of the coastal carbon sink. Nature Communications, 14, 1137. 
    https://doi.org/10.1038/s41467-023-36803-7.
