"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 17 August 2021
"""

import numpy as np
import math
from yaml import full_load, dump

from barrier3d import Barrier3dBmi
from bmftc import Bmftc


def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)


def init_equal(bmftc, datadir, input_file):
    """Initialize Barrier3D and set identical parameters equal to values in PyBMFT-C"""

    fid = datadir + input_file

    barrier3d = Barrier3dBmi()

    barrier3d.initialize(fid)

    # Equalize Values of Identical Parameter
    set_yaml("TMAX", bmftc.dur + 1, fid)  # [yrs] Duration of simulation
    barrier3d.model._RSLR = np.ones([len(barrier3d.model._RSLR) + 1]) * (bmftc.RSLRi / 1000) / 10 # [m/yr] Relative sea-level rise rate, converted units
    barrier3d.model._BayDepth = bmftc.Bay_depth[bmftc.startyear - 1] / 10  # [yrs] Initial depth of bay

    return barrier3d


class BarrierBMFT:
    """
    Couples Barrier3D and PyBMFT-C
    """

    def __init__(
            self,
            coupling_on=True,
    ):
        """ Initialize Barrier3D and PyBMFT-C """

        self._coupled = coupling_on

        # ===========================================
        # Initialize Model Classes
        # ===========================================

        # Initialize PyBMFT-C
        self._bmftc = Bmftc(
            name="default",
            time_step=1,
            time_step_count=30,
            relative_sea_level_rise=1,
            reference_concentration=50,
            slope_upland=0.005,
        )

        # Specify data directory with Barrier3D initial conditions
        datadir = "Input/Barrier3D/"
        input_file = "barrier3d-parameters.yaml"

        # Initialize Barrier3D and set matching parameters equal
        self._barrier3d = init_equal(self._bmftc, datadir, input_file)

    # ===========================================
    # Time Loop
    # ===========================================

    def update(self, time_step):
        """Update BarrierBMFT by one time step"""

        # Advance Barrier3D
        self._barrier3d.update()

        # Advance PyBMFT-C
        self._bmftc.update()

        if self._coupled:

            # Calculate back-barrier shoreline change from Barrier3D
            delta_x_b = (self._barrier3d.model._x_b_TS[-2] - self._barrier3d.model._x_b_TS[-1]) * 10  # Convert from dam to m

            # Extract overwash bay deposition from subaqueous portion of B3D
            barrier_transect = np.mean(self._barrier3d.model._InteriorDomain, axis=1) * 10
            overwash_bay_deposition = barrier_transect[np.where(barrier_transect <= 0)[0][0]: np.where(barrier_transect <= 0)[0][-1] + 1]  # [m] Depths of subaqueous cells
            overwash_bay_deposition -= self._bmftc.elevation[self._bmftc.startyear + time_step - 1, 0: len(overwash_bay_deposition)]  # Get deposition for this time step
            overwash_bay_deposition[overwash_bay_deposition < 0] = 0  # Can't have negative deposition

            # Adjust fetch in PyBMFT-C according to back-barrier shoreline change
            self._bmftc._bfo = self._bmftc.bfo + int(round(delta_x_b))  # IR 19 Aug 21: Is this an OK way to do it? Or should the left (i.e. first) cell in the bay be deleted in elevation, mineral, organic_dep, etc???

            # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
            self._barrier3d.model._BayDepth = self._bmftc.db / 10

            # Add overwash deposition in the bay from Barrier3D to bay of PyBMFT-C
            self._bmftc.elevation[self._bmftc.startyear + time_step, 0: len(overwash_bay_deposition)] += overwash_bay_deposition  # IR 17 Aug 21: PyBMFTC erases this addition each time-step, so this currently does effectively nothing

            # Adjust bay routing width in Barrier3D if fetch in PyBMFTC becomes narrower than 100 m
            if self._bmftc.bfo < 100:
                self._barrier3d.model.bay_routing_width = int(round(self._bmftc.bfo / 10))


    @property
    def bmftc(self):
        return self._bmftc

    @property
    def barrier3d(self):
        return self._barrier3d
