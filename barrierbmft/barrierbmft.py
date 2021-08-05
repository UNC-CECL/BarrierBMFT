"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 6 July 2021
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
    set_yaml("TMAX", bmftc.dur + 1, fid)
    # set_yaml("RSLR_const", bmftc.RSLRi / 1000, fid)
    # set_yaml("BayDepth", bmftc.Bay_depth[bmftc.startyear], fid)
    # barrier3d.model._TMAX = bmftc.dur  # [yrs] Duration of simulation
    barrier3d.model._RSLR = np.ones([len(barrier3d.model._RSLR) + 1]) * (bmftc.RSLRi / 1000) / 10 # [m/yr] Relative sea-level rise rate, converted units
    barrier3d.model._BayDepth = bmftc.Bay_depth[bmftc.startyear - 1] / 10  # [yrs] Initial depth of bay

    return barrier3d


class BarrierBMFT:
    """
    Couples Barrier3D and PyBMFT-C
    """

    def __init__(
            self,
    ):
        """ Initialize Barrier3D and PyBMFT-C """

        # ===========================================
        # Initialize Model Classes
        # ===========================================

        # Initialize PyBMFT-C
        self._bmftc = Bmftc(
            name="default",
            time_step=1,
            time_step_count=70,
            relative_sea_level_rise=1,
            reference_concentration=10,
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

        # Only run Barrier3D after PyBMFTC has completed spin up
        # if time_step >= self._bmftc.startyear:

        # Advance Barrier3D
        self._barrier3d.update()

        # Calculate back-barrier shoreline change from Barrier3D
        delta_x_b = (self._barrier3d.model._x_b_TS[-2] - self._barrier3d.model._x_b_TS[-1]) * 10  # Convert from dam to m

        # Advance PyBMFT-C
        self._bmftc.update()

        # Adjust fetch in PyBMFT-C according to back-barrier shoreline change
        # self._bmftc._bfo = self._bmftc.bfo + int(round(delta_x_b))

        # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
        self._barrier3d.model._BayDepth = self._bmftc.db / 10


    @property
    def bmftc(self):
        return self._bmftc

    @property
    def barrier3d(self):
        return self._barrier3d




