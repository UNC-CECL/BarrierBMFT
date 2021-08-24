"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 24 August 2021
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
            bay_unit_slope=0.0004,  # Hog Bay (Finkelstein & Ferland, 1987)
            bay_unit_thickness=4.3,  # Hog Bay (Finkelstein & Ferland, 1987)
            bay_unit_fine=0.5,  # Hog Bay (Finkelstein & Ferland, 1987)
    ):
        """ Initialize Barrier3D and PyBMFT-C """

        # ===========================================
        # Initialize Model Classes

        # Initialize PyBMFT-C
        self._bmftc = Bmftc(
            name="default",
            time_step=1,
            time_step_count=50,
            relative_sea_level_rise=1,
            reference_concentration=50,
            slope_upland=0.005,
            seagrass_on=False,
        )

        # Specify data directory with Barrier3D initial conditions
        datadir = "Input/Barrier3D/"
        input_file = "barrier3d-parameters.yaml"

        # Initialize Barrier3D and set matching parameters equal
        self._barrier3d = init_equal(self._bmftc, datadir, input_file)

        # ===========================================
        # Initialize variables
        self._coupled = coupling_on
        self._bay_unit_slope = bay_unit_slope  # Slope of underlying bay straigraphic unit
        self._bay_unit_depth = bay_unit_thickness  # [m] Depth (stratigraphic height) of underlying bay unit
        self._bay_unit_fine = bay_unit_fine  # Sand proportion, relative to mud, of underlying bay unit
        self._x_b_TS_bmft = np.zeros([self._bmftc.dur])

    # ===========================================
    # Time Loop

    def update(self, time_step):
        """Update BarrierBMFT by one time step"""

        # Advance Barrier3D
        self._barrier3d.update()

        # Advance PyBMFT-C
        self._bmftc.update()

        if self._coupled:

            # # Calculate back-barrier shoreline change from Barrier3D
            # delta_x_b = (self._barrier3d.model.x_b_TS[-1] - self._barrier3d.model.x_b_TS[-2]) * 10  # Convert from dam to m; (+) landward movement, (-) seaward movement

            # # Extract overwash bay deposition from subaqueous portion of B3D
            # barrier_transect = np.mean(self._barrier3d.model.InteriorDomain, axis=1) * 10
            # overwash_bay_deposition = barrier_transect[np.where(barrier_transect <= 0)[0][0]: np.where(barrier_transect <= 0)[0][-1] + 1]  # [m] Depths of subaqueous cells
            # overwash_bay_deposition -= self._bmftc.elevation[self._bmftc.startyear + time_step - 1, 0: len(overwash_bay_deposition)]  # Get deposition for this time step
            # overwash_bay_deposition[overwash_bay_deposition < 0] = 0  # Can't have negative deposition

            # ===========================================
            # Adjust fetch and barrier-bay shoreline position in PyBMFT-C according to back-barrier shoreline change
            delta_x_b = (self._barrier3d.model.x_b_TS[-1] - self._barrier3d.model.x_b_TS[-2]) * 10  # [m] Calculate back-barrier shoreline change from Barrier3D; (+) landward movement, (-) seaward movement
            self._bmftc._bfo = self._bmftc.bfo - int(round(delta_x_b))  # Fetch
            self._bmftc._x_b = self._bmftc.x_b + int(round(delta_x_b))  # Bay-barrier shoreline
            self._x_b_TS_bmft[time_step - 1] = self._bmftc.x_b  # Store in array

            # ===========================================
            # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
            self._barrier3d.model._BayDepth = self._bmftc.db / 10

            # ===========================================
            # Adjust bay routing width in Barrier3D if fetch in PyBMFTC becomes narrower than 100 m
            if self._bmftc.bfo < 100:
                self._barrier3d.model.bay_routing_width = int(round(self._bmftc.bfo / 10))

            # ===========================================
            # Adjust proportion of fine sediment ouctropping on Barrier3D shoreface according to bay stratigraphy from PyBMFT-C

            # Find thickness of 100% fine sediment deposited since model start
            percent_fine = 1  # IR 24Aug21 To Do: Determine and account for coarse overwash sediment deposition
            if self._bmftc.x_b < 0:
                bay_dep_thickness = 0  # Assume eroding only barrier material on shoreface (i.e., 100% sand)
            else:
                bay_dep_thickness = (self._bmftc.elevation[self._bmftc.startyear + time_step, self._bmftc.x_b] - self._bmftc.elevation[0, self._bmftc.x_b]) * percent_fine  # [m] Sums total deposition since model start

            # Find thickness of underlying bay unit
            bay_unit_thickness = self._bay_unit_depth - self._bmftc.x_b * self._bay_unit_slope  # [m]

            # Calculate total effective thickness of underlying bay unit and newly deposited bay sediment
            fine_thickness = bay_unit_thickness * self._bay_unit_fine + bay_dep_thickness  # [m]
            if fine_thickness + self._bmftc.db > self._barrier3d.model.DShoreface * 10:  # If the base of the fine unit extends below the shoreface toe
                fine_thickness = (self._barrier3d.model.DShoreface * 10) - self._bmftc.db  # Thickness of estuarine unit impacting migration can't be greater than the depth of shoreface erosion (i.e., shoreface toe depth)

            # Set fine sediment fraction on the shoreface in Barrier3D (follows Nienhuis & Lorenzo-Trueba, 2019, Geosci. Model Dev.)
            self._barrier3d.model._shoreface_fine_fraction = fine_thickness / (fine_thickness + self._bmftc.db + self._barrier3d.model.h_b_TS[-1] * 10)

            # ===========================================
            # Add overwash deposition in the bay from Barrier3D to bay of PyBMFT-C
            barrier_transect = np.mean(self._barrier3d.model.InteriorDomain, axis=1) * 10
            overwash_bay_deposition = barrier_transect[np.where(barrier_transect <= 0)[0][0]: np.where(barrier_transect <= 0)[0][-1] + 1]  # [m] Depths of subaqueous cells
            overwash_bay_deposition -= self._bmftc.elevation[self._bmftc.startyear + time_step - 1, 0: len(overwash_bay_deposition)]  # Extracts overwash bay deposition from subaqueous portion of B3D for this time step
            overwash_bay_deposition[overwash_bay_deposition < 0] = 0  # Can't have negative deposition
            self._bmftc.elevation[self._bmftc.startyear + time_step, 0: len(overwash_bay_deposition)] += overwash_bay_deposition  # IR 17 Aug 21: PyBMFTC erases this addition each time-step, so this currently does effectively nothing


    @property
    def bmftc(self):
        return self._bmftc

    @property
    def barrier3d(self):
        return self._barrier3d

    @property
    def x_b_TS_bmft(self):
        return self._x_b_TS_bmft
