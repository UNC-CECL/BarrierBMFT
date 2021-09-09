"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 24 August 2021
"""
import bisect
import numpy as np
import math
import matplotlib.pyplot as plt
from yaml import full_load, dump

from barrier3d import Barrier3dBmi
from bmftc import Bmftc


def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)


def init_equal(bmftc_ML, bmftc_BB, datadir, input_file):
    """Initialize Barrier3D and set identical parameters equal to values in PyBMFT-C"""

    fid = datadir + input_file

    barrier3d = Barrier3dBmi()

    barrier3d.initialize(fid)

    # Check if PyBMFT-C Parameters Are Equal
    if bmftc_ML.dur != bmftc_BB.dur:
        print("time_step_count parameters not equal")
        return
    if bmftc_ML.RSLRi != bmftc_BB.RSLRi:
        print("relative_sea_level_rise parameters not equal")
        return
    if bmftc_ML.Co != bmftc_BB.Co:
        print("reference_concentration parameters not equal")
        return
    if bmftc_ML.slope != bmftc_BB.slope:
        print("slope_upland parameters not equal")
        return
    if bmftc_ML.bfo != bmftc_BB.bfo:
        print("bay_fetch_initial parameters not equal")
        return
    if bmftc_ML.mwo != bmftc_BB.mwo:
        print("marsh_width_initial parameters not equal")
        return
    if bmftc_ML.wind != bmftc_BB.wind:
        print("wind_speed parameters not equal")
        return
    if bmftc_ML.seagrass_on != bmftc_BB.seagrass_on:
        print("seagrass_on parameters not equal")
        return
    if bmftc_ML.tcr != bmftc_BB.tcr:
        print("critical_shear_mudflat parameters not equal")
        return

    # Equalize Barrier3D/PyBMFT-C Values of Identical Parameters
    set_yaml("TMAX", bmftc_ML.dur + 1, fid)  # [yrs] Duration of simulation
    barrier3d.model._RSLR = np.ones([len(barrier3d.model._RSLR) + 1]) * (bmftc_ML.RSLRi / 1000) / 10 # [m/yr] Relative sea-level rise rate, converted units
    barrier3d.model._BayDepth = bmftc_ML.Bay_depth[bmftc_ML.startyear - 1] / 10  # [yrs] Initial depth of bay

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
        # Mainland shoreline
        self._bmftc_ML = Bmftc(
            name="mainland",
            time_step_count=30,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=5,
            seagrass_on=False,
            forest_on=True,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
        )

        # Back-barier shoreline
        self._bmftc_BB = Bmftc(
            name="back-barrier",
            time_step_count=30,
            relative_sea_level_rise=3,
            reference_concentration=50,
            slope_upland=0.005,
            bay_fetch_initial=3000,
            wind_speed=5,
            seagrass_on=False,
            forest_on=True,
            critical_shear_mudflat=0.2,
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
        )

        # Specify data directory with Barrier3D initial conditions
        datadir = "Input/Barrier3D/"
        input_file = "barrier3d-parameters.yaml"

        # Initialize Barrier3D and set matching parameters equal
        self._barrier3d = init_equal(self._bmftc_ML, self._bmftc_BB, datadir, input_file)

        # ===========================================
        # Initialize variables
        self._coupled = coupling_on
        self._bay_unit_slope = bay_unit_slope  # Slope of underlying bay straigraphic unit
        self._bay_unit_depth = bay_unit_thickness  # [m] Depth (stratigraphic height) of underlying bay unit
        self._bay_unit_fine = bay_unit_fine  # Sand proportion, relative to mud, of underlying bay unit
        self._x_b_TS_ML = np.zeros([self._bmftc_ML.dur])
        self._x_b_TS_BB = np.zeros([self._bmftc_BB.dur])

    # ===========================================
    # Time Loop

    def update(self, time_step):
        """Update BarrierBMFT by one time step"""

        # Advance Barrier3D
        self._barrier3d.update()

        if self._coupled:

            # ===========================================
            # Adjust fetch and barrier-bay shoreline position in PyBMFT-C according to back-barrier shoreline change (i.e. overwash deposition)
            # delta_x_b_b3d = (self._barrier3d.model.x_b_TS[-1] - self._barrier3d.model.x_b_TS[-2]) * 10  # [m] Calculate back-barrier shoreline change from Barrier3D; (+) landward movement, (-) seaward movement
            # self._bmftc_BB._x_f = self._bmftc_BB.x_f - int(math.ceil(delta_x_b_b3d))  # Bay-barrier shoreline  # IR 9Sep21: ceil or round?
            # Add overwash elevation change to bmftc elev and recalculate x_f and x_m

            # ===========================================
            # Update PyBMFT-C transect elevation based on Barrier3D interior morphology
            barrier_transect = np.mean(self._barrier3d.model.InteriorDomain, axis=1) * 10  # Take average across alongshore dimension, convert to m (vertical dimension)
            boundary_elevation = self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_f - 1]
            x_b_effective = np.where(barrier_transect < boundary_elevation)[0][0]  # Cross shore location of marsh-barrier boundary
            barrier_interior_transect = barrier_transect[:x_b_effective]  # Transect of subaerial barrier interior from Barrier3D
            x = np.linspace(1, len(barrier_interior_transect) * 10, num=len(barrier_interior_transect) * 10)
            xp = np.linspace(1, len(barrier_interior_transect), num=len(barrier_interior_transect)) * 10
            barrier_interior_transect = np.interp(x, xp, barrier_interior_transect)  # Interpolate from dam to m (horizontal dimension)
            previous_interior_width = len(self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_f:])  # [m] Width of barrier interior from last timestep
            if len(barrier_interior_transect) > previous_interior_width:
                barrier_interior_transect = barrier_interior_transect[-previous_interior_width]  # Trim length to fit PyBMFT-C array length
            elif len(barrier_interior_transect) < previous_interior_width:
                add_zeros = np.zeros([previous_interior_width - len(barrier_interior_transect)])
                barrier_interior_transect = np.append(add_zeros, barrier_interior_transect)  # Add zeros (ocean) to fit PyBMFT-C array length
            self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_f:] = np.flip(barrier_interior_transect)  # Replace previous barrier interior with current (flipped cross-shore)

        print()
        print("     x_f_BB:", self._bmftc_BB.x_f)

        # Advance PyBMFT-C mainland and back-barrier marshes
        self._bmftc_ML.update()
        self._bmftc_BB.update()

        if self._coupled:

            """ 
            8 Sep 21
            To do:
                - update x_b (ML) and xm (BB) based on b3d bb shoreline change
                - incorporate effect of depositional thickness on veg survival (Walters and Kirwan 2016)
                - update bmftc "foreset" transect section based on barrier morphology
                - fix overwash deposition within bay
            """

            # ===========================================
            # Update fetch and marsh point locations from PyBMFT-C bay erosion/deposition processes

            # Calculate change in fetch from erosion of both marshes
            delta_fetch_ML = self._bmftc_ML.bfo - self._bmftc_ML.fetch[self._bmftc_ML.startyear + time_step - 1]  # [m] Mainland marsh
            delta_fetch_BB = self._bmftc_BB.bfo - self._bmftc_BB.fetch[self._bmftc_BB.startyear + time_step - 1]  # [m] Back-barrier marsh

            # Determine change in x_b location
            self._bmftc_ML._x_b = self._bmftc_ML.x_b - delta_fetch_BB
            self._bmftc_BB._x_b = self._bmftc_BB.x_b - delta_fetch_ML
            self._x_b_TS_ML[time_step] = self._bmftc_ML.x_b  # Save to array
            self._x_b_TS_BB[time_step] = self._bmftc_BB.x_b  # Save to array

            # Determine new fetch based on change in opposite marsh - both fetches should be exactly the same!
            self._bmftc_ML._bfo = self._bmftc_ML.bfo + delta_fetch_BB
            self._bmftc_BB._bfo = self._bmftc_BB.bfo + delta_fetch_ML
            self._bmftc_ML.fetch[self._bmftc_ML.startyear + time_step] = self._bmftc_ML.bfo  # Save to array
            self._bmftc_BB.fetch[self._bmftc_BB.startyear + time_step] = self._bmftc_BB.bfo  # Save to array

            # ===========================================
            # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
            self._barrier3d.model._BayDepth = np.mean([self._bmftc_ML.db, self._bmftc_BB.db]) / 10

            # ===========================================
            # Adjust bay routing width in Barrier3D if fetch in PyBMFTC becomes narrower than 100 m
            if self._bmftc_ML.bfo < 100:
                self._barrier3d.model.bay_routing_width = int(round(self._bmftc_ML.bfo / 10))

            # ===========================================
            # Adjust proportion of fine sediment ouctropping on Barrier3D shoreface according to bay stratigraphy from PyBMFT-C

            # Find thickness of 100% fine sediment deposited since model start
            percent_fine = 1  # IR 24Aug21 To Do: Determine and account for coarse overwash sediment deposition
            if self._bmftc_ML.x_b < 0:
                bay_dep_thickness = 0  # Assume eroding only barrier material on shoreface (i.e., 100% sand)
            else:
                bay_dep_thickness = (self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step, self._bmftc_ML.x_b] - self._bmftc_ML.elevation[0, self._bmftc_ML.x_b]) * percent_fine  # [m] Sums total deposition since model start

            # Find thickness of underlying bay unit
            bay_unit_thickness = self._bay_unit_depth - self._bmftc_ML.x_b * self._bay_unit_slope  # [m]

            # Calculate total effective thickness of underlying bay unit and newly deposited bay sediment
            fine_thickness = bay_unit_thickness * self._bay_unit_fine + bay_dep_thickness  # [m]
            if fine_thickness + self._bmftc_ML.db > self._barrier3d.model.DShoreface * 10:  # If the base of the fine unit extends below the shoreface toe
                fine_thickness = (self._barrier3d.model.DShoreface * 10) - self._bmftc_ML.db  # Thickness of estuarine unit impacting migration can't be greater than the depth of shoreface erosion (i.e., shoreface toe depth)

            # Set fine sediment fraction on the shoreface in Barrier3D (follows Nienhuis & Lorenzo-Trueba, 2019, Geosci. Model Dev.)
            self._barrier3d.model._shoreface_fine_fraction = fine_thickness / (fine_thickness + self._bmftc_ML.db + self._barrier3d.model.h_b_TS[-1] * 10)

            # ===========================================
            # Add overwash deposition in the bay from Barrier3D to bay of PyBMFT-C -- IR 17 Aug 21: Currently does nothing...
            overwash_bay_deposition = barrier_transect[np.where(barrier_transect <= 0)[0][0]: np.where(barrier_transect <= 0)[0][-1] + 1]  # [m] Depths of subaqueous cells
            overwash_bay_deposition -= self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step - 1, 0: len(overwash_bay_deposition)]  # Extracts overwash bay deposition from subaqueous portion of B3D for this time step
            overwash_bay_deposition[overwash_bay_deposition < 0] = 0  # Can't have negative deposition
            self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step, 0: len(overwash_bay_deposition)] += overwash_bay_deposition  # IR 17 Aug 21: PyBMFTC erases this addition each time-step, so this currently does effectively nothing


    @property
    def bmftc(self):
        return self._bmftc_ML

    @property
    def barrier3d(self):
        return self._barrier3d

    @property
    def x_b_TS_bmft(self):
        return self._x_b_TS_bmft
