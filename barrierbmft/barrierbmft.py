"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 4 April 2022
"""

import numpy as np
import math
import warnings
from yaml import full_load, dump
from barrier3d import Barrier3dBmi
from bmftc import Bmftc
import matplotlib.pyplot as plt

warnings.simplefilter("ignore", category=RuntimeWarning)


def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)
    doc[var_name] = new_vals
    with open(file_name, "w") as f:
        dump(doc, f)


def init_equal(bmftc_ML, bmftc_BB, datadir, parameter_file, storm_file):
    """Initialize Barrier3D and set identical parameters equal to values in PyBMFT-C"""

    fid = datadir + parameter_file

    # Set storm series and duration
    set_yaml("storm_file", storm_file, fid)
    set_yaml("TMAX", bmftc_ML.dur + 1, fid)  # [yrs] Duration of simulation

    # Initialize
    barrier3d = Barrier3dBmi()
    barrier3d.initialize(config_file=fid)

    # Check if PyBMFT-C Parameters Are Equal
    if bmftc_ML.dur != bmftc_BB.dur:
        print(" - time_step_count parameters not equal")
        return
    if bmftc_ML.RSLRi != bmftc_BB.RSLRi:
        print(" - relative_sea_level_rise parameters not equal")
        return
    if bmftc_ML.Co != bmftc_BB.Co:
        print(" - reference_concentration parameters not equal")
        return
    if bmftc_ML.slope != bmftc_BB.slope:
        print(" - slope_upland parameters not equal")
    if bmftc_ML.bfo != bmftc_BB.bfo:
        print(" - bay_fetch_initial parameters not equal")
        return
    if bmftc_ML.forest_width_initial != bmftc_BB.forest_width_initial:
        print(" - forest_width_initial parameters not equal")
        return
    if bmftc_ML.forest_width_initial_fixed != bmftc_BB.forest_width_initial_fixed:
        print(" - forest_width_initial_fixed parameters not equal")
        return
    if bmftc_ML.mwo != bmftc_BB.mwo:
        print(" - marsh_width_initial parameters not equal")
    if bmftc_ML.wind != bmftc_BB.wind:
        print(" - wind_speed parameters not equal")
        return
    if bmftc_ML.tcr != bmftc_BB.tcr:
        print(" - critical_shear_mudflat parameters not equal")
        return

    # Equalize Barrier3D/PyBMFT-C Values of Identical Parameters
    barrier3d.model._TMAX = bmftc_ML.dur + 1  # [yrs] Duration of simulation
    barrier3d.model._RSLR = np.ones([len(barrier3d.model.RSLR) + 1]) * (bmftc_ML.RSLRi / 1000) / 10  # [m/yr] Relative sea-level rise rate, converted units
    barrier3d.model._BayDepth = bmftc_ML.Bay_depth[bmftc_ML.startyear - 1] / 10  # [yrs] Initial depth of bay

    return barrier3d


class BarrierBMFT:
    """
    Couples Barrier3D and PyBMFT-C
    """

    def __init__(
            self,
            name="BarrierBMFT",
            time_step_count=50,
            relative_sea_level_rise=12,
            reference_concentration=60,
            slope_upland=0.005,
            storm_file="StormTimeSeries_1000yr.npy",  # "StormSeries_VCR_Berm1pt9m_Slope0pt04.npy",
            parameter_file="barrier3d-parameters.yaml",
    ):
        """ Initialize Barrier3D and PyBMFT-C """

        # ===========================================
        # Initialize Model Classes

        # Initialize PyBMFT-C
        # Mainland shoreline
        self._bmftc_ML = Bmftc(
            name="mainland",
            time_step_count=time_step_count,
            relative_sea_level_rise=relative_sea_level_rise,
            reference_concentration=reference_concentration,
            slope_upland=slope_upland,
            bay_fetch_initial=5000,
            forest_width_initial_fixed=False,
            forest_width_initial=5000,
            wind_speed=6,
            forest_on=True,
            filename_equilbaydepth="Input/PyBMFT-C/Equilibrium Bay Depth.mat",  # "Input/PyBMFT-C/EquilibriumBayDepth_f3000_w5.mat",
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50_width500.mat",  # "Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
            marsh_width_initial=500,
        )

        # Back-barier shoreline
        self._bmftc_BB = Bmftc(
            name="back-barrier",
            time_step_count=time_step_count,
            relative_sea_level_rise=relative_sea_level_rise,
            reference_concentration=reference_concentration,
            slope_upland=slope_upland,
            bay_fetch_initial=5000,
            forest_width_initial_fixed=False,
            forest_width_initial=5000,  # 5000 accomodates 250 yrs at R=15 and S=0.001
            wind_speed=6,
            forest_on=False,
            filename_equilbaydepth="Input/PyBMFT-C/Equilibrium Bay Depth.mat",  # "Input/PyBMFT-C/EquilibriumBayDepth_f3000_w5.mat",
            filename_marshspinup="Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50_width500.mat",  # "Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50.mat",
            marsh_width_initial=500,
        )

        # Initialize Barrier3D
        datadir = "Input/Barrier3D/"  # Specify data directory
        self._parameter_file = parameter_file  # Specify file with Barrier3D initial conditions
        self._barrier3d = init_equal(self._bmftc_ML, self._bmftc_BB, datadir, self._parameter_file, storm_file)  # Initialize Barrier3D and set matching parameters equal

        # ===========================================
        # Initialize break variables
        self._Barrier3D_Break = False
        self._BMFTC_Break = False

        # ===========================================
        # Add initial barrier topography from Barrier3D to initial "forest" (i.e., subaerial) portion of PyBMFT-C transect
        b3d_transect = np.mean(self._barrier3d.model.InteriorDomain, axis=1) * 10  # Take average across alongshore dimension, convert to m (vertical dimension)
        x = np.linspace(1, len(b3d_transect) * 10, num=len(b3d_transect) * 10)
        xp = np.linspace(1, len(b3d_transect), num=len(b3d_transect)) * 10
        xp = xp - 5
        b3d_transect = np.interp(x, xp, b3d_transect)  # Interpolate from dam to m (horizontal dimension)
        x_f = np.where(b3d_transect < (self._barrier3d.model.SL * 10))[0][0] - 1  # [m] Distance of first interior (subaerial) cell from B3D ocean shoreline (excluding dunes/beach)
        b3d_transect = b3d_transect[:x_f]
        b3d_transect = np.flip(b3d_transect)
        b3d_transect = b3d_transect + self._bmftc_BB.msl[self._bmftc_BB.startyear - 1] + self._bmftc_BB.amp + (self._bmftc_BB.RSLRi / 1000)  # Convert vertical datums

        # Adjust size of Barrier3D topo to fit PyBMFT-C "forest" section
        BB_forest_len = len(self._bmftc_BB.elevation[self._bmftc_BB.startyear, self._bmftc_BB.x_f:])
        if len(b3d_transect) > BB_forest_len:
            subtract = len(b3d_transect) - BB_forest_len
            b3d_transect = b3d_transect[:-subtract]
        elif len(b3d_transect) < BB_forest_len:
            add = np.ones([BB_forest_len - len(b3d_transect)]) * (self._bmftc_BB.msl[self._bmftc_BB.startyear] + self._bmftc_BB.amp)
            b3d_transect = np.append(b3d_transect, add)

        # Replace initial subaerial elevation in PyBMFT-C with Barrier3D initial barrier elevation
        self._bmftc_BB.elevation[self._bmftc_BB.startyear - 1, self._bmftc_BB.x_f:] = b3d_transect  # Replace!

        # ===========================================
        # Initialize variables
        self._name = name
        self._x_b_TS_ML = np.zeros([self._bmftc_ML.dur])
        self._x_b_TS_BB = np.zeros([self._bmftc_BB.dur])
        self._LandscapeTypeWidth_TS = np.zeros([self._bmftc_BB.dur, 7])
        self._bay_overwash_carryover = 0  # [m^3] Volume of overwash deposition into back-barrier bay from previous year that did not fill new cell up to sea level; is added to overwash bay dep in following year
        initial_BB_subaerial_width = self._bmftc_BB.B - self._bmftc_BB.x_f
        self._x_s_offset = initial_BB_subaerial_width - (self._barrier3d.model.InteriorWidth_AvgTS[-1] * 10)  # Initial location of B in PyBMFT-C relative to x_s_initial in Barrier3D
        self._cumul_len_change = [0]
        self._delta_fetch_BB_TS = []
        self._delta_fetch_ML_TS = []
        self._OWspread = 0  # [%] Percentage of overwash past marsh edge that is spread across bay

    # =======================================================================================================================================================================================================================================
    # =======================================================================================================================================================================================================================================
    # Time Loop

    def update(self, time_step):

        """Update BarrierBMFT by one time step"""

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Advance PyBMFT-C mainland and back-barrier marshes
        self._bmftc_ML.update()
        self._bmftc_BB.update()

        # Set OCb equal
        self._bmftc_ML.OCb[self._bmftc_ML.startyear + time_step] = (self._bmftc_ML.OCb[self._bmftc_ML.startyear + time_step] + self._bmftc_BB.OCb[self._bmftc_BB.startyear + time_step]) / 2
        self._bmftc_BB.OCb[self._bmftc_BB.startyear + time_step] = (self._bmftc_ML.OCb[self._bmftc_ML.startyear + time_step] + self._bmftc_BB.OCb[self._bmftc_BB.startyear + time_step]) / 2

        # Check if marsh has completely drowned or basin is completely full
        if self._bmftc_ML.drown_break == 1 or self._bmftc_BB.drown_break == 1:
            self._bmftc_ML._dur = time_step
            self._bmftc_ML._endyear = self._bmftc_ML.startyear + time_step
            self._bmftc_BB._dur = time_step
            self._bmftc_BB._endyear = self._bmftc_BB.startyear + time_step
            self._BMFTC_Break = True
            print("PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full")
            return  # If so, end simulation

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Update fetch and marsh point locations from PyBMFT-C bay erosion/deposition processes

        # Calculate change in fetch from erosion of both marshes
        delta_fetch_ML = self._bmftc_ML.bfo - self._bmftc_ML.fetch[self._bmftc_ML.startyear + time_step - 1]  # [m] Mainland marsh
        delta_fetch_BB = self._bmftc_BB.bfo - self._bmftc_BB.fetch[self._bmftc_BB.startyear + time_step - 1]  # [m] Back-barrier marsh

        self._delta_fetch_BB_TS.append(delta_fetch_BB)
        self._delta_fetch_ML_TS.append(delta_fetch_ML)

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

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
        self._barrier3d.model._BayDepth = np.mean([self._bmftc_ML.db, self._bmftc_BB.db]) / 10

        # Adjust bay routing width in Barrier3D if fetch in PyBMFTC becomes narrower than 100 m
        if self._bmftc_ML.bfo < 100:
            self._barrier3d.model._bay_routing_width = int(math.floor(self._bmftc_ML.bfo / 10))
        else:
            self._barrier3d.model._bay_routing_width = 15

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Add marsh from PyBMFT-C to Barrier3D

        # Extract and convert marsh elevation from PyBMFT-C
        marsh_transect = self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m: self._bmftc_BB.x_f + 1]  # Marsh elevation from PyBMFT-C

        if len(marsh_transect) >= 1:
            len_marsh_transect = 10 * ((len(marsh_transect) + 5) // 10)  # Cross-shore length of marsh rounded to nearest dam
            self._cumul_len_change.append(self._cumul_len_change[-1] + (len_marsh_transect - len(marsh_transect)))
            x = np.linspace(1, len(marsh_transect) / 10, num=int((len_marsh_transect / 10)))
            xp = np.linspace(1, len(marsh_transect) / 10, num=int(len(marsh_transect)))
            marsh_transect = np.interp(x, xp, marsh_transect)  # Interpolate marsh elevation from m to dam in the horizontal dimension
            marsh_transect = marsh_transect - (self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step - 1] + self._bmftc_BB.amp)  # Make marsh elevation relative to MHW datum
            marsh_transect = marsh_transect / 10  # Convert from m to dam in the vertial dimension
            marsh_transect = np.flip(marsh_transect)

        StartDomainWidth = np.shape(self._barrier3d.model.InteriorDomain)[0]  # Width of interior domain from last time step

        # Find barrier interior widths for each dam alongshore
        InteriorWidth = [0] * self._barrier3d.model.BarrierLength
        for bl in range(self._barrier3d.model.BarrierLength):
            width = next((index for index, value in enumerate(self._barrier3d.model.InteriorDomain[:, bl]) if value <= self._barrier3d.model.SL), StartDomainWidth)
            width = width - 1
            if width < 0:
                width = 0
            InteriorWidth[bl] = width

        # Update Barrier3D Domain Sizes
        Target_width_barriermarsh = self._bmftc_BB.B - self._bmftc_BB.x_m - self._x_s_offset  # [m] Target width of barrier-marsh
        Target_width_barriermarsh = math.ceil(Target_width_barriermarsh / 10)  # [dam]

        addRows = Target_width_barriermarsh - StartDomainWidth + 1  # Number of rows to add (if positive) or subtract (if negative) from Barrier3D domain

        if addRows > 0:
            # Update interior domain size
            Marsh_Addition = np.ones([addRows, self._barrier3d.model.BarrierLength]) * -self._barrier3d.model._BayDepth
            Zero_Addition = np.zeros([addRows, self._barrier3d.model.BarrierLength])
            NewDomain = np.vstack([self._barrier3d.model.InteriorDomain, Marsh_Addition])
            # Update size of shrub domains, too
            self._barrier3d.model._ShrubDomainFemale = np.vstack([self._barrier3d.model.ShrubDomainFemale, Zero_Addition])
            self._barrier3d.model._ShrubDomainMale = np.vstack([self._barrier3d.model.ShrubDomainMale, Zero_Addition])
            self._barrier3d.model._ShrubDomainDead = np.vstack([self._barrier3d.model.ShrubDomainDead, Zero_Addition])
            self._barrier3d.model._ShrubPercentCover = np.vstack([self._barrier3d.model.ShrubPercentCover, Zero_Addition])
            self._barrier3d.model._DeadPercentCover = np.vstack([self._barrier3d.model.DeadPercentCover, Zero_Addition])
            self._barrier3d.model._BurialDomain = np.vstack([self._barrier3d.model.BurialDomain, Zero_Addition])
            self._barrier3d.model._ShrubDomainAll = self._barrier3d.model.ShrubDomainFemale + self._barrier3d.model.ShrubDomainMale
        elif addRows < 0:
            # Update interior domain size
            NewDomain = self._barrier3d.model.InteriorDomain[:addRows, :]
            # Update size of shrub domains, too
            self._barrier3d.model._ShrubDomainFemale = self._barrier3d.model.ShrubDomainFemale[:addRows, :]
            self._barrier3d.model._ShrubDomainMale = self._barrier3d.model.ShrubDomainMale[:addRows, :]
            self._barrier3d.model._ShrubDomainDead = self._barrier3d.model.ShrubDomainDead[:addRows, :]
            self._barrier3d.model._ShrubPercentCover = self._barrier3d.model.ShrubPercentCover[:addRows, :]
            self._barrier3d.model._DeadPercentCover = self._barrier3d.model.DeadPercentCover[:addRows, :]
            self._barrier3d.model._BurialDomain = self._barrier3d.model.BurialDomain[:addRows, :]
            self._barrier3d.model._ShrubDomainAll = self._barrier3d.model._ShrubDomainFemale + self._barrier3d.model._ShrubDomainMale
        else:
            NewDomain = self._barrier3d.model.InteriorDomain  # Domains stay same size

        if len(marsh_transect) >= 1:
            # Update Marsh In Barrier3D
            x_marsh = Target_width_barriermarsh + 1  # [dam] Cross-shore location of marsh edge relative to interior domain
            for w in range(self._barrier3d.model.BarrierLength):
                width_diff = x_marsh - (InteriorWidth[w] + len(marsh_transect))
                if width_diff < 0:
                    MarshTransect = marsh_transect[:-int(abs(width_diff))]  # [dam]
                elif width_diff > 0:
                    add = np.ones([int(abs(width_diff))]) * marsh_transect[-1]  # Set additional marsh cells to elevation of last marsh
                    MarshTransect = np.append(marsh_transect, add)  # [dam]
                else:
                    MarshTransect = marsh_transect  # [dam]

                InteriorTransect = NewDomain[:InteriorWidth[w], w]  # [dam]
                BarrierMarshTransect = np.append(InteriorTransect, MarshTransect)  # Combine interior and marsh

                NewDomain[:len(BarrierMarshTransect), w] = BarrierMarshTransect
                NewDomain[len(BarrierMarshTransect):, w] = (self._barrier3d.model.SL - np.mean([self._bmftc_ML.db, self._bmftc_BB.db])) / 10

        self._barrier3d.model.InteriorDomain = NewDomain

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Advance Barrier3D
        self._barrier3d.update()

        # Check if barrier has drowned
        if self._barrier3d.model.drown_break == 1:
            self._bmftc_ML._dur = time_step
            self._bmftc_ML._endyear = self._bmftc_ML.startyear + time_step
            self._bmftc_BB._dur = time_step
            self._bmftc_BB._endyear = self._bmftc_BB.startyear + time_step
            self._Barrier3D_Break = True
            print("Barrier3D Simulation Break: barrier drowned")
            return  # If so, end simulation

        # ===================================================================================================================================================================================================================================
        # ===================================================================================================================================================================================================================================
        # Update PyBMFT-C transect elevation based on Barrier3D elevation change

        shoreline_change = self._barrier3d.model.x_s_TS[-1] - self._barrier3d.model.x_s_TS[-2]
        self._x_s_offset = self._x_s_offset + (shoreline_change * 10)

        start_b3d = np.mean(NewDomain, axis=1) * 10  # Barrier3D domain before update, averaged across alongshore dimension, converted to m (vertical dimension)
        end_b3d = np.mean(self._barrier3d.model.InteriorDomain, axis=1) * 10  # Barrier3D domain after update, averaged across alongshore dimension, converted to m (vertical dimension)

        # Update start domain size to match end domain
        sc_b3d = self._barrier3d.model.ShorelineChangeTS[-1]  # Shoreline change [dam] from Barrier3D model update (this timestep)
        if sc_b3d < 0:  # Shoreline erosion
            start_b3d = start_b3d[abs(sc_b3d):]  # Trim off front
        elif sc_b3d > 0:  # Shoreline progradation
            add = np.zeros([sc_b3d])
            start_b3d = np.append(add, start_b3d)  # Add zeros to front

        if len(start_b3d) < len(end_b3d):
            add = np.ones([len(end_b3d) - len(start_b3d)]) * np.mean([self._bmftc_ML.db, self._bmftc_BB.db]) * -1  # Add bay cells
            start_b3d = np.append(start_b3d, add)
        elif len(start_b3d) > len(end_b3d):
            subtract = len(end_b3d) - len(start_b3d)
            start_b3d = start_b3d[:subtract]

        # Calculate change in elevation from Barrier3D update
        end_b3d = end_b3d + (self._barrier3d.model.RSLR[time_step] * 10)  # Offset sea-level rise from Barrier3D so that it isn't counted twice (i.e. RSLR already taken into account in PyBMFT-C)
        elevation_change_b3d = end_b3d - start_b3d  # Change in elevation across transect after Barrier3d update; [dam] horizontal dimension, [m] vertical dimentsion

        # Interpolate from dam to m (horizontal dimension)
        x = np.linspace(1, len(elevation_change_b3d) * 10, num=len(elevation_change_b3d) * 10)
        xp = np.linspace(1, len(elevation_change_b3d), num=len(elevation_change_b3d)) * 10
        xp = xp - 5
        elevation_change_b3d = np.interp(x, xp, elevation_change_b3d)

        off = int(abs(math.floor(self._x_s_offset)))  # [m] Offset of barrier shoreline and B

        # Incorporate elevation change from Barrier3D into back-barrier instance of PyBMFT-C
        if int(math.floor(self._x_s_offset)) < 0:
            elevation_change_b3d = np.flip(elevation_change_b3d[off:])  # Flip orientation
            marsh_barrier_width = (self._bmftc_BB.B - self._bmftc_BB.x_m)
            x_m_change = abs(math.floor(len(elevation_change_b3d) - marsh_barrier_width))  # Location of marsh edge within elevation_change_b3d

            # Add subaerial elevation change
            self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, -len(elevation_change_b3d[x_m_change:]):] += elevation_change_b3d[x_m_change:]
            # Store mass of overwash mineral sediment deposited across transect
            self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, -len(elevation_change_b3d[x_m_change:]):] += (elevation_change_b3d[x_m_change:] * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Determine volume of sed deposited past initial marsh edge and into bay
            sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
            sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
            self._bmftc_BB._Fow_min = max(0, sum_bay_dep * self._bmftc_BB.rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass
            self._bmftc_ML._Fow_min = max(0, sum_bay_dep * self._bmftc_ML.rhos)  # [kg/yr] Same as above for ML instance

            # Add volume of carryover from last time step
            sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

            # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
            # new_marsh_height = self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_m] - (self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.db)  # [m]
            new_marsh_height = self._bmftc_BB.db
            # new_marsh_height = self._bmftc_BB.db - self._bmftc_BB.amp

            # Determine distance of marsh progradation from overwash deposition
            progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
            progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
            self._bay_overwash_carryover = (progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

            if progradation > 0:
                # Add subaqueous elevation change
                self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += new_marsh_height
                # Store mass of overwash mineral sediment deposited across transect
                self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += (new_marsh_height * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Spread overwash bay flux evenly across bay bottom
            bay_accrete = sum_bay_dep / (self._bmftc_BB.bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
            self._bmftc_BB._db = self._bmftc_BB.db + bay_accrete  # Update bay depth
            self._bmftc_ML._db = self._bmftc_ML.db + bay_accrete  # Update bay depth

        elif int(math.floor(self._x_s_offset)) > 0:
            elevation_change_b3d = np.flip(elevation_change_b3d)
            marsh_barrier_width = (self._bmftc_BB.B - self._bmftc_BB.x_m)
            x_m_change = abs(math.floor(len(elevation_change_b3d) - (marsh_barrier_width - off)))  # Location of marsh edge within elevation_change_b3d

            # Add subaerial elevation change
            self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.B - off - len(elevation_change_b3d[x_m_change:]): self._bmftc_BB.B - off] += elevation_change_b3d[x_m_change:]
            # Store mass of overwash mineral sediment deposited across transect
            self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, self._bmftc_BB.B - off - len(elevation_change_b3d[x_m_change:]): self._bmftc_BB.B - off] += (elevation_change_b3d[x_m_change:] * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Determine volume of sed deposited past initial marsh edge and into bay
            sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
            sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
            self._bmftc_BB._Fow_min = max(0, sum_bay_dep * self._bmftc_BB.rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass
            self._bmftc_ML._Fow_min = max(0, sum_bay_dep * self._bmftc_ML.rhos)  # [kg/yr] Same as above for ML instance

            # Add volume of carryover from last time step
            sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

            # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
            # new_marsh_height = self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_m] - (self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.db)  # [m]
            new_marsh_height = self._bmftc_BB.db
            # new_marsh_height = self._bmftc_BB.db - self._bmftc_BB.amp

            # Determine distance of marsh progradation from overwash deposition
            progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
            progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
            self._bay_overwash_carryover = (progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

            if progradation > 0:
                # Add subaqueous elevation change
                self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += new_marsh_height
                # Store mass of overwash mineral sediment deposited across transect
                self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += (new_marsh_height * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Spread 50% of overwash bay flux evenly across bay bottom
            bay_accrete = sum_bay_dep / (self._bmftc_BB.bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
            self._bmftc_BB._db = self._bmftc_BB.db + bay_accrete  # Update bay depth
            self._bmftc_ML._db = self._bmftc_ML.db + bay_accrete  # Update bay depth

            # Remove barrier at front and set to msl to account for shoreline change
            self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, -off:] = self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp

        else:
            elevation_change_b3d = np.flip(elevation_change_b3d)
            marsh_barrier_width = (self._bmftc_BB.B - self._bmftc_BB.x_m)
            x_m_change = abs(math.floor(len(elevation_change_b3d) - marsh_barrier_width))  # Location of marsh edge within elevation_change_b3d

            # Add subaerial elevation change
            self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, -len(elevation_change_b3d[x_m_change:]):] += elevation_change_b3d[x_m_change:]
            # Store mass of overwash mineral sediment deposited across transect
            self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, -len(elevation_change_b3d[x_m_change:]):] += (elevation_change_b3d[x_m_change:] * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Determine volume of sed deposited past initial marsh edge and into bay
            sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
            sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
            self._bmftc_BB._Fow_min = max(0, sum_bay_dep * self._bmftc_BB.rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass
            self._bmftc_ML._Fow_min = max(0, sum_bay_dep * self._bmftc_ML.rhos)  # [kg/yr] Same as above for ML instance

            # Add volume of carryover from last time step
            sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

            # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
            # new_marsh_height = self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step - 1, self._bmftc_BB.x_m] - (self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.db)  # [m]
            new_marsh_height = self._bmftc_BB.db
            # new_marsh_height = self._bmftc_BB.db - self._bmftc_BB.amp

            # Determine distance of marsh progradation from overwash deposition
            progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
            progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
            self._bay_overwash_carryover = (progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

            if progradation > 0:
                # Add subaqueous elevation change
                self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += new_marsh_height
                # Store mass of overwash mineral sediment deposited across transect
                self._bmftc_BB.mineral_dep[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m - progradation: self._bmftc_BB.x_m] += (new_marsh_height * self._bmftc_BB.rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

            # Spread 50% of overwash bay flux evenly across bay bottom
            bay_accrete = sum_bay_dep / (self._bmftc_BB.bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
            self._bmftc_BB._db = self._bmftc_BB.db + bay_accrete  # Update bay depth
            self._bmftc_ML._db = self._bmftc_ML.db + bay_accrete  # Update bay depth

        # Calculate new marsh and "forest" edge positions after overwash
        self._bmftc_BB._x_m = self._bmftc_BB.x_m - progradation
        try:
            self._bmftc_BB._x_f = max(self._bmftc_BB.x_m + 1, np.where(self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, :] > self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.Dmin + 0.03)[0][0])
        except IndexError:
            self._bmftc_BB._x_f = self._bmftc_BB.B
            # If x_f can't be found, barrier has drowned
            self._bmftc_ML._dur = time_step
            self._bmftc_ML._endyear = self._bmftc_ML.startyear + time_step
            self._bmftc_BB._dur = time_step
            self._bmftc_BB._endyear = self._bmftc_BB.startyear + time_step
            self._BMFTC_Break = True
            print("PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full")
            return  # End simulation

        # Store new positions
        self._bmftc_BB.Marsh_edge[self._bmftc_BB.startyear + time_step] = self._bmftc_BB.x_m  # Save to array
        self._bmftc_BB.Forest_edge[self._bmftc_BB.startyear + time_step] = self._bmftc_BB.x_f  # Save to array

        # Determine change in x_b location
        self._bmftc_ML._x_b = self._bmftc_ML.x_b + progradation
        self._x_b_TS_ML[time_step] = self._bmftc_ML.x_b  # Save to array
        self._x_b_TS_BB[time_step] = self._bmftc_BB.x_b  # Save to array

        # Determine new fetch based on change in opposite marsh - both fetches should be exactly the same!
        self._bmftc_ML._bfo = self._bmftc_ML.bfo - progradation
        self._bmftc_BB._bfo = self._bmftc_BB.bfo - progradation
        self._bmftc_ML.fetch[self._bmftc_ML.startyear + time_step] = self._bmftc_ML.bfo  # Save to array
        self._bmftc_BB.fetch[self._bmftc_BB.startyear + time_step] = self._bmftc_BB.bfo  # Save to array

        # Update marsh scarp heigh parameter
        self._bmftc_BB._dmo = self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m]

        # Store landscape type widths for this time step
        if int(math.floor(self._x_s_offset)) < 0:
            barrier_width = len(self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_f:]) + off
        elif int(math.floor(self._x_s_offset)) > 0:
            barrier_width = len(self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_f:]) - off
        else:
            barrier_width = len(self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_f:])
        BB_marsh_width = (self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m: self._bmftc_BB.x_f] > self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.Dmax).sum()
        ML_marsh_width = (self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step, self._bmftc_ML.x_m: self._bmftc_ML.x_f] > self._bmftc_ML.msl[self._bmftc_ML.startyear + time_step] + self._bmftc_ML.amp - self._bmftc_ML.Dmax).sum()
        BB_marsh_pond_width = (self._bmftc_BB.elevation[self._bmftc_BB.startyear + time_step, self._bmftc_BB.x_m: self._bmftc_BB.x_f] < self._bmftc_BB.msl[self._bmftc_BB.startyear + time_step] + self._bmftc_BB.amp - self._bmftc_BB.Dmax).sum()
        ML_marsh_pond_width = (self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step, self._bmftc_ML.x_m: self._bmftc_ML.x_f] < self._bmftc_ML.msl[self._bmftc_ML.startyear + time_step] + self._bmftc_ML.amp - self._bmftc_ML.Dmax).sum()
        forest_width = len(self._bmftc_ML.elevation[self._bmftc_ML.startyear + time_step, self._bmftc_ML.x_f:])
        self._LandscapeTypeWidth_TS[time_step, :] = [barrier_width, BB_marsh_width, self._bmftc_BB.bfo, ML_marsh_width, forest_width, BB_marsh_pond_width, ML_marsh_pond_width]

        # ===================================================================================================================================================================================================================================

    @property
    def bmftc(self):
        return self._bmftc_ML

    @property
    def barrier3d(self):
        return self._barrier3d

    @property
    def bmftc_BB(self):
        return self._bmftc_BB

    @property
    def bmftc_ML(self):
        return self._bmftc_ML

    @property
    def cumul_len_change(self):
        return self._cumul_len_change

    @property
    def delta_fetch_BB_TS(self):
        return self._delta_fetch_BB_TS

    @property
    def delta_fetch_ML_TS(self):
        return self._delta_fetch_ML_TS

    @property
    def x_b_TS_ML(self):
        return self._x_b_TS_ML

    @property
    def BMFTC_Break(self):
        return self._BMFTC_Break

    @property
    def Barrier3D_Break(self):
        return self._Barrier3D_Break

    @property
    def LandscapeTypeWidth_TS(self):
        return self._LandscapeTypeWidth_TS

    @property
    def name(self):
        return self._name
