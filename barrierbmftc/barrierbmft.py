"""
BarrierBMFT: Coupled Barrier-Bay-Marsh-Forest Model

Couples Barrier3D (Reeves et al., 2021) with the BMFT-C model (python version)

Copyright Ian RB Reeves
Last updated: 6 July 2021
"""

import numpy as np

from barrier3d import Barrier3d
from bmftc import Bmftc

class BarrierBMFT:
    """

    """
    def __init__(
            self,
            name,
            relative_sea_level_rise=1,
            reference_concentration=10,

    ):