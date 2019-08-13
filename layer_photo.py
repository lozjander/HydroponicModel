# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:29:28 2019

@author: ljander
"""
import numpy as np


def layer_photosyn(eps, P_gm, DIS_1, DIS_2, DIS_3, WT_1, WT_2, WT_3, LAI, PAR,
                   k_ext):
    """Calculates the photosynthesis rate of a certain layer of the canopy.

    Parameters
    ----------
    eps :
      leaf initial light use efficiency [g.CO2.j^−1]
    P_gm :
      the maximum photosynthesis rate
    DIS_i :int
      the distance coefficient of gauss integral (i is the layer number of the canopy layers)
    WT_i : float
      Weight factor for Gaussian integration
    LAI : float
      leaf area index - one-sided green leaf area per unit ground surface area
    PAR :
      photosynthetically active radiation [J m^-2 h^-1]
    k_ext:
      the extinction coefficient of canopy

    Variables
    ---------
    LGUSS_i : float
      canopy depth of gauss layer
    L_i : float
      the quantity of PAR arriving at the ith layer of the canopy
    P_gi : float
      the photosynthesis rate at the ith layer of the canopy

    Returns
    -------
    P_g: float
      photosythesis rate in the plant

    Notes
    -----
    # gauss integral - its the integral of a gaussian bell curve  e^(-x^2)
    # WT = 0.2778, 0.4444, 0.2778 for 1st, 2nd, 3rd canopy layer respectively
    # DIS = 0.1127, 0.5, 0.8873 for 1st, 2nd, 3rd canopy layer respectively
    # k_ext = 0.8
    # are all the parameters just lists???????
            """

    # is this section irrelevant due to using lettuce, there will be some shading??? - loz
    LGUSS_1 = DIS_1 * LAI
    LGUSS_2 = DIS_2 * LAI
    LGUSS_3 = DIS_3 * LAI
    # LGUSS_i is the canopy depth of gauss layer

    L_1 = PAR * k_ext * np.exp(-k_ext * LGUSS_1)
    L_2 = PAR * k_ext * np.exp(-k_ext * LGUSS_2)
    L_3 = PAR * k_ext * np.exp(-k_ext * LGUSS_3)
    # L_i is the quantity of PAR arriving at the ith layer of the canopy

    P_g1 = P_gm * (1 - np.exp(-eps * L_1 / P_gm))
    P_g2 = P_gm * (1 - np.exp(-eps * L_2 / P_gm))
    P_g3 = P_gm * (1 - np.exp(-eps * L_3 / P_gm))

    P_g = (P_g1 * WT_1 + P_g2 * WT_2 + P_g3 * WT_3) * LAI
    return P_g


