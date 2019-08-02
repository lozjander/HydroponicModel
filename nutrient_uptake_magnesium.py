# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 17:16:56 2019

@author: ljander
"""

# simulate magnesium uptake by tomato plant
def mg_uptake(t, RSA, PAR, EC,a,b,c,d,DSR):
    """Calculates Mg uptake.

    Parameters
    ----------
    t : int
      hour of the simulation
    RSA : float
      root surface area of tomato plant [m^2]

    Returns
    -------
    uptake: float
      magnesium uptake in summit or other


    Notes
    -----
    # simulate magnesium uptake using Michaelis-Menten active uptake model
    # t is the hour of the simulation
    # RSA is the root surface area of tomato plant [m^2]
    # PAR is photosynthetically active radiation [J m^-2 h^-1]
    # EC is the electrical conductivity of nutrient solution
    # U_max is the maximum rate of magnesium uptake [mg m^-2 h]
    # K_m is the Michaelis-Menten constant (the photosynthetic photon flux density at 1/2 J_max)
    # DSR is daily solar radiation [MJ m^-2]
    # a,b,c,d are regression parameters for calculating U_max and K_m
    # assuming the maximum and minimum temperature of greenhouse equals 26 and 18 degrees
    # the growing degree days is simulated using t

    """


    T_max = 26
    T_min = 18
    T_base = 8
    GDD = ((T_max+T_min)/2-T_base)*t/24

    # simulate photo thermal unit
    PTU = GDD*DSR
    # compute Michaelis-Menten parameters
    U_max = a*EC*PTU+b
    K_m = c*EC*PTU+d
    # converting [umol m^-2 s] to [J m^-2 h^-1] assuming an average wavelength of 550 nm
    wl = 550*1*10**-9 # [m]
    c = 3*10**8 # [m/s]
    h = 6.626*10**-34
    NA = 6.022*10**23
    # unit conversion
    K_m = c/wl*h*NA*1*10**-6*3600*K_m
    # calculate uptake based on Michaelis-Menten type model
    uptake = RSA*PAR*U_max/(K_m+PAR)
    return uptake
