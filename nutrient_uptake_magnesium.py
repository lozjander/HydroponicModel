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
    PAR : int
      photosynthetically active radiation [J m^-2 h^-1]
    EC : float
      electrical conductivity of the nutrient solution
    a : float
      regression parameter for calulcating U_max and K_m
    b : float
      regression parameter for calulcating U_max and K_m
    c : float
      regression parameter for calulcating U_max and K_m
    d : float
      regression parameter for calulcating U_max and K_m
    DSR : int
      daily solar radiation [MJ m^-2]

    Variables
    ---------
    U_max : float
      maximum rate of magnesium uptake [mg m^-2 h]
    K_m : float
      michaelis-menten constant [PPFD @ 0.5 J_max]
    PTU : float
      photo thermal unit
    GDD : float
      Growing Degree Days

    Constants
    ---------
    NA : float
      avogadro's constant [6.022*10**23]
    h : float
      planck constant [6.626*10**-34]
    c : float
      speed of light in a vacuum [3*10**8]


    Returns
    -------
    uptake: float
      magnesium uptake in summit or other


    Notes
    -----
    # assuming the maximum and minimum temperature of greenhouse equals 26 and 18 degrees
    # the growing degree days is simulated using t
    #  K_m is used twice, its calculated at the beginning for each time step then calculated to convert the unit it
    is expressed in from PPFD [umol m^-2 s] to [J m^-2 h^-1] Radiant Exposure per hour
    # are all the parameters just lists???????

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
    # unit conversion --- 3600 seconds in an hour, 1e-6 as PPFD is in micromol?
    K_m = c/wl*h*NA*1*10**-6*3600*K_m
    # calculate uptake based on Michaelis-Menten type model
    uptake = RSA*PAR*U_max/(K_m+PAR)
    return uptake
