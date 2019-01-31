# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:29:28 2019

@author: ljander
"""
import numpy as np
from scipy.integrate import odeint

# Calculate the photosynthesis rate of a certain layer of the canopy
def layer_photo(eps, P_gm, DIS_1, DIS_2, DIS_3, WT_1, WT_2, WT_3, LAI, PAR,
                k_ext):
    # DIS_i is the distance coefficient of gauss integral (i is the layer number
    # of the canopy layers)
    # LAI is leaf area index
    # PAR is photosynthetical active radiation
    # k_ext is the extinction coefficient of canopy
    # P_gm is the maximum photosynthesis rate
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


# Simulate the dry matter production of tomato plant
def dry_matter_production(W_p0, LAI0, leaf0, root0, stem0, fruit0, t, t_start,
                          T_average, C_ppm, DIS_1, DIS_2, DIS_3, PAR, TEP,
                          k_ext, WT_1, WT_2, WT_3, C_f, f):
    # W_p0, LAI0, leaf0, root0, stem0, fruit0 are the intial value for each state variables
    # W_p is the dry weight of whole plant
    # LAI is leaf area index
    # leaf is the dry weight of leaf
    # root is the dry weight of root
    # stem is the dry weight of stem
    # fruit is the dry weight of fruit
    # t is the time step
    # t_start is the planting data (t_start_th day of the year)
    # C_ppm is the CO2 concentration
    # T_average is the average temperature inside greenhouse
    def diff(y, t_):
        # define state variables
        # TEP is the cumulative PAR
        # C_f is the conversion factor from assimilates to dry matter
        # f is a regression parameter
        W_p = y[0]  # [g DM]
        LAI = y[1]
        leaf = y[2]  # [g DM]
        root = y[3]  # [g DM]
        stem = y[4]  # [g DM]
        fruit = y[5]  # [g DM]
        # calculate the eps and P_gm based on temperature and co2 concentration
        C_a = C_ppm
        ghe = 42.7 + 1.68 * (T_average - 25) + 0.012 * (T_average - 25)**2
        eps = 1.544 * 10**-5 * (C_a - ghe) / (C_a + 2 * ghe)
        P_gm = 0.013 * (C_a - ghe)
        # calculate the plant photosynthetic rate [g CO2 m**-2 h**-1]
        P_g = layer_photo(eps, P_gm, DIS_1, DIS_2, DIS_3, WT_1, WT_2, WT_3,
                          LAI, PAR, k_ext)

        # simulate growth rate (dry matter production rate) of plant [g DM h**-1]
        # maintenance rate R_m is calculated based on the four type of organs of tomato plant
        R_m = (0.03 * leaf + 0.015 * stem + 0.01 * root +
               0.01 * fruit) / 24 * 2**((T_average - 25) / 10)

        # CO2 released by growth respiration is calculated as a factor times dWdt, for implementation
        # , the equation is divided by 1.191
        dWdt = C_f * (P_g * 30 / 44 - R_m *
                      (1 - np.exp(-f *
                                  (-0.0008 * np.log(W_p) + 0.0055)))) / 1.191
        # RGR is computed by equation -0.0008*np.log(W_p)+0.0055 obtained using linear regression analysis (Appendix 2)
        # 30/44 converts from CO2 to assimilates [g CH2O g**-1 CO2]

        # tomato dry matter partition (first partition to root and shoot)
        # PIS is the partition index of shoot
        # PIR is the partition index of root
        PIS = 1 - 0.12 * np.exp(-TEP / 100)
        PIR = 1 - PIS
        # simulate root growth
        d_rootdt = dWdt * PIR

        # partitioning of leaf, stem and fruit is based on the dry matter partitioned to the shoot
        # PIL is the partition index of leaf
        # PIST is the partition index of stem
        # PIF is the partition index of fruit
        PIL = PIS * (0.23 + 0.59 * np.exp(-TEP / 110))

        # simulate leaf growth rate
        d_leafdt = dWdt * PIL

        PIST = np.piecewise(TEP, [0 <= TEP <= 21, 21 < TEP],
                            [0.02 * TEP, 0.2 + 0.3 * np.exp(-TEP / 108)]) * PIS
        # simulate stem growth rate
        d_stemdt = dWdt * PIST

        PIF = PIS * (1 - PIL - PIST)
        PIF_p = np.piecewise(PIF, [PIF < 0, 0 <= PIF], [0, PIF])

        # simulate fruit growth rate
        d_fruitdt = dWdt * PIF_p
        # simulate specific leaf area based on t_start
        SLA = (266 + 88 * np.sin(2 * np.pi *
                                 (t_ / 24 + 68 + t_start) / 365)) * 1 * 10**-4
        # simulate leaf area index using leaf dry matter production
        d_laidt = np.piecewise(LAI, [0 <= LAI <= 3, LAI > 3],
                               [SLA * d_leafdt, 0])
        val = [dWdt, d_laidt, d_leafdt, d_rootdt, d_stemdt, d_fruitdt]

        return val
    # set initial state for the differential equation
    y0 = np.asarray([W_p0, LAI0, leaf0, root0, stem0, fruit0])
    # solve the differential equation using ‘odeint’ (integrate a system of ordinary differential equations)
    rez = odeint(diff, y0, t, full_output=0)
    return rez


W_p0 = ???
LAIO = ???
T_average = 28.0
dmp = dry_matter_production(W_p0, LAI0, leaf0, root0, stem0, fruit0, t, t_start,
                            T_average, C_ppm, DIS_1, DIS_2, DIS_3, PAR, TEP,
                            k_ext, WT_1, WT_2, WT_3, C_f, f)

print("Dry matter production is: {}".format(dmp))