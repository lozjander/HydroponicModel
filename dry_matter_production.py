import numpy as np
from scipy.integrate import odeint


def dry_matter_production(W_p0, LAI0, leaf0, root0, stem0, fruit0, t, t_start,
                          T_average, C_ppm, DIS_1, DIS_2, DIS_3, PAR, TEP,
                          k_ext, WT_1, WT_2, WT_3, C_f, f):
    """Simulates the dry matter production of tomato plant.

    Parameters
    ----------
    W_p0 : float
      Initial dry weight of the whole tomato plant [g DM]
    LAI_0 : float
      Initial leaf area index of the tomato plant [g DM]
    leaf0 : float
      Initial dry weight of the tomato plant leaves [g DM]
    root0 : float
      Initial dry weight of the tomato plant roots [g DM]
    stem0 : float
      Initial dry weight of the tomato plant stem [g DM]
    fruit0 : float
      Initial dry weight of the tomato plant fruits [g DM]
    t : ???
      timestep
    t_start: int
      is the planting data (t_start_th day of the year)
    T_average : float
      the average temperature inside greenhouse
    C_ppm : int
      CO2 concentration
    DIS_i :int
      the distance coefficient of gauss integral (i is the layer number of the canopy layers)
    PAR : int
      photosynthetically active radiation
    TEP : int
      the cumulative PAR
    k_ext: float
      the extinction coefficient of canopy
    WT_i : float
      Weight factor for Gaussian integration
    C_f : float
      Conversion factor from assimilates to dry matter
    f : float
     Regression parameter

    Returns
    -------
    Given a set of initial values (y0) calculated by diff , it returns an array containing the value of y for each
    desired time in t, with the initial value `y0` in the first row.

    Array of values for Leaf Area Index, Whole plant dry weight, & Leaf, stem, fruit, root dry weights.


    Notes
    -----
    # are all the parameters just lists???????
            """

    def diff(y, t):
        """Define state variables (y) at a timestep (t).

        Parameters
        ----------
        y : ???
          list of state variables
        t : ???
          timestep

        Variables
        ---------
        C_a : int
          CO2 concentraion [ppm]
        ghe : float
         'Gamma' the CO2 compensation point in absence of dark respiration
        eps :
          leaf initial light use efficiency [g.CO2.j^−1]
        P_gm : float
          the maximum photosynthesis rate
        R_m : float
          maintenance rate - of respiration?
        dWdt : float
          Whole Plant growth rate
        RGR : float
          Relative Growth Rate
        PIS : float
          Partition Index Shoot
        PIR : float
          Partition Index Root
        d_rootdt : float
          root growth rate
        PIL : float
          Partition Index Leaf
        PIST : float
          Partition Index Stem
        PIF : float
          Partition Index Fruit
        d_leafdt : float
          leaf growth rate
        d_stemdt : float
          stem growth rate
        d_fruitdt : float
          fruit growth rate
        SLA : float
          Specific Leaf Area



        Returns
        -------
        Val - a list of values (rate of change - ODE) calculated for timestep (t):
        Whole plant dry weight
        Leaf Area Index
        Leaf dry weight
        Root dry weight
        Stem dry weight
        fruit dry weight

        Notes
        -----
        # are all the parameters just lists???????

                """
        W_p = y[0]  # [g DM]
        LAI = y[1]  # LAI has no unit
        leaf = y[2]  # [g DM]
        root = y[3]  # [g DM]
        stem = y[4]  # [g DM]
        fruit = y[5]  # [g DM]
        # calculate the eps and P_gm based on temperature and co2 concentration
        C_a = C_ppm
        ghe = 42.7 + 1.68 * (T_average - 25) + 0.012 * (T_average - 25) ** 2
        eps = 1.544 * 10 ** -5 * (C_a - ghe) / (C_a + 2 * ghe)
        P_gm = 0.013 * (C_a - ghe)
        # calculate the plant photosynthetic rate [g CO2 m**-2 h**-1]
        P_g = layer_photosyn(eps, P_gm, DIS_1, DIS_2, DIS_3, WT_1, WT_2, WT_3,
                             LAI, PAR, k_ext)

        # simulate growth rate (dry matter production rate) of plant [g DM h**-1]
        # maintenance rate R_m is calculated based on the four type of organs of tomato plant
        R_m = (0.03 * leaf + 0.015 * stem + 0.01 * root +
               0.01 * fruit) / 24 * 2 ** ((T_average - 25) / 10)

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
        # simulate specific leaf area based on t_start - t_ could be 't' or 't_start' try both
        t_maybe = t
        SLA = (266 + 88 * np.sin(2 * np.pi *
                                 (t_maybe / 24 + 68 + t_start) / 365)) * 1 * 10 ** -4
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
