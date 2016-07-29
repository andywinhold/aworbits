# This document contains functions for calculating different times and dates. Specifically Julian day number and J2000.
# Written in code by Andrew Winhold, equation credit to Orbital Mechanics for Engineering Students by Howard Curtis

import numpy as np


# To calculate the Julian Day
def JD(J_0, UT):
    JD = J_0 + UT/24
    return JD


# To calculate initial Julian reference day.
def J_0(y, mo, d):
    """
    Inputs
    ------

    y : year, integer in range 1901 <= y <= 2099

    mo : month, integer in range 1 <= m <= 12

    d : day, integer in range 1 <= d <= 31

    Outputs
    -------

    J_0 : integer, initial julian date based on provided y, m, and d.
    """

    J_0 =  367 * y - np.int((7/4) * (y + np.int((mo + 9)/12))) + np.int((275 * mo )/9) + d + 1721013.5
    J_0 = round(J_0, 3)
    return J_0

# Calculate UT
def UT(hr, m, sec):
    UT = hr + (m/60) + (sec/3600)
    UT = round(UT, 3)
    return UT