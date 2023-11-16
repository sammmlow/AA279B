# -*- coding: utf-8 -*-
"""
AA 279B Final Project: Generation of Porkchop Plots to determine the total
time of flight and back-dating the time of departure using a fixed time of
arrival during Saturn-Titan summer solstice.
"""

# First step, porkchop plot for Earth to Venus in the time range 2040 - 2050.
# Determine the most feasible step to get to Venus for a gravity assist.
import matplotlib.pyplot as plt
from astropy import units as u

from poliastro.bodies import Earth, Mars, Saturn, Jupiter
from poliastro.plotting.porkchop import PorkchopPlotter
from poliastro.util import time_range

plt.close('all')

# Original Cassini mission:
# 1. Earth (15 Oct 1997) to Venus (26 Apr 1998)
# 2. Venus (26 Apr 1998) stay in current ellipse around Sun
# 3. Current Ellipse to a thrust (3 Dec 1998)
# 4. Thrust (3 Dec 1998) to Venus again (24 Jun 1999)
# 5. Earth (15 Oct 1997) to Venus (26 Apr 1998)
# 6. Earth (15 Oct 1997) to Venus (26 Apr 1998)

# Interplanetary Trajectory:
# Create porkchop plots for each of the following transfers:

# 1. Earth - Saturn
# 2. Earth - Mars - Saturn
# 3. Earth - Mars - Jupiter - Saturn
# 4. Earth - Jupiter - Saturn

# Arrival date at Saturn should be: Around 2048 because summer solstice
# is Jan 2051 and we need time to transfer to Titan. Then we back out
# the departure dates from there to minimize Delta-V

# Total time of flight takes reference from Cassini - 6.7 years, so we'll
# use a window of 6 to 8 years to get to Saturn, and similarly back to Earth.

# Earth to Saturn
c3_e2s = 250 * u.km**2 / u.s**2
depart_e2s = time_range("2037-09-15", end="2038-01-01")
arrive_e2s = time_range("2039-07-01", end="2045-01-01")
porkchop_e2s = PorkchopPlotter(Earth, Saturn, depart_e2s, arrive_e2s,
                                max_c3 = c3_e2s, vhp = False, tfl = False)
porkchop_e2s.porkchop()

# Earth to Mars
c3_e2m = 50 * u.km**2 / u.s**2
depart_e2m = time_range("2032-11-01", end="2033-08-01")
arrive_e2m = time_range("2033-07-01", end="2034-09-01")
porkchop_e2m = PorkchopPlotter(Earth, Mars, depart_e2m, arrive_e2m,
                                max_c3 = c3_e2m, vhp = False, tfl = False)
porkchop_e2m.porkchop()

# Earth to Jupiter
c3_e2j = 250 * u.km**2 / u.s**2
depart_e2j = time_range("2035-05-01", end="2035-10-01")
arrive_e2j = time_range("2037-03-01", end="2038-03-01")
porkchop_e2j = PorkchopPlotter(Earth, Jupiter, depart_e2j, arrive_e2j,
                                max_c3 = c3_e2j, vhp = False, tfl = False)
porkchop_e2j.porkchop()

# Mars to Saturn
c3_m2s = 150 * u.km**2 / u.s**2
depart_m2s = time_range("2035-08-01", end="2036-04-01")
arrive_m2s = time_range("2038-01-01", end="2045-01-01")
porkchop_m2s = PorkchopPlotter(Mars, Saturn, depart_m2s, arrive_m2s,
                                max_c3 = c3_m2s, vhp = False, tfl = False)
porkchop_m2s.porkchop()

# Mars to Jupiter
c3_m2j = 100 * u.km**2 / u.s**2
depart_m2j = time_range("2035-03-01", end="2035-11-01")
arrive_m2j = time_range("2037-01-01", end="2039-01-01")
porkchop_m2j = PorkchopPlotter(Mars, Jupiter, depart_m2j, arrive_m2j,
                                max_c3 = c3_m2j, vhp = False, tfl = False)
porkchop_m2j.porkchop()

# Jupiter to Saturn
c3_j2s = 30 * u.km**2 / u.s**2
depart_j2s = time_range("2033-01-01", end="2041-01-01")
arrive_j2s = time_range("2044-01-01", end="2047-01-01")
porkchop_j2s = PorkchopPlotter(Jupiter, Saturn, depart_j2s, arrive_j2s,
                                max_c3 = c3_j2s, vhp = False, tfl = False)
porkchop_j2s.porkchop()