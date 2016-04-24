# Equations for getting to state vector of orbital elements
# initial calculations from page 389 of Orbital Mechanics for Engineering Students, by Howard Curtis
# 
#
#
import numpy as np
import rama.py

#gravitational parameter of the Sun (p.388 of OMES)
mu = 1.327e11  #in km^3 / s^2

#Angular Momentum
def ang_mom(mu, a, e)
""" Define angular momentum for state vector calculation.

Inputs
------

mu : gravitational parameter

a  : semi-major axis

e  : eccentricity

Outputs
-------

h  : scalar, angular momentum of body
"""
h = np.sqrt(mu*a*(1-e**2))
return h

# True Anomaly (p.112 OMES, eqn 3.10)
def t_anom(E, e)
""" Calculates the true anomaly of a body

angular parameter that defines the position of a body moving along a Keplerian orbit. 
It is the angle between the direction of periapsis (closest approach to cent. body) 
and the current position of the body, as seen from the main focus of the ellipse 
(the point around which the object orbits).

Inputs
------

E : eccentric anomaly

e : eccentricity of orbiting body

Outputs
-------

theta : angle
"""
theta = 2*np.arctan(np.tan(E/2)*np.sqrt((1-e)/(1+e)))
return theta

#Rotational matrix about the z-axis through the angle w
def z_rot(angle)
""" evaluate the z axis Euler rotation matrix for omega (lowercase) the argument of perihelion and Omega (uppercase) the right ascenscion of ascending node.
Taken from eqns. 4.39 and 4.41 p. 173

** In using make sure the omegas are not confused at output.

Inputs
------

omega : angle

Outputs
-------

r_omega : rotation matrix used in calculation

from p.611 of OMES
 cos(omega)  sin(omega  0
-sin(omega)  cos(omega) 0
     0           0      1

"""
z_rotmat = np.zeros([3,3], dtype=np.float)
z_rotmat[0,0] = z_rotmat[1,1] = np.cos(angle)
z_rotmat[0,1] = np.sin(angle)
z_rotmat[1,0] = -np.sin(angle)
z_rotmat[0,2] = z_rotmat[1,2] = z_rotmat[2,0] = z_rotmat[2,1] = 0
z_rotmat[2,2] = 1

return z_rotmat

#Rotational matrix about the x-axis through the angle i (inclination)
def x_rot(angle)
""" evaluate the x axis Euler rotation matrix for i, inclination.

taken from eqn. 4.40 p. 173

Inputs
------

i : angle

Outputs
-------

r_incl : rotation matrix used in calculation

from p.611 of OMES
  1      0            0
  0   cos(incl)    sin(incl)
  0   -sin(incl)   cos(incl) 
"""
x_rotmat = np.zeros([3,3], dtype=np.float)
x_rotmat[0,0] = 1
x_rotmat[0,1] = x_rotmat[0,2] = x_rotmat[1,0] = x_rotmat[2,0] = 0
x_rotmat[1,1] = x_rotmat[2,2] = np.cos(angle)
x_rotmat[1,2] = np.sin(angle)
x_rotmat[2,1] = -np.sin(angle)

return x_rotmat


def state_from_coe() 
""" Calculate the state vector for planets under consideration
sourced from algorithm 4.2 in OMES p. 610

Outputs
-------

r : position vector

v : velocity vector
"""
# Angular Momentum
h = ang_mom(mu, system["Semi-Major Axis"], system["Eccentricity"])
# True Anomaly
theta = t_anom(E, system["Eccentricity"])

# 4.37 position array
r_ar = np.zeros([3,1], dtype=np.float)
r_ar[0] = np.cos(theta)
r_ar[1] = np.sin(theta)
r_ar[2] = 0

#4.38 velocity array
v_ar = np.zeros([3,1], dtype=np.float)
v_ar[0] = -np.sin(theta)
v_ar[1] = e + np.cos(theta)
v_ar[2] = 0

# equations 4.37 and 4.38 in OMES p. 173
rp = (h**2/mu) * (1/(1 + system["Eccentricity"]*np.cos(theta))) * r_ar
vp = (mu/h) * v_ar

# p. 173 eqn. 4.42, transition of planes
Q = np.zeros([3,3], dtype=np.float64)
Q = z_rot(omega) * x_rot(system["Inclination"]) * z_rot(system["Ascending Longitude"])

#eqn. 4.46, column vectors
r = Q*rp
v = Q*vp

# Transpose to row vectors
rT = np.transpose(r)
vT = np.transpose(v)

return rT, vT









