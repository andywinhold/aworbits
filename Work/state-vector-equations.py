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
def rot_omega(omega)
""" evaluate the rotation matrix for omega the argument of perihelion

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
r_omega = np.zeros([3,3], dtype=np.float)
r_omega[0,0] = r_omega[1,1] = np.cos(omega)
r_omega[0,1] = np.sin(omega)
r_omega[1,0] = -np.sin(omega)
r_omega[0,2] = r_omega[1,2] = r_omega[2,0] = r_omega[2,1] = 0
r_omega[2,2] = 1

return r_omega

#Rotational matrix about the x-axis through the angle i (inclination)
def rot_incl(i)
""" evaluate the rotation matrix for i, inclination

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
r_incl = np.zeros([3,3], dtype=np.float)
r_incl[0,0] = 1
r_incl[0,1] = r_incl[0,2] = r_incl[1,0] = r_incl[2,0] = 0
r_incl[1,1] = r_incl[2,2] = np.cos(i)
r_incl[1,2] = np.sin(i)
r_incl[2,1] = -np.sin(i)

return r_incl

#Rotational matrix about the z-axis through the angle Ohm (uppercase omega, longitude of Ascending Node)
def rot_along(ohm)
""" evaluate the rotation matrix for Ohm

Inputs
------

ohm : angle

Outputs
-------

r_ohm : rotation matrix used in calculation

from p.611 of OMES
  1      0            0
  0   cos(incl)    sin(incl)
  0   -sin(incl)   cos(incl) 
"""
r_ohm = np.zeros([3,3], dtype=np.float)
r_ohm[0,0] = r_ohm[1,1] = np.cos(ohm)
r_ohm[0,1] = np.sin(ohm)
r_ohm[1,0] = -np.sin(ohm)
r_ohm[0,2] = r_ohm[1,2] = r_ohm[2,0] = r_ohm[2,1] = 0
r_ohm[2,2] = 1

return r_incl


def state_from_class() 
""" Calculate the state vector for planets under consideration
sourced from algorithm 4.2 in OMES p. 610

Outputs
-------

r : position vector

v : velocity vector
"""
h = ang_mom(mu, a, e)
theta = t_anom(E, e)

# 4.37 array
r_ar = np.zeros([3,1], dtype=np.float)
r_ar[0] = np.cos(theta)
r_ar[1] = np.sin(theta)
r_ar[2] = 0
#4.38 array
v_ar = np.zeros([3,1], dtype=np.float)
v_ar[0] = -np.sin(theta)
v_ar[1] = e + np.cos(theta)
v_ar[2] = 0

# equations 4.37 and 4.38 in OMES p. 173
# components of the state vector of a body relative to its perifocal reference
# The perifocal coordinate (PQW) system is a frame of reference for an orbit. The frame is centered at the focus of the orbit, i.e. the celestial body about which the orbit is centered.
rp = (h**2/mu) * (1/(1 + e*np.cos(theta))) * r_ar
vp = (mu/h) * v_ar













