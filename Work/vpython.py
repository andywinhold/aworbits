import getcoordinates
import rama
from visual import *
import numpy as np

tmax = 1 #units of years
dt = 0.001 #units of years
t = np.linspace(0, tmax, int(tmax/dt))

r, v = getcoordinates.get_perifocal()
er, ev = getcoordinates.get_heliocentric(r, v)

er = np.append(er, [[4,0,1]], axis=0)
ev = np.append(ev, [[-10,1.1,-1.8]], axis = 0)


sun = sphere(pos = er[0], radius = 0.15, color = color.yellow, make_trail=True)
mercury = sphere(pos = er[1], radius = 0.1, make_trail=True)
venus = sphere(pos = er[2], radius =0.1, make_trail=True)
earth = sphere(pos = er[3], radius =0.1, make_trail=True)
mars = sphere(pos = er[4], radius =0.1, make_trail=True)
jupiter = sphere(pos = er[5], radius =0.1, make_trail=True)
rama_sphere = sphere(pos = er[10], radius=0.05, make_trail=True)
#saturn = sphere(pos = er[6], radius =0.1, make_trail=True)
#uranus = sphere(pos = er[7], radius =0.1, make_trail=True)
#neptune = sphere(pos = er[8], radius =0.1, make_trail=True)


xt, vt, KE, UE, TE = rama.dynamics(er, ev, dt, tmax)

for i in range(len(xt)):
    sun.pos = xt[i,0]
    mercury.pos = xt[i,1]
    venus.pos = xt[i,2]
    earth.pos = xt[i,3]
    mars.pos = xt[i,4]
    jupiter.pos = xt[i,5]
    #saturn.pos = xt[i,6]
    #uranus.pos = xt[i,7]
    #neptune.pos = xt[i,8]
    rama_sphere.pos = xt[i,10]

    rate(10)
