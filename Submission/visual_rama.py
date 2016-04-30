import getcoordinates
import rama
from visual import *
import numpy as np

# based on getcoordinates
tmax = 2 #units of years
dt = 0.001 #units of years
t = np.linspace(0, tmax, int(tmax/dt))

r, v = getcoordinates.get_perifocal()
er, ev = getcoordinates.get_heliocentric(r, v)

rama_r =[[-2.0002,1.99998,.19994]]
rama_v = [[7.97836,-3.667243,-0.005509]]

er, ev = getcoordinates.add_Rama(er, ev, rama_r, rama_v)

sun = sphere(pos = er[0], radius = 0.15, color = color.yellow, make_trail=True, material=materials.emissive)
mercury = sphere(pos = er[1], radius = 0.02, color = (1,0.7,0.2), make_trail=True, material=materials.rough)
venus = sphere(pos = er[2], radius =0.05, color = color.yellow, make_trail=True, material=materials.marble)
earth = sphere(pos = er[3], radius =0.05, make_trail=True, material=materials.BlueMarble)
mars = sphere(pos = er[4], radius =0.03, color = color.red, make_trail=True, material=materials.bricks)
jupiter = sphere(pos = er[5], radius =0.1, make_trail=True, color = (.7,.3,0.1), material=materials.marble)
rama_sphere = cylinder(pos = er[10], axis=ev[10]/20, radius=0.01, make_trail=True, color = color.cyan)
#saturn = sphere(pos = er[6], radius =0.1, make_trail=True)
#uranus = sphere(pos = er[7], radius =0.1, make_trail=True)
#neptune = sphere(pos = er[8], radius =0.1, make_trail=True)
scene = display(center=rama_sphere.pos, forward=rama_sphere.pos)


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
    rama_sphere.axis = .1*vt[i,10]/np.linalg.norm(vt[i,10])
    rate(70)
    scene.center=rama_sphere.pos
