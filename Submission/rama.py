#Copyright (C) 2016  Jordan Boyd, Ashley Mascareno, Andrew Winhold
# Released under the GNU Public Licence, v3 or any higher version
# Please cite where necessary.

import numpy as np


#============================================================
# Parameters
#------------------------------------------------------------
G_gravity = 3.964016 * 9.945193 #in astronomical units^3/year^2, close to 4(pi^2)

#Masses in solar masses
mass = np.array([
            1.,                 #Sun
            0.330104/(1.989e6), #Mercury
            4.86732/(1.989e6),  #Venus
            5.97219/(1.989e6),  #Earth
            0.641693/(1.989e6), #Mars
            1898.19/(1.989e6),  #Jupiter
            568.319/(1.989e6),  #Saturn
            86.8103/(1.989e6),  #Uranus
            102.410/(1.989e6),  #Neptune
            0.0166/(1.989e6),   #Eris
            0.0166/(1.989e8)    #Rama, 1/100 mass of Eris
            ])

def gravity(coordinates):
    """Calculates total force on each solar system object using Newtonian gravity.
    
    Parameters
    ----------
    coordinates : array
        [x,y,z] coordinates for every solar system object
        
    Returns
    -------
    force : array
        Returns net force on each solar system object as (x,y,z) array and total potential as scalar.
    """
    
    N = len(coordinates) #number of objects
    r = np.zeros((N, 3)) #distance between objects
    force = np.zeros((N, 3)) #force between objects
    U = 0 #potential energy
    
    for i in range(0, N - 1):
        for j in range(i + 1, N):
            #calculates distance vector between each object
            r[i] = coordinates[j] - coordinates[i]
            
            #force calculation
            rmag = np.sqrt(r[i, 0]**2 + r[i, 1]**2 + r[i, 2]**2)
            rhat = r[i] / rmag
            dF =  rhat * G_gravity * mass[i] * mass[j] / (rmag**2)
            force[j] += -dF
            force[i] += dF
            U += -G_gravity * mass[i] * mass[j] / rmag    
            
    return force, U

def dynamics(x0, v0, dt=0.01, tmax=10):
    """Integrate equations of motions
    Parameters
    ----------
    x0 : array
         Nx3 array containing the starting coordinates of the planetary objects.
         (Note that x0 is changed and at the end of the run will
         contain the last coordinates.)
    v0 : array
         Nx3 array containing the starting velocities, eg np.zeros((N,3))
     
    masses : array
        masses for every object in solar system

    dt : float
         integration timestep in units?
         
    tmax : int, optional
         total run time
         
    Returns array of coordinates of planetary objects, their velocities, total kinetic, potential, and overall energies.
    """
    
    N = len(x0) #number of objects
    nsteps = int(tmax/dt)
    #x = np.zeros((nsteps/10,N,3))
    x = np.zeros((nsteps,N,3))
    vf = np.zeros((nsteps,N,3))
    dx = np.copy(x0)
    v = np.copy(v0)
    vhalf = np.zeros((N,3))
    Ut = np.zeros(nsteps)
    kinetic = np.zeros(nsteps)
    totalE = np.zeros(nsteps)
    
    Ft, Ut[0] = gravity(dx)
    
    for i in range(nsteps):
        for j in range(N):
            vhalf[j] = v[j] + 0.5 * dt * Ft[j] / mass[j]
            dx[j] += dt * vhalf[j]
        Ft, Ut[i]= gravity(dx)
        for j in range(N):
            v[j] = vhalf[j] + 0.5 * dt * Ft[j] / mass[j]
            kinetic[i] += 0.5 * mass[j] * np.sum(v[j]**2) 
        #if i%10 == 0:
            #x[int(i/10)] = dx
        x[i] = dx
        vf[i] = v
            
    totalE = kinetic + Ut        
    #If rama is present, get position of earth and rama and determine distance between the two.
    #----------------------------------------------------------------------
    if len(x[0,:]) == 11:
        earth_pos = np.zeros(len(x[:]))
        rama_pos = np.zeros_like(earth_pos)
        dist = np.zeros_like(earth_pos)
        dist_mag = np.zeros(len(earth_pos))

        earth_pos = x[:,3]
        rama_pos = x[:,10]
        #distance between the two
        dist = np.abs(earth_pos - rama_pos)
        #array to store the closer values
        close = np.zeros((nsteps,), dtype=np.float64)
        #if the distance is less than 0.01 AU then it prints the distance between Earth and Rama.
        for i in range(len(dist)):
            dist_mag[i] = np.linalg.norm(dist[i]) - 4.25875e-5 #The second value is the radius of Earth.
            if dist_mag[i] < .01:
                print("Iteration:",i,",",
                      "Rama distance from Earth (au):", dist_mag[i])


    return x, vf, kinetic, Ut, totalE
