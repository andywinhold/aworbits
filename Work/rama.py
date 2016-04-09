import numpy as np


#============================================================
# Parameters of the problem
#------------------------------------------------------------
#

G_gravity = 4*np.pi**2

# mass in AU
mass = {'Sun': 1.,
        'Uranus': 4.366244e-5,
        'Neptune': 5.151389e-5,
}

# orbital period in Earth years
period = {'Uranus': 84.0110,
          'Neptune': 164.7901,
}

def gravity(coordinates, masses):
    """Calculates total force on each solar system object using Newtonian gravity.
    
    Parameters
    ----------
    coordinates : array
        [x,y,z] coordinates for every solar system object
        
    masses : array
        masses for every object in solar system
        
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
            dF =  rhat * G_gravity * masses[i] * masses[j] / (rmag**2)
            force[j] += dF
            force[i] += -dF
            U += G_gravity * masses[i] * masses[j] / rmag    
            
    return force, U

def dynamics(x0, v0, dt, masses, nsteps=10):
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
         
    nsteps : int, optional
         number of integrator time steps
         
    Returns array of coordinates of planetary objects, their velocities, total kinetic, potential, and overall energies.
    """
    
    N = len(x0) #number of objects
    x = np.copy(x0)
    v = np.copy(v0)
    vhalf = np.zeros((N,3))
    Ut = np.zeros(nsteps)
    kinetic = np.zeros(nsteps)
    totalE = np.zeros(nsteps)
    
    Ft, Ut[0]= gravity(x, masses)
    
    for i in range(nsteps):
        for j in range(N):
            vhalf[j] = v[j] + 0.5 * dt * Ft[j] / masses[j]
            x[j] += dt * vhalf[j]
        Ft, Ut[i]= gravity(x, masses)
        for j in range(N):
            v[j] = vhalf[j] + 0.5 * dt * Ft[j] / masses[j]
            kinetic[i] += 0.5 * masses[j] * np.sum(v[j]**2) 
        totalE[i] = kinetic[i] + Ut[i] 
            
    return x, v, kinetic, Ut, totalE