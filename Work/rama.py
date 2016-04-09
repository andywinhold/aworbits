import numpy as np

G_gravity = 4*np.pi**2

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
        Returns net force on each solar system object as '(N, 3)' array and total potential as scalar.
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
            U += -G_gravity * masses[i] * masses[j] / rmag    
            
    return force, U