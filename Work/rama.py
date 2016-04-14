import numpy as np


#============================================================
# Parameters of the problem
#------------------------------------------------------------
#

G_gravity = 4*np.pi**2

#Masses in solar masses
mass = np.array([
            1.,                 #Sun
            0.330104/(1.989e6), #Mercury
            4.86732/(1.989e6),  #Venus
            5.97219/(1.989e6),  #Earth
            0.641693/(1.989e6), #Mars
            1898.13/(1.989e6),  #Jupiter
            568.319/(1.989e6),  #Saturn
            86.8103/(1.989e6),  #Uranus
            102.410/(1.989e6),  #Neptune
            ])

#Semi-Major Axis in AU
smAxis = np.array([
            0.,          #Sun
            0.38709927,  #Mercury
            0.72333566,  #Venus
            1.00000261,  #Earth
            1.52371034,  #Mars
            5.20288700,  #Jupiter
            9.53667594,  #Saturn
            19.18916464, #Uranus
            30.06992276  #Neptune
            ])

# Eccentricity
eccentricity = np.array([
            0.,          #Sun
            0.20563593,  #Mercury
            0.00677672,  #Venus
            0.01671123,  #Earth
            0.09339410,  #Mars
            0.04838624,  #Jupiter
            0.05386179,  #Saturn
            0.04725744,  #Uranus
            0.00859048   #Neptune
            ])

#Inclination in degrees
inclination = np.array([
            0.,          #Sun
            7.00497902,  #Mercury
            3.39467605,  #Venus
            -0.00001531, #Earth
            1.84969142,  #Mars
            1.30439695,  #Jupiter
            2.48599187,  #Saturn
            0.77263783,  #Uranus
            1.77004347   #Neptune
            ])
inclination = np.deg2rad(inclination) #converting to radians

#Mean Longitude in degrees
mLong = np.array([
            0.,            #Sun
            252.25032350,  #Mercury
            181.97909950,  #Venus
            100.46457166,  #Earth
            -4.55343205,   #Mars
            34.39644051,   #Jupiter
            49.95424423,   #Saturn
            313.23810451,  #Uranus
            -55.12002969   #Neptune
            ])
mLong = np.deg2rad(mLong) #converting to radians

#Longitude of perihelion in degrees
pLong = np.array([
            0.,            #Sun
            77.45779628,   #Mercury
            131.60246718,  #Venus
            102.93768193,  #Earth
            -23.94362959,  #Mars
            14.72847983,   #Jupiter
            92.59887831,   #Saturn
            170.95427630,  #Uranus
            44.96476227    #Neptune
            ])
pLong = np.deg2rad(pLong) #converting to radians

#Longitude of ascending node in degrees
aLong = np.array([
            0.,            #Sun
            48.33076593,   #Mercury
            76.67984255,  #Venus
            0.0,  #Earth
            49.55953891,  #Mars
            100.47390909,   #Jupiter
            113.66242448,   #Saturn
            74.01692503,  #Uranus
            131.78422574    #Neptune
            ])
aLong = np.deg2rad(aLong) #converting to radians

#Dictionary of solar system properties
system = {"Mass":mass, "Semi-Major Axis":smAxis, "Eccentricity":eccentricity, "Inclination":inclination, "Mean Longitude":mLong,
          "Perihelion Longitude":pLong, "Ascending Longitude":aLong}

def getCoordinates():
    """Calculates cartesian coordinates of Heliocentric solar system for each solar system object.

    Returns
    -------
    coordinates : array
        Returns coordinates of every object as (x,y,z) array, with the sun at the center.
    """
    tol = 1e-5
    coordinates = np.zeros((len(system["Mass"]), 3))
    
    omega = np.degrees(system["Perihelion Longitude"]) - np.degrees(system["Ascending Longitude"])
    meanAnomaly = (system["Mean Longitude"] - system["Perihelion Longitude"]) % np.pi
    E = meanAnomaly + np.deg2rad(system["Eccentricity"]) * np.sin(meanAnomaly)
    deltaE = np.zeros_like(E)
                           
    count = 0
    while np.abs(np.linalg.norm(deltaE)) > tol:
        deltaE = (M - (E - np.deg2rad(system["Eccentricity"]) * np.sin(E))) / (1 - system["Eccentricity"] * np.cos(E))
        E = E + deltaE
        count += 1
        if count > 1000:
            print("E did not converge to a solution after 1000 iterations.")
            break
    
    #calculates heliocentric coordinates in orbital plane
    coordinates[:, 0] = system["Semi-Major Axis"] * (np.cos(E) - system["Eccentricity"])
    coordinates[:, 1] = system["Semi-Major Axis"] * np.sin(E) * np.sqrt(1 - system["Eccentricity"]**2)
                  
    #calculates heliocentric coordinates in ecliptic plane
    coordinates[:, 0] = (np.cos(omega) * np.cos(system["Ascending Longitude"]) - 
                         np.sin(omega) * np.sin(system["Ascending Longitude"]) * 
                         np.cos(system["Inclination"])) * coordinates[:, 0] + (np.sin(omega) * np.cos(system["Ascending Longitude"]) - np.cos(omega) * np.sin(system["Ascending Longitude"]) * np.cos(system["Inclination"])) * coordinates[:, 1]
    coordinates[:, 1] = (np.cos(omega) * np.sin(system["Ascending Longitude"]) + 
                         np.sin(omega) * np.cos(system["Ascending Longitude"]) * 
                         np.cos(system["Inclination"])) * coordinates[:, 0] + (-np.sin(omega) * np.sin(system["Ascending Longitude"]) + np.cos(omega) * np.cos(system["Ascending Longitude"]) * np.cos(system["Inclination"])) * coordinates[:, 1]
    coordinates[:, 2] = (np.sin(omega) * np.sin(system["Inclination"])
                        ) * coordinates[:, 0] + (np.cos(omega) * np.sin(system["Inclination"])) * coordinates[:, 1]
    
    return coordinates
    

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
            dF =  -rhat * G_gravity * masses[i] * masses[j] / (rmag**2)
            force[j] += dF
            force[i] += -dF
            U += -G_gravity * masses[i] * masses[j] / rmag    
            
    return force, U

def dynamics(x0, v0, dt, masses, tmax=10):
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
    x = np.zeros((nsteps/10,N,3))
    dx = np.copy(x0)
    v = np.copy(v0)
    vhalf = np.zeros((N,3))
    Ut = np.zeros(nsteps)
    kinetic = np.zeros(nsteps)
    totalE = np.zeros(nsteps)
    
    Ft, Ut[0]= gravity(dx, masses)
    
    for i in range(nsteps):
        for j in range(N):
            vhalf[j] = v[j] + 0.5 * dt * Ft[j] / masses[j]
            dx[j] += dt * vhalf[j]
        Ft, Ut[i]= gravity(dx, masses)
        for j in range(N):
            v[j] = vhalf[j] + 0.5 * dt * Ft[j] / masses[j]
            kinetic[i] += 0.5 * masses[j] * np.sum(v[j]**2) 
        totalE[i] = kinetic[i] + Ut[i]
        if i%10 == 0:
            x[int(i/10)] = dx
            
    return x, v, kinetic, Ut, totalE