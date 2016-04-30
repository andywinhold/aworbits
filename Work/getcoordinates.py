import numpy as np

#Units of mass are in solar mass
#Units of time are in years
#Units of position are in AU
#Units of velocity are in AU/year
#

mu = 3.964016 * 9.945193 #AU^3/(solarmass * yr^2)

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
            ])

#Semi-Major Axis in AU
smAxis = np.array([
            0.,          #Sun
            0.38709927,  #Mercury
            0.72333566,  #Venus
            1.00000261,  #Earth
            1.52371034,  #Mars
            5.200441899,  #Jupiter
            9.53667594,  #Saturn
            19.18916464, #Uranus
            30.06992276, #Neptune
            67.68871444,# Eris
            ])

# Eccentricity, dimensionless
eccentricity = np.array([
            0.,          #Sun
            0.20563593,  #Mercury
            0.00677672,  #Venus
            0.01671123,  #Earth
            0.09339410,  #Mars
            0.04838624,  #Jupiter
            0.05386179,  #Saturn
            0.04725744,  #Uranus
            0.00859048,  #Neptune
            0.44068      #Eris
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
            1.77004347,   #Neptune
            44.0445      #Eris
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
            -55.12002969,  #Neptune
            204.16         #Eris
            ])
mLong = np.deg2rad(mLong) #converting to radians

#Argument of perihelion in degrees
pLong = np.array([
            0.,            #Sun
            77.45779628,   #Mercury
            131.60246718,  #Venus
            102.93768193,  #Earth
            -23.94362959,  #Mars
            14.72847983,   #Jupiter
            92.59887831,   #Saturn
            170.95427630,  #Uranus
            44.96476227,   #Neptune
            187.1498689    #Eris ** Found by adding argument of perihelion
                           #and longitude of the ascending node, per JPL details.
            ])
pLong = np.deg2rad(pLong) #converting to radians

#Longitude of ascending node in degrees
aLong = np.array([
            0.,            #Sun
            48.33076593,   #Mercury
            76.67984255,   #Venus
            0.0,           #Earth
            49.55953891,   #Mars
            100.47390909,  #Jupiter
            113.66242448,  #Saturn
            74.01692503,   #Uranus
            35.906450258,  #Neptune
            131.78422574   #Eris
            ])
aLong = np.deg2rad(aLong) #converting to radians

system = {"Mass":mass, "Semi-Major Axis":smAxis, "Eccentricity":eccentricity, "Inclination":inclination, "Mean Longitude":mLong,"Perihelion Longitude":pLong, "Ascending Longitude":aLong}

def ang_mom(mu, a, e):
    """ Define angular momentum for state vector calculation.

    Inputs
    ------

    mu : gravitational parameter, units of AU^3/(solarmass * yr^2)

    a  : semi-major axis, units of AU

    e  : eccentricity, dimensionless
    
    Outputs
    -------

    h  : scalar, angular momentum of body, units in solarmass * AU/year^2
    """
    h = np.sqrt(mu*a*(1-e**2))
    return h

# True Anomaly (p.112 OMES, eqn 3.10)
def t_anom(E, e):
    """ Calculates the true anomaly of a body

    angular parameter that defines the position of a body moving along a Keplerian orbit. 
    It is the angle between the direction of periapsis (closest approach to cent. body) 
    and the current position of the body, as seen from the main focus of the ellipse 
    (the point around which the object orbits).

    Inputs
    ------

    E : eccentric anomaly, dimensionless

    e : eccentricity of orbiting body, dimensionless

    Outputs
    -------

    theta : angle, in radians
    """
    theta = 2*np.arctan(np.tan(E/2)*np.sqrt((1-e)/(1+e)))
    return theta

#Calculate angular Momentum, amomentum
h = ang_mom(mu, system["Semi-Major Axis"],system["Eccentricity"])

#Calculate mean Anomaly, m_anomaly
m_anomaly = system["Mean Longitude"] - system["Perihelion Longitude"]
   
#Calculate eccentric Anomaly, E
tol = 1e-5
E = m_anomaly + system["Eccentricity"] * np.sin(m_anomaly)
deltaE = np.zeros_like(E)
count = 0
while np.abs(np.linalg.norm(deltaE)) > tol:
    deltaE = (M - (E - system["Eccentricity"] * np.sin(E))) / (1 - system["Eccentricity"] * np.cos(E))
    E = E + deltaE
    count += 1
    if count > 1000:
        print("E did not converge to a solution after 1000 iterations.")
        break
        
#Calculate true Anomaly, real_anomaly
real_anomaly = t_anom(E, system["Eccentricity"])

#Calculate argument of perihelion, a_long
a_long = system["Perihelion Longitude"] - system["Ascending Longitude"]

def get_perifocal():
    """Calculates perifocal coordinates from orbital elements.
    
    Perifocal coordinates do not take the z-axis into account, only the x and y orientation of the object.
    
    Outputs
    -------

    rp : N x 3 array
        array of sun-centered coordinates for each object in the system dictionary, units in AU
        
    vp : N x 3 array
        velocities of each object in the solar system, units in AU/year
    """
    
    # 4.37 position array
    r_ar = np.zeros((len(system["Mass"]),3))
    r_ar[:,0] = np.cos(a_long)
    r_ar[:,1] = np.sin(a_long)
    
    #4.38 velocity array
    v_ar = np.zeros((len(system["Mass"]),3))
    v_ar[:,0] = -np.sin(a_long)
    v_ar[:,1] = system["Eccentricity"] + np.cos(a_long)
    
    # equations 4.37 and 4.38 in OMES p. 173
    rp = np.zeros((len(system["Mass"]),3))
    vp = np.zeros((len(system["Mass"]),3))

    rp[:,0] = (h**2/mu) * (1/(1 + system["Eccentricity"]*np.cos(a_long))) * r_ar[:,0]
    rp[:,1] = (h**2/mu) * (1/(1 + system["Eccentricity"]*np.cos(a_long))) * r_ar[:,1]
    vp[1:,0] = (mu/h[1:]) * v_ar[1:,0]
    vp[1:,1] = (mu/h[1:]) * v_ar[1:,1]

    return rp, vp

def get_heliocentric(r, v):
    """Transforms perifocal coordinates into heliocentric cordinates.
    
    Heliocentric coordinates are oriented with respect to the ecliptic plane of the solar system 
    and correctly model the solar system.
    
    Outputs
    -------

    ecliptic_r : N x 3 array
        array of sun-centered coordinates for each object in the system dictionary 
        in heliocentric frame, units of AU
        
    ecliptic_v : N x 3 array
        velocities of each object in the solar system in heliocentric frame, units of AU/year
    """
    ecliptic_r = np.zeros((len(system["Mass"]),3))
    ecliptic_v = np.zeros((len(system["Mass"]),3))
    
    cosw = np.cos(a_long) #small omega
    sinw = np.sin(a_long)
    cosO = np.cos(system["Ascending Longitude"]) #big omega
    sinO = np.sin(system["Ascending Longitude"])
    cosI = np.cos(system["Inclination"]) #i
    sinI = np.sin(system["Inclination"])
    
    #Equations derived from rotation matrix
    ecliptic_r[:,0] = (cosO * cosw - sinO * sinw * cosI) * r[:,0] + (-cosO * sinw - sinO * cosw * cosI) * r[:,1]
    ecliptic_r[:,1] = (sinO * cosw + cosO * sinw * cosI) * r[:,0] + (-sinO * sinw + cosO * cosw * cosI) * r[:,1]
    ecliptic_r[:,2] = sinw * sinI * r[:,0] + cosw * sinI * r[:,1]
    
    #Equations derived from rotation matrix
    ecliptic_v[:,0] = (cosO * cosw - sinO * sinw * cosI) * v[:,0] + (-cosO * sinw - sinO * cosw * cosI) * v[:,1]
    ecliptic_v[:,1] = (sinO * cosw + cosO * sinw * cosI) * v[:,0] + (-sinO * sinw + cosO * cosw * cosI) * v[:,1]
    ecliptic_v[:,2] = sinw * sinI * v[:,0] + cosw * sinI * v[:,1]

    return ecliptic_r, ecliptic_v

    
    
