#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# Molecular Dynamics of the Lennard Jones Fluid
# Skeleton code --- incomplete.
#
# Written by Oliver Beckstein for ASU PHY494
# http://asu-compmethodsphysics-phy494.github.io/ASU-PHY494/
# Placed into the Public Domain

import time
import numpy as np
import mdIO
import mdInit

# constants (note: our units for energy are kJ/mol!)
kBoltzmann = 1.3806488e-23   # J/K
RgasConst = 8.3144621e-3     # J/(mol*K)  
Rgas = 0.0083144621          # KJ/(mol*K)  <--- this is what we use
epsilon = 0.99774            # 120 K times kb
sigma = 0.34                 # nm

# parameters for noble gases
parameters = {
    'Ar': {
        'mass': 39.948,       # u
        'epsilon': 0.99774,   # kJ/mol = 120 K * kB
        'sigma': 0.34,        # nm
        }
    }

def minimage(x, box):
    """Minimum image coordinates of x in orthorhombic box."""
    # as in Frenkel and Smit
    return x - box*np.rint(x/box)
        
def initial_velocities(atoms, T0):
    """Generate initial velocities for *atoms*.
    - random velocities
    - total momentum zero
    - kinetic energy corresponds to temperature T0
    Parameters
    ----------
    atoms : list
         list of atom names, e.g. `['Ar', 'Ar', ...]`
    T0 : float
         initial temperature in K
    Returns
    -------
    velocities : array
         Returns velocities as `(N, 3)` array.
    """
    Natoms = len(atoms)
    masses = np.array([parameters[name]['mass'] for name in atoms])
    v = mdInit.random_velocities(Natoms)
    v[:] = mdInit.remove_linear_momentum(v, masses)
    return mdInit.rescale(v, T0, masses), masses

def lj(coordinates, box):
    """Calculates total force on each atom using periodic boundary conditions.
    
    Parameters
    ----------
    coordinates : array
        [x,y,z] coordinates for every atom
    box : array
        array of [x,y,z] lengths for lattice area
        
    Returns
    -------
    force : array
        Returns net force on each atom as '(N, 3)' array and total potential as scalar.
    """
    
    N = len(coordinates) #number of atoms
    r = np.zeros((N, N, len(box))) #distance between atoms
    force = np.zeros((N, len(box))) #force between atoms
    U = 0 #potential energy
    rcut = 3 * parameters['Ar']['sigma']
    Ucut = 4 * parameters['Ar']['epsilon'] * ((parameters['Ar']['sigma']/rcut)**12 
                                             - (parameters['Ar']['sigma']/rcut)**6)
    
    for i in range(0, N - 1):
        for j in range(i + 1, N):
            #calculates distance vector between each atom
            r[i,j] = coordinates[j] - coordinates[i]
            r[i,j] = minimage(r[i,j], box)
            
            #force calculation
            rmag = np.sqrt(r[i,j,0]**2 + r[i,j,1]**2 + r[i,j,2]**2)
            if rmag > rcut: #if distance is farther than 3 sigma then force is zero.
                continue
            else:
                virial = 0 #initialize virial sum
                rhat = r[i,j] / rmag
                ratio = (parameters['Ar']['sigma'] / rmag)**6
                dF =  (rhat / rmag) * 48 * parameters['Ar']['epsilon'] * ratio * (ratio - 0.5)
                force[j] += dF
                force[i] += -dF
                U += 4 * parameters['Ar']['epsilon'] * ratio * (ratio - 1) - Ucut
                virial += np.dot(r[i, j], dF) #this is the sum for the pressure estimator
    return force, U, virial

def dynamics(atoms, x0, v0, dt, box, V, nsteps=10, filename="trajectory.xyz"):
    """Integrate equations of motions
    Parameters
    ----------
     atoms : list
         list of atom names
     x0 : array
         Nx3 array containing the starting coordinates of the atoms.
         (Note that x0 is changed and at the end of the run will
         contain the last coordinates.)
     v0 : array
         Nx3 array containing the starting velocities, eg np.zeros((N,3))
         or velocities that generate a given temperature
     dt : float
         integration timestep in ps (e.g. 0.001 ps = 1 fs)
     V : float
         Total volume
     nsteps : int, optional
         number of integrator time steps
     filename : string
         filename of trajectory output in xyz format
    Writes coordinates to file `filename`.
    """
    N = len(atoms)
    x = np.copy(x0)
    v = np.copy(v0)
    vhalf = np.zeros((N,len(box)))
    Ut = np.zeros(nsteps)
    kinetic = np.zeros(nsteps)
    totalE = np.zeros(nsteps)
    temp = np.zeros(nsteps)
    pressure = np.zeros(nsteps)
    virial = np.zeros(nsteps)
    masses = np.array([parameters[name]['mass'] for name in atoms])
    
    Ft, Ut[0], virial[0] = lj(x0, box) #initialize force and potential energy
    
    mdIO.write_xyz(filename, atoms, x, box) #write xyz file using provided code
    
    for i in range(nsteps):
        for j in range(N):
            vhalf[j] = v[j] + 0.5 * dt * Ft[j] / parameters[atoms[j]]['mass']
            x[j] += dt * vhalf[j]
            x[j] = minimage(x[j], box)
        for j in range(nsteps - (nsteps % 10)): # this for loop writes only every 10 steps to the file
            if i == (10*j):
                mdIO.write_xyz_frame(open(filename, "a"), atoms, x, box, i) #append new frame to xyz file
        Ft, Ut[i], virial[i] = lj(x, box)
        for j in range(N):
            v[j] = vhalf[j] + 0.5 * dt * Ft[j] / parameters[atoms[j]]['mass']
            kinetic[i] += 0.5 * parameters[atoms[j]]['mass'] * np.sum(v[j]**2) 
            temp[i] = mdInit.kinetic_temperature(v, masses)
            pressure[i] = (Rgas * N * (10**12) * temp[i] - virial[i]/3)/V
        totalE[i] = kinetic[i] + Ut[i] 
            
    return x, v, kinetic, Ut, totalE, temp, pressure
        
if __name__ == "__main__":
    import lattice
    import time
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    #------------------------------------------------------------
    # initial conditions
    #------------------------------------------------------------
    density = 1.374   #g/cm**3
    initial_temperature = 94.4  #Kelvin
    num_atoms = 64
    atom_name = "Ar"
    dt = 0.01  # timestep in ps, 10 fs is a good starting value
    runtime = 10 # total run time in ps
    nsteps = int(runtime/dt)

    #------------------------------------------------------------
    # initialization
    #------------------------------------------------------------
    atoms, coordinates, box, V = lattice.generate_lattice(density, num_atoms, atomname=atom_name,
                                                       lattice="cubic")
    velocities, masses = initial_velocities(atoms, initial_temperature)


    #------------------------------------------------------------
    # MD
    #------------------------------------------------------------
    start = time.time()
    finalx, finalv, k, U, E, temp, pressure = dynamics(atoms, coordinates, velocities, dt, box, V, nsteps=nsteps, filename="trajectory.xyz")
    stop = time.time()
    wallT = stop - start
    
    #------------------------------------------------------------
    # Analysis
    #------------------------------------------------------------
    #Time averaged values
    avg_Temp = np.sum(temp) / nsteps
    avg_U = (np.sum(U) / nsteps) / num_atoms
    avg_E = ((np.sum(E)) / nsteps) / num_atoms
    avg_pressure = np.sum(pressure) / nsteps
    
    #Standard deviations
    sumT = np.sum((temp - avg_Temp)**2) / nsteps
    sd_T = np.sqrt(sumT)
    sumU = np.sum((U / num_atoms-avg_U)**2) / nsteps
    sd_pot = np.sqrt(sumU)
    sumE = np.sum((E / num_atoms - avg_E)**2) / nsteps
    sd_tot = np.sqrt(sumE)
    sumP = np.sum((pressure - avg_pressure)**2) / nsteps
    sd_P = np.sqrt(sumP)
    
    #Energy drift
    initalE = U[0] + k[0]
    e_sum = np.sum(abs(1 - ((E) / initalE)))
    e_drift = e_sum / nsteps
    
    print()
    print("--- Final Results " + 42 * "-")
    print("Average Temperature (K):               %f" % avg_Temp)
    print("   Standard Deviation (K):             %f" % sd_T)
    print("Average Pressure (N/(mol*nm**2)):               %f" % avg_pressure)
    print("   Standard Deviation (N/(mol*nm**2)):             %f" % sd_P)
    print("Average Potential per atom (kJ/mol):  %f" % avg_U)
    print("   Standard Deviation (kJ/mol):        %f" % sd_pot)
    print("Average Energy per atom (kJ/mol):      %f" % avg_E)
    print("   Standard Deviation (kJ/mol):        %f" % sd_tot)
    print("The energy drift is (kJ/mol):          %f" % e_drift)
    print(60*"-")
    
    
    #plots
    
    t = np.linspace(0, runtime, nsteps)
        
    #Initial Positions
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2], color='blue')
    plt.title("Initial Frame")
    ax1.set_xlabel('Length (nm)')
    ax1.set_ylabel('Length (nm)')
    ax1.set_zlabel('Length (nm)')
    ax1.figure.savefig("Initial Frame.png")
    
    #Final Positions
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(finalx[:,0],finalx[:,1], finalx[:,2], color='red')
    plt.title("Final Frame")
    ax2.set_xlabel('Length (nm)')
    ax2.set_ylabel('Length (nm)')
    ax2.set_zlabel('Length (nm)')
    ax2.figure.savefig("Final Frame.png")
    
    #Total Energy Graph
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    plt.plot(t, k, color='red', label=r"$Kinetic$")
    plt.plot(t, U, color='blue', label=r"$Potential$")
    plt.plot(t, E, color='green', label=r"$H(t) = T_{kin}+U$") # H is Hamiltonian
    plt.legend(loc='best')
    plt.title("System Energies", color='black')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (kJ/mol)')
    ax3.figure.savefig("System Energy.png")

    #Energy per particle
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    plt.plot(t, k/num_atoms, color='red', label=r"$Kinetic$")
    plt.plot(t, U/num_atoms, color='blue', label=r"$Potential$")
    plt.plot(t, E/num_atoms, color='green', label=r"$H(t) = T_{kin}+U$") # H is Hamiltonian
    plt.legend(loc='best')
    plt.title("Energy per Particle", color='black')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (kJ/mol)')
    ax4.figure.savefig("Energy per particle.png")
    
    #Instantaneous temp
    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111)
    plt.plot(t, temp)
    plt.title("Temperature (K)", color='black')
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    ax5.figure.savefig("Temperature.png")
    
    #Instantaneous pressure
    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111)
    plt.plot(t, pressure)
    plt.title("Pressure (N/(mol*nm**2))", color='black')
    plt.xlabel('Time (ps)')
    plt.ylabel('Pressure (N/(mol*nm**2))')
    ax6.figure.savefig("Pressure.png")

    #------------------------------------------------------------
    # Performance
    #------------------------------------------------------------
    stepDay = (nsteps/wallT) * 3600 * 24 #steps per day
    nanoDay = (runtime / wallT) * 3600 * 24 / 1000 #nanoseconds per day
    
    print()
    print("--- Performance " + 44 * "-") 
    print("Wall time (s):                               %f" % wallT)
    print("Real time per simulated time step (s/step):  %f" % (wallT / nsteps))
    print("Steps per wall time (steps):                 %i" % nsteps)
    print("Simulated picoseconds per wall time (ps):    %d" % runtime)
    print("Steps per day (steps/day):                   %d" % stepDay)
    print("Simulated nanoseconds per day (ns/day):      %f" % nanoDay)
    print(60*"-")

    
    
    
    