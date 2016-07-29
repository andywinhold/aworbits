# Rendezvous with Ramageddon #

As outlined in the project proposal, this project sought to accurately simulate parts of the solar system along with a foreign object in a hyperbolic orbit passing by the sun. To run our code the files getcoordinates.py and rama.py need to be imported. In addition the four commands below are needed to run the simulation.

```getcoordinates.get_perifocal()``` - returns two arrays containing the positions and velocities of the sun and all the planets in the solar system. These are in perifocal coordinates, meaning that they are centered at the focus of the orbit but do not extend into the z-axis.

```getcoordinates.get_heliocentric(r, v)``` - requires the positions and velocites returned by the previous command. This functions to transform perifocal coordinates in to heliocentric coordinates. Heliocentric coordinates have the sun at their focus and accurately model the angle of inclination of each object in the solar system. The function returns two arrays containg the positions and velocities respectively.

```getcoordinates.add_Rama(r, v)``` - This function adds Rama, our alien object to the position and velocity arrays that were returned from ```getcoordinates.get_heliocentric(r, v)```. The position and velocity of Rama can be any number but must be in the format of 1x3 arrays, where x = [0], y = [1], and z = [2]. The function returns 2 arrays, one for position and velocity, that contain the sun, all of the planets, and rama.

```rama.dynamics(r, v, dt=0.01, tmax=10)``` - This function runs the simulation. It requires the arrays generated above but it's flexible enough that it can be run with and without Rama. It also requires a time step (default is 0.01 years) and a maximum length of time (default is 10 years). It returns five arrays, the positions of each object, the velocities, the kinetic energy, the potential energy, and the total energies of the system. Each array is populated for every timestep.


There are 4 files in the submission directory. The function and purpose of each of the four files is listed below.

### rama.py
This file contains our force and velocity verlet functions.

### getcoordinates.py
This file contains the orbital elements of each solar system object other than Rama. It also contains the functions necessary to convert these orbital elements to perifocal and heliocentric coordinates. Lastly, it also contains the function to add rama to the arrays for heliocentric coordinates.

### visual_rama.py
This is the file that was used to visulize the simulation. It requires VPython to run properly. Once this is installed it can be run in a jupyter notebook or using the Vpython VIDLE.


### Run Rama.ipynb
This is a jupyter notebook that has been included in the submission to make running the simulation easier. It imports the necessary modules and contains the code to set up and run the simulation as outlined above. It prints a 3d graph with the orbits of the planets over the specified time period.
