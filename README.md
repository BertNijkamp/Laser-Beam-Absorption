# Laser-Beam-Absorption

This project calculates the laser propagation and absorption for a powder bed and a particle doublet. It also includes the code for the random walk.

Most files are functions which are used in the calculations, except for three files:

* *Propagation3Dbed.m*
* *RandomWalkPropagation.m*
* *PG_60_1.txt*

*Propagation3Dbed.m* is the file which is used to calculate the propagation and absorption in a particle bed, as well as rendering the result. It has several settings, which can be found in the file itself. Be aware that you should only set ShowLaser to 1 if you use a low NumberOfRays.

*RandomWalkPropagation.m* simulates a semi-random walk and it's variables can again be changed in the file itself.

*PG_60_1.txt* is a file used in *Propagation3Dbed.m*. It consists of a list of x, y and z coordinates of the particles in a bed and their radii, and is generated using MercuryDPM. The bed is monodisperse with a particle diameter of 60 μm, has bed dimensions of 1x1x0.5 mm and a volume fraction of 0.6.

The code for calculating the propagation and absorption for a particle doublet is *Propagation3Ddoublet.m*. This is a function and has the NumberOfRays, the distance between the particles and the original diameter as an input and outputs the absorption for each particle.

The other files do the following:

* *DiameterIncrease.m* calculates the increased diameter of two overlapping particles based on ther diameter and distance from eachother
* *DrawAbsorption3D.m* renders the spherical particles with a colour based on their absorption
* *DrawParticles3D.m* renders the spherical particles with a colour based on their diameter
* *GaussianIntensity.m* calculates the starting energy for all rays in a curcular pattern with a total energy of 1
* *InParticle3D.m* checks if a new laser position is inside one or more particles
* *MatrixToNumber.m* converts a list of particles to a single number
* *NumberToMatrix.m* converts a single number to a list of particles
* *PointsInCircle.m* randomly distributes points in a circle
* *Refraction3D.m* calculates everything that happens during refraction and reflection
