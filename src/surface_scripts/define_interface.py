#Find the interface of the water surface/air surface using the last frame from simulation
#Then deletes surfactant molecules that are > certain distance from interface 

#!/usr/bin/env python 

import sys
import numpy as np
import readDCD

#read in xyz coordinates of water from last step in simulation 
water_coordinates = open('xyz_water.txt', 'r')
surf_coordinates = open('xyz.txt', 'r')

dum=water_coordinates.readlines()
lfile=[]
for i in dum:
    lfile.append(i.split())
lammps_dump=np.array(lfile[9:],float)
xyz=lammps_dump[:,2:5]
mass=lammps_dump[:,1]

z_max = 60.0
z_min = -150 
delta_x = 60
delta_y = 60
delta_z = 4

#z_value = np.array(xyz[:,2],float)

w_beads = []
i=0

while z_max > 25:
    z_max= z_max-delta_z
    print z_max
    for z_value in xyz[:,2]:
        if (z_max-delta_z) < z_value < z_max:
            print z_value
            w_beads.append(z_value)
