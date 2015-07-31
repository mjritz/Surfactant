#Attempt to find interface using each frame from dcd file 

#!/usr/bin/env python 

import sys
import numpy as np
import readDCD
import variables
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

f1= open ('out_file.txt','wa')
out_stream = f1

f2 = open('trimmed_system.xyz', 'wa') #Open in write and append format
outstream = f2

#Read in dcd file, define local water surface, calculate distance from surfactant 
xbox, ybox, zbox, w_filename, w_mol, pack_tol, bead_dist, boxx, boxy, boxz, bead_per, equil_time, surf_x, surf_y, surf_zlow, surf_zhi, per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd, prod_dcd = variables.read_variables('input_file.txt')

nAt,BoxD_frames,BoxAngle_frames,xyz_frames = readDCD.read_DCD(equil_dcd)

#define total number of molecules, number of water molecules, and number of surfactants
total_molecules = nAt
water_molecules = int(w_mol)
bead_per_surfactant = int(bead_per)
surf_molecules = (total_molecules-water_molecules)/bead_per_surfactant
xyz_final_frame = xyz_frames[-1]

water = []
w_molecule=0
for w_molecule in range(w_mol):
    water.append(xyz_final_frame[w_molecule,:])
    w_molecule+=1
water_array =np.array(water, float)

surfactants = []
surf_mol= w_mol+1
#print surf_mol
for surf_mol in range(w_mol, total_molecules):
        surfactants.append(xyz_final_frame[surf_mol,:])
        surf_mol+=1
        surf_array =np.array(surfactants, float)

#FOR DEBUGGING

#print xyz_final_frame
#print xyz_final_frame.shape
#print water_array
#print water_array.shape

#create array for water molecules and surfactant molecules 
#x = []
#y = []
#z = []
#x.append(xyz_final_frame[:,0])
#y.append(xyz_final_frame[:,1])
#z.append(xyz_final_frame[:,2])
#x_array = np.array(x, float)
#y_array = np.array(y, float)
#z_array = np.array(z, float)
#print x_array.shape


#concatenate the four [N,1] arrays together to form and [N,4] array 
#xyz_all_array = np.concatenate((x_array, y_array, z_array), axis=0)

#print xyz_all_array.shape
#print xyz_all_array

g = 'original xyz file\n' 
for (i,(x,y,z)) in enumerate(water_array):
    g += "%5i %8.3f%8.3f%8.3f\n" % (i+1, x, y, z)
for (i,(x,y,z)) in enumerate(surf_array):
    g += "%5i %8.3f%8.3f%8.3f\n" % (i+1, x, y, z)
g+= "End\n"

#out_stream.write(g)
#out_stream.close()

def find_interface(boxx, boxy, boxz, water_array):

    def wrapcoord(xyzcom,boxD,box_center=[0,0,0]):
        #shift box to center
        xyzcom-=box_center
        xyzcom-=boxD*np.around(xyzcom/boxD)
        return xyzcom
    
    D=3
    boxD=np.zeros(D)
    boxD[0]=float(boxx)*2
    boxD[1]=float(boxy)*2
    boxD[2]=float(boxz)

    #wrap water coordinates to box
    xyz_wat=wrapcoord(water_array,boxD,box_center=[0,0,0])

    #calculate distribution of waters as a function of z
    delta_z =0.50
    wat_freq,bins=np.histogram(xyz_wat[:,2],bins=np.linspace(-float(boxD[2]/2.0), float(boxD[2]/2.0)+1, ((float(boxD[2]/2.0)+1)-(-float(boxD[2]/2.0)))/delta_z))
    #print wat_freq

    #find min and max of water surface, cutoff used arbitrarily to avoid pure   vapor
    #assumes that there is a value of more than 20 waters in an 80x80x1 A slab (from ND, fine assumption), 20/(80x80x1) = 0.003125
    #Using density, NA, g/mol, 1bead/molecule => 1bead/59.88Ang => should be more than 30 molecules in 60x60x1 slab, 30/(60x60x1)
    vapor_cut=.01388*boxD[0]*boxD[1]
    minl=np.argmax(wat_freq>vapor_cut)
    minr=np.size(wat_freq)-np.argmax(wat_freq[::-1]>vapor_cut)
    #print minl, minr

    #estimate bulk water number density from middle 20 A of water slab
    slabcenter=(minl+minr)/2.0
    wat_density=np.average(wat_freq[slabcenter-10:slabcenter+10])
    #water density should be around 60 for 60x60x1 slab 
    #print wat_density

    #find interface
    interface_l=np.argmax(wat_freq>wat_density/2.0)
    interface_r=np.size(wat_freq)-np.argmax(wat_freq[::-1]>wat_density/2.0)-1
    #print interface_l, interface_r

    lower_interface = -int(boxD[2]/2.0)+(delta_z*interface_l)
    upper_interface = -int(boxD[2]/2.0)+(delta_z*interface_r)
    return lower_interface, upper_interface    

lower_interface, upper_interface = find_interface(boxx, boxy, boxz, water_array) 

#print surf_array.shape
accepted_surf = []
mol_array = surf_array.reshape(-1, bead_per_surfactant, 3)
#print mol_array.shape

#search through the first beads to find  z values inside  of the range and copy
#molecule to a new list 
acceptable_range = float(bead_dist)*(float(bead_per_surfactant)-2) 
maximum = upper_interface+acceptable_range
minimum = lower_interface-acceptable_range
imol = 0
for zval in  mol_array[:,0,2]:
    if zval<maximum and zval>minimum:
        accepted_surf.append(mol_array[imol,:,:])
    else:
        pass
    imol +=1

print "\n\nlower interface:",lower_interface,";", "upper interface:",upper_interface
print "accepted surfactants", acceptable_range,"angstroms away from lower and upper interface"
accepted_surf_array = np.array(accepted_surf, float)
num_surf_mol = accepted_surf_array.shape[0]
good_surf_array = accepted_surf_array.reshape(-1, 3)
deleted = surf_molecules - good_surf_array.shape[0]/bead_per_surfactant
print  "started with:",surf_molecules,"surfactants;", "deleted:",surf_molecules - good_surf_array.shape[0]/bead_per_surfactant,"surfactants\n\n"

pdb_array = np.concatenate((water_array,good_surf_array), axis=0)
number_of_atoms = str(pdb_array.shape[0])
comment = "Created after trimming surfactants"

#Writing out xyz file -----------------------------------------------------
s=number_of_atoms
s+= "\n"
s+=comment
s+="\n"
for (i, (x, y, z)) in enumerate(pdb_array):
    s+="%3i%12.3f%12.3f%12.3f\n" % (i+1, x, y, z)

outstream.write(s)
outstream.close()

p = "-------------------------------------------------------------------"
p += "\n\nlower interface: %5.3f; upper interface: %5.3f\n"  %(lower_interface ,upper_interface)
p += "accepted surfactants %5.3f angstroms away from lower and upper interface\n" %(acceptable_range)
p += "started with: %i surfactants; deleted: %i surfactants\n"  % (surf_molecules, deleted)
p+= "-------------------------------------------------------------------\n\n"
out_stream.write(p)
out_stream.close()
