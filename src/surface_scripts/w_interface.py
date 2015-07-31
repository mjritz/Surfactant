#Script written by Nathan Duff
#Useful example of how to use definitions and executing a definition 
#Also shows his way of defining the interface, wrap coordinates, write out lammps file, read lammps file 


import numpy as np
import sys
import time
import bisect
import os


def write_lammps_dump(xyz,lammps_dump,header,dump_file_name):
  s=header
  for i,(pos) in enumerate(xyz):
    s+= ('%d %20.15g %20.15g %20.15g\n') % (lammps_dump[i,0],pos[0],pos[1],pos[2])
    f=open(dump_file_name,'w')
    f.write(s)
    f.close()
def read_lammps_dump(dump_file_name):
  f=open(dump_file_name,'r')
  dum=f.readlines()
  lfile=[]
  for i in dum:
    lfile.append(i.split())
  lammps_dump=np.array(lfile[9:],float)
  xyz=lammps_dump[:,2:5]
  mass=lammps_dump[:,1]
  return xyz,mass,lammps_dump
#wrap coordinates
def wrapcoord(xyzcom,boxD,box_center=[0,0,0]):
  #shift box to center
  xyzcom-=box_center
  xyzcom-=boxD*np.around(xyzcom/boxD)
  return xyzcom

#def main():
D=3
boxD=np.zeros(D)
  #set variables from command line inputs
dump_file=sys.argv[1]
mol_size=int(sys.argv[2])
boxD[0]=float(sys.argv[3])
boxD[1]=float(sys.argv[4])
boxD[2]=float(sys.argv[5])
    
accepted_surf = [] 

  #read lammps dump file of surfactant molecules
xyz,mass,lammps_dump=read_lammps_dump(dump_file)

  #make lammps header for output. Everything except number of atoms and ITEM ATOMS has arbitrary values
lammps_dump_header=("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"+str(np.size(lammps_dump,axis=0))+"\nITEM: BOX BOUNDS pp pp pp\n"
"-40 40\n-40 40\n-150 150\nITEM: ATOMS id xu yu zu\n")
  
  #find interface of water slab

  #read in coordinates of water from lammps dump file
xyz_wat,mass_wat,lammps_dump_wat=read_lammps_dump('xyz_water.txt')


#calculate distribution of waters as a function of z
wat_freq,bins=np.histogram(xyz_wat[:,2],range(-int(boxD[2]/2.0),int(boxD[2]/2.0)+1,1))
print wat_freq

  #find min and max of water surface, cutoff used arbitrarily to avoid pure vapor
  #assumes that there is a value of more than 20 waters in an 80x80x1 A slab
vapor_cut=.003125*boxD[0]*boxD[1]
print vapor_cut
minl=np.argmax(wat_freq>vapor_cut)
minr=np.size(wat_freq)-np.argmax(wat_freq[::-1]>vapor_cut)

  #estimate bulk water number density from middle 20 A of water slab
slabcenter=(minl+minr)/2.0
wat_density=np.average(wat_freq[slabcenter-10:slabcenter+10])

  #find interface
interface_l=np.argmax(wat_freq>wat_density/2.0)
interface_r=np.size(wat_freq)-np.argmax(wat_freq[::-1]>wat_density/2.0)-1

  #create array for surfactant molecules
xyz_surf_array = np.array(xyz, float)
mol_array = xyz_surf_array.reshape(-1, mol_size, 3)
  
  #search through the first beads to find  z values inside  of the range and copy
  #molecule to a new list 
imol = 0
for zval in  mol_array[:,0,2]:
    if zval<minl and zval>-minr:
        accepted_surf.append(mol_array[imol,:,:])
    else:
        pass
    imol +=1

print interface_l 
print interface_r
accepted_surf_array = np.array(accepted_surf, float)
good_surf_array = accepted_surf_array.reshape(-1, 3)
xyz_water = np.array(xyz_wat, float)

array = np.concatenate((xyz_water,good_surf_array),axis=0)

  #write out shifted molecules to lammps dump files
write_lammps_dump(array,lammps_dump,lammps_dump_header,dump_file)
   
