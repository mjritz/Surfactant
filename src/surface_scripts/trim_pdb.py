#Trims surfactants above a specified distance from the interface. The z value of interface is hard coded 
#as well as the distance from the interface and the number of beads in a molecule (easy to read this value 
#in using raw_input) *python script_name input(pdb) output(pdb)
#Writes out new trimmed pdb file 


#/usr/bin/env python

import sys
import numpy as np

#Open input file and output file from command line 
f1 = open(sys.argv[1],'r') #Open in read format
f2 = open(sys.argv[2], 'wa') #Open in write and append format
in_stream = f1 
out_stream = f2


#Initialize variables-----------------------------------------------
iline = 0
atom_name = []
x = []
y = []
z = []
accepted_surf = []


#Read in data from the pdb file-------------------------------------
for line in in_stream.readlines(): 
  if line.startswith("CRYST1"):
    comment = line 
  else:
    data = line.split()
    if data[0] == 'End': 
      lastline = iline 
    if data[0] == 'ATOM': 
      atom_name.append([line[13:16]])
      x.append([line[31:38]]) 
      y.append([line[39:46]])
      z.append([line[47:54]]) 
  iline +=1
in_stream.close()

#   1 - 6 Record name "ATOM " 
#   7 - 11 Integer serial Atom serial number. 
#   13 - 16 Atom name Atom name. 
#   17 Character altLoc Alternate location indicator. 
#   18 - 20 Residue name resName Residue name. 
#   22 Character chainID Chain identifier. 
#   23 - 26 Integer resSeq Residue sequence number. 
#   27 AChar iCode Code for insertion of residues. 
#   31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms. 
#   39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms. 
#   47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms. 
#   55 - 60 Real(6.2) occupancy Occupancy. 
#   61 - 66 Real(6.2) tempFactor Temperature factor. 
#   77 - 78 LString(2) element Element symbol, right-justified. 
#   79 - 80 LString(2) charge Charge on the atom.

#Create arrays, concatenate to form one large array, reshape, and delete----
atom_name_array = np.array(atom_name, int)
x_array = np.array(x, float)
y_array = np.array(y, float) 
z_array = np.array(z, float) 

#concatenate the four [N,1] arrays together to form and [N,4] array 
name_xyz_array = np.concatenate((atom_name_array, x_array, y_array, z_array), axis=1)

remove_w = (name_xyz_array == 3).sum(1)
surf_array = name_xyz_array[remove_w == 0, :] #create an array without water
w_array = name_xyz_array[remove_w != 0, :]    #create an array with only water
#reshape array so that that tail beads in a molecule are all in the same dimension 
mol_array = surf_array.reshape(-1, 10, 4)

#search through the first beads to find  z values inside  of the range and copy
#molecule to a new list 
imol = 0
for zval in  mol_array[:,0,3]:
  if zval<45 and zval>-45: 
    accepted_surf.append(mol_array[imol,:,:])    
  else:
    pass
  imol +=1

accepted_surf_array = np.array(accepted_surf, float) 
good_surf_array = accepted_surf_array.reshape(-1, 4)

#--------------------------------------------------------------------------
#writing out the molecule numbers to concatenate onto the surfactant array
w_mol_num_array = np.array(range(1,4801))

molecule_num = []
mol_num = 4801
mol_max = 4801 + accepted_surf_array.shape[0]
x = 4801

for x in range(mol_num,mol_max):
  for i in range(10):
    molecule_num.append(x)
  i+=1
x+=1
surf_mol_num_array = np.array(molecule_num, int)

arr1 = np.concatenate((w_mol_num_array, surf_mol_num_array),axis=0)
mol_num = arr1[:,np.newaxis]
#--------------------------------------------------------------------------

#Create the pdb array that will be used to write out file -----------------
array = np.concatenate((w_array,good_surf_array),axis=0)  
pdb_array = np.concatenate((array, mol_num), axis=1)

#Writing out pdb file -----------------------------------------------------
s=comment

iatom = 0 
for (i, (ID, x, y, z, mol)) in enumerate(pdb_array):
  s += "ATOM %5i %4s    %3s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i+1, ID, 'X', mol, x, y, z, 0.00, 0.00) 
s+= "End\n"

out_stream.write(s)
out_stream.close()







