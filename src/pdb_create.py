#Reads in structure of surfactant and outputs xyz file, bond information, and
#modifies system.lt file for running moltemplate
#Assumption: No dendrimerization and surfactants are symmetric


#!/usr/bin/env python 

import sys
import numpy as np
import variables 

#Read in the surfactant structure and input file for moltemplate 
in_structure = open(sys.argv[1], 'r') 
out_structure = open('surf.xyz', 'wa')
packmol_file = open('surf_water.inp', 'wa')

#
with in_structure as f:                                                                                  
    lines = f.readlines()                                                                                
    num_lines = len([l for l in lines if l.strip(' \n') != ''])   

xbox, ybox, zbox, w_filename, w_mol, pack_tol, bead_dist, boxx, boxy, boxz,     bead_per, main_branch, equil_time, surf_x, surf_y, surf_zlow, surf_zhi,         per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd,        prod_dcd, cutoff = variables.read_variables('input_file.txt')

#Define initial variables and lists
number_of_atoms = 0
i = 0 
main_branch= str()
branch = str()
j = 0
k = 0
zero = 0.00000

#Reading structure file and writing out structure as xyz file
#structure file cannot have any comment lines (only the structure)
in_structure_again = open(sys.argv[1], 'r')
for lines in in_structure_again.readlines():
  if 'B' in lines or 'T' in lines or 'H' in lines:
    atoms_on_each_line = len(lines.split())
    number_of_atoms += atoms_on_each_line
  for x in lines.split():                                                                                 
    if 'T'in x:                                                                                         
      main_branch += "CM         %9.5f         %9.5f          %9.5f\n" %(zero,zero,i)                                     
      i+=float(bead_dist)       
    if x == 'H':
      main_branch += "OA         %9.5f         %9.5f          %9.5f\n" %(zero,zero,i)  
      i+= float(bead_dist)
  if lines.count('B')==0:
    pass
  if lines.count('B') > 0:
    number_side_beads = len(lines.split())-1
    z_value =float(bead_dist)*(j) 
    for bead in lines.split():
      y_value = float(bead_dist)*-(number_side_beads/2)+k
      if y_value != 0:
        branch += "CM         %9.5f         %9.5f          %9.5f\n" % (zero, y_value ,z_value)   
      k+=float(bead_dist)
    k=0
  j +=1


s = str(number_of_atoms)                                                                                 
s += "\nSU1.xyz\n" 
s += main_branch
s += branch

#print "Input file for moltemplate"
in_lt_file = open('system.lt', 'w')

print "Number of Surfactants on Interface:"
i = int(raw_input('> '))

g = 'import "wat.lt" \n'
g+= 'import "surf.lt" \n'

g+='''
write_once("Data Boundary") {
-%4.6f   %12.6f xlo xhi
-%4.6f   %12.6f ylo yhi
-%4.5f   %12.5f zlo zhi
} \n\n''' % (float(boxx), float(boxx), float(boxy), float(boxy), float(boxz), float(boxz))

g+= "wat = new Wat [%i] \n\n" % (w_mol)

g+= "surf = new Surf [%i]" % (i*2)

in_lt_file.write(g)
in_lt_file.close()
out_structure.write(s)
out_structure.close()

#writing out packmol input file 
h = "tolerance %f\n" % (pack_tol)                                                                                           
h+="filetype xyz\n"                                                                                             
h+="output surf_water.xyz\n\n"""                                                                                                         

h += "structure %s\n" % (w_filename)                                                                              
h += "  fixed 0. 0. 0. 0. 0. 0.\n"                                                                              
h += "end structure\n\n"                                                                                                          

h += "structure surf.xyz\n"                                                                                       
h += "  resnumber 3\n"

h += "  number %i\n" % (i)                                                                                             
h += "  inside box -%4.1f -%4.1f %4.1f %4.1f %4.1f %4.1f \n" % (surf_x, surf_y, surf_zlow, surf_x, surf_y, surf_zhi)                                                             
h += "  end atoms\n"  
h += "end structure\n\n"                                                                                            
                                                                                                                   
h += "structure surf.xyz\n"                                                                                       
h += "  resnumber 3\n"                                                                                            
h += "  number %i\n" %(i)                                                                                             
h += "  inside box -%4.1f -%4.1f -%4.1f %4.1f %4.1f -%4.1f \n" % (surf_x, surf_y, surf_zhi, surf_x, surf_y, surf_zlow)                                                          
h += "  end atoms\n"                                                                                              
h += "end structure\n\n"                                                                                            

packmol_file.write(h)
packmol_file.close()
