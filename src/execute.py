#!/usr/bin/env python

#This script reads in a surfactant structure file and outputs a data file for
#LAMMPS: 'python execute.py surfactant_structure.py' 
#In the same directory must have surfactant structure, wat.lt, mie_surf.lt, 
#water_60xy80z.xyz, pdb_create.py, xyz_coord_bonds.py 


#Assumptions: no dendrimerization, no extra lines at the end of the structure file, T (main tail), B (side branches), H (head groups)
#Name of water file and tolerance are written in pdb_create.py near the bottom
#Initial distance of beads written in xyz_coord_bonds.py
#Size of system is written near the bottom of pdb_create.py

import write_lammps_input
import xyz_coord_bonds
import os
from pdb_create import i
import sys
from xyz_coord_bonds import head_groups
from xyz_coord_bonds import main_tail
from xyz_coord_bonds import descript

version = 'auto1.4'
description = descript 
length = main_tail
head = head_groups

os.system('cp input_file.txt src/.' )

if description != 'mb':
    print "Description of branched surfactant: inc_length, dec_length, eq_length, random"
    further_description = str(raw_input('>'))
    
print "Date (dd_mm_yy)?"
lammps_input_files = str(raw_input('> '))

os.system('/gpfs_share/santiso/SOFTWARE/packmol/packmol < surf_water.inp') 

os.system('/gpfs_share/santiso/SOFTWARE/moltemp_src.2013-3-03/moltemplate.sh  -nocheck -xyz surf_water.xyz  system.lt')

if description == 'sb':
    os.system('mkdir ../%s' %(description))
    os.system('mkdir ../%s/sb%i%i' %(description, head, length))
    os.system('mkdir ../%s/sb%i%i/%s' %(description, head, length,further_description))
    os.system('mkdir ../%s/sb%i%i/%s/%i' %(description, head, length,further_description, i))
    os.system('mkdir ../%s/sb%i%i/%s/%i/%s' %(description, head, length,further_description, i, lammps_input_files))
    os.system('mkdir ../%s/sb%i%i/%s/%i/%s/trim_system' %(description, head, length,further_description, i, lammps_input_files))
    os.system('cp system.data ../%s/sb%i%i/%s/%i/%s/.' %(description,head, length, further_description, i,lammps_input_files))
    os.system('cp %s ../%s/sb%i%i/%s/%i/%s/.' %(sys.argv[1], description,head, length, further_description, i,lammps_input_files))
    os.system('cp %s ../%s/sb%i%i/%s/%i/%s/.' %('system_all.in', description,head, length, further_description, i,lammps_input_files))
    os.system('cp %s ../%s/sb%i%i/%s/%i/%s/.' %('surf.lt', description,head, length, further_description, i,lammps_input_files))
    os.system('cp src/copied_files/* ../%s/sb%i%i/%s/%i/%s/.' %(description,head, length, further_description, i,           lammps_input_files))
    os.system('cp input_file.txt ../%s/sb%i%i/%s/%i/%s/.' %(description,head, length, further_description, i, lammps_input_files))
    os.system('cp surf.xyz ../%s/sb%i%i/%s/%i/%s/.' %(description,head, length, further_description, i, lammps_input_files))
    os.system('rm -rf output_ttree')
    os.system('rm system.i*')
    os.system('clear')
    os.chdir('../%s/sb%i%i/%s/%i/%s/' % (description,head, length, further_description, i,lammps_input_files))
    os.system('bsub < run_lammps')

if description == 'mb': 
    os.system('mkdir ../%s' %(description))
    os.system('mkdir ../%s/mb%i%i' %(description, head, length))
    os.system('mkdir ../%s/mb%i%i/%i' %(description, head, length, i))
    os.system('mkdir ../%s/mb%i%i/%i/%s' %(description, head, length, i,lammps_input_files))
    os.system('mkdir ../%s/mb%i%i/%i/%s/trim_system' %(description, head, length, i,lammps_input_files))
    os.system('cp system.data ../%s/mb%i%i/%i/%s/.' %(description, head, length, i,lammps_input_files))
    os.system('cp %s ../%s/mb%i%i/%i/%s/.' %(sys.argv[1], description, head, length, i,lammps_input_files))
    os.system('cp %s ../%s/mb%i%i/%i/%s/.' %('system_all.in', description, head, length, i,lammps_input_files))
    os.system('cp src/copied_files/* ../%s/mb%i%i/%i/%s/.' %(description, head, length, i,lammps_input_files))
    os.system('cp input_file.txt ../%s/mb%i%i/%i/%s/.' %(description, head, length, i,lammps_input_files))
    os.system('cp surf.xyz ../%s/mb%i%i/%i/%s/.'    %(description, head, length, i,lammps_input_files))
    os.system('cp %s ../%s/mb%i%i/%i/%s/.' %('surf.lt', description, head, length, i,lammps_input_files))
    os.system('rm -rf output_ttree')
    os.system('rm system.i*')
    os.system('clear')
    os.chdir('../%s/mb%i%i/%i/%s/' %(description, head, length, i, lammps_input_files))
    os.system('bsub < run_lammps')

    


