#!/usr/bin/env python  

import os  
import trim_system

os.system('/gpfs_share/santiso/mjritz/LAMMPS/Surf/automated/moltemp_src/moltemplate.sh -nocheck -xyz trimmed_system.xyz system_trimmed.lt')    

os.system('rm -rf output_ttree')
os.system('rm system_trimmed.i*')

print "\n\nlower interface:",trim_system.def_surface.lower_interface,";", "upper interface:",trim_system.def_surface.upper_interface
print "accepted surfactants", trim_system.def_surface.acceptable_range,"angstroms away from lower and upper interface"
print  "started with:",trim_system.def_surface.surf_molecules,"surfactants;", "deleted:",trim_system.def_surface.surf_molecules - trim_system.def_surface.good_surf_array.shape[0]/trim_system.def_surface.bead_per_surfactant,"surfactants\n\n"
