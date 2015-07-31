#!/usr/bin/env python  

import os  
import trim_system
import sys

f1 = open(sys.argv[3],'a') #Open in read format                                                          
out_stream = f1

lt_file = sys.argv[2]

os.system('/gpfs_share/santiso/SOFTWARE/moltemp_src.2013-3-03/moltemplate.sh -nocheck -xyz               trimmed_system.xyz %s'%(lt_file))

os.system('rm -rf output_ttree')
os.system('rm system_trimmed.i*')

#Writing an output file 
p = "-------------------------------------------------------------------"
p += "\nlower interface: %5.3f; upper interface: %5.3f\n"  %(trim_system.def_surface.lower_interface ,trim_system.def_surface.upper_interface)
p += "accepted surfactants %5.3f angstroms away from lower and upper interface\n" %(trim_system.def_surface.acceptable_range)
p += "started with: %i surfactants; deleted: %i surfactants\n"  % (trim_system.def_surface.surf_molecules, trim_system.def_surface.deleted)
p+= "top layer has %i surfactants and bottom layer has %i surfactants\n" %(trim_system.def_surface.total_top, trim_system.def_surface.total_bottom)
p+= "-------------------------------------------------------------------\n\n"

out_stream.write(p)
out_stream.close()

