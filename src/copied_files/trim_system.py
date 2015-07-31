#!/usr/bin/env python

import  def_surface
import  variables
import sys

#print "Input file for moltemplate"
in_lt_file = open(sys.argv[2], 'w')

xbox, ybox, zbox, w_filename, w_mol, pack_tol, bead_dist, boxx, boxy, boxz,     bead_per, main_branch, equil_time, surf_x, surf_y, surf_zlow, surf_zhi,         per_loop, steps,temp, z_distance, epsilon, sigma, r_distance, equil_dcd,        prod_dcd, cutoff = variables.read_variables('input_file.txt')

g = 'import "wat.lt" \n'
g+= 'import "surf.lt" \n'

g+='''
write_once("Data Boundary") {
  -%4.6f   %12.6f xlo xhi
    -%4.6f   %12.6f ylo yhi
    -%4.5f   %12.5f zlo zhi
} \n\n''' % (boxx, boxx, boxy, boxy, boxz, boxz)

g+= "wat = new Wat [%i] \n\n" % (w_mol)

g+= "surf = new Surf [%i]" % (int(def_surface.num_surf_mol))

in_lt_file.write(g)
in_lt_file.close()
