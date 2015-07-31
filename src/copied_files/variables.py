#!/usr/bin/env python 

import struct
import numpy as np

def read_variables(filename):
    var_file=open(filename,"rb")
    for vars in var_file.readlines():
        if '#' not in vars:
            variable = vars.split(' ')
            if variable[0] == 'water_box':
                xbox = float(variable[1])
                ybox = float(variable[2])
                zbox = float(variable[3])
            if variable[0] == 'water_filename':
                w_filename = variable[1]
            if variable[0] == 'water_molecules':
                w_mol = int(variable[1])
            if variable[0] == 'packmol_tolerance':
                pack_tol = float(variable[1])
            if variable[0] == 'bead_distance':
                bead_dist = variable[1]
            if variable[0] == 'box_size':
                boxx = float(variable[1])
                boxy = float(variable[2])
                boxz = float(variable[3])
            if variable[0] == 'bead_per_molecule':
                bead_per = variable[1]
            if variable[0] == "bead_in_main_branch":
                main_branch = variable[1]
            if variable[0] == 'equilibrate':
                equil_time = variable[1]
            if variable[0] == 'surf_space':
                surf_x = float(variable[1])
                surf_y = float(variable[2])
                surf_zlow = float(variable[3])
                surf_zhi = float(variable[4])
            if variable[0] == 'num_loops':
                per_loop = variable[1]
            if variable[0] == 'step_per_loop':
                steps = variable[1]
            if variable[0] == 'temperature':
                temp = float(variable[1])
            if variable[0] == 'wall':
                z_distance = variable[1]
                epsilon = variable[2]
                sigma = variable[3]
                r_distance = variable[4]
            if variable[0] == 'DCD_filename':
                equil_dcd = variable[1]
                prod_dcd = variable[2]
            if variable[0] == 'cutoff_distance':
                cutoff = variable[1]
    return xbox, ybox, zbox, w_filename, w_mol, pack_tol, bead_dist, boxx, boxy, boxz, bead_per, main_branch, equil_time, surf_x, surf_y, surf_zlow, surf_zhi, per_loop, steps, temp, z_distance, epsilon, sigma, r_distance, equil_dcd, prod_dcd, cutoff
