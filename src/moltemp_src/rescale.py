import numpy as np
import sys
import time
#This code is not robust and a
#this script assumes the following:
# 1 type of molecule
# molecules are in a continuous block
# assumes box center at 0,0,0
# assumes reading unwrapped coordinates

def mol_trans(xyz,delta,mol_com):
# translate com of molecules
	xyz_s=np.empty_like(xyz)
	scale=np.array([np.sqrt(1.0+delta),np.sqrt(1.0+delta),1.0/(1.0+delta)])
	#print scale
	
	#print xyz[0]
	#print mol_com[0]
	transcom=(mol_com*scale)-mol_com
	#print transcom[0]
	for i in range(np.size(xyz,axis=0)):
		xyz_s[i]=xyz[i]+transcom[i]
	#print xyz_s[0]
	return np.array(xyz_s.reshape((-1,3)))

def write_lammps_dump(xyz,lammps_dump,header,dump_file_name):
	s=header
	for i,(pos) in enumerate(xyz):
		s+= ('%d %20.15g %20.15g %20.15g\n') % (lammps_dump[i,0],pos[0],pos[1],pos[2])
	f=open(dump_file_name,'w')
	f.write(s)
	f.close()

def main():
	D=3
	boxD=np.zeros(D)
	
	#set variables from command line inputs
	dump_file=sys.argv[1]
	xyz_p_file=sys.argv[2]
	xyz_m_file=sys.argv[3]
	delta=float(sys.argv[4])
	mol_size=int(sys.argv[5])
	boxD[0]=float(sys.argv[6])
	boxD[1]=float(sys.argv[7])
	boxD[2]=float(sys.argv[8])
	
	#read lammps dump file
	f=open(dump_file,'r')
	dum=f.readlines()
	lfile=[]
	for i in dum:
		lfile.append(i.split())
	lammps_dump=np.array(lfile[9:],float)
	
	#make lammps header for output. Everything except number of atoms and ITEM ATOMS has arbitrary values
	lammps_dump_header="ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"+str(np.size(lammps_dump,axis=0))+"\nITEM: BOX BOUNDS pp pp pp\n\
	-40 40\n-40 40\n-150 150\nITEM: ATOMS id xu yu zu"

	xyz=lammps_dump[:,2:5]
	mass=lammps_dump[:,1]

	#compute center of mass of each molecule
	mxyz=np.zeros_like(xyz)
	mol_mass=np.sum(mass[0:mol_size])
	for k in range(3):
		mxyz[:,k]=xyz[:,k]*mass
	mxyz=mxyz.reshape((-1,mol_size,D))
	mol_com=np.sum(mxyz,axis=1)/(mol_mass)
	#wrap center of mass to minimum image
	#this should not be necessary with unwrapped coordinates
	#mol_com-=boxD*np.around(mol_com/boxD)


	# translate com of molecules for plus delta and minus delta shift
	#create a nmolecules x mol_size x D array
	xyz=xyz.reshape((-1,mol_size,D))

	xyz_p=mol_trans(xyz,delta,mol_com)
	xyz_m=mol_trans(xyz,-delta,mol_com)
	
	
	
	write_lammps_dump(xyz_p,lammps_dump,lammps_dump_header,xyz_p_file)
	write_lammps_dump(xyz_m,lammps_dump,lammps_dump_header,xyz_m_file)
	return 0

a=main()

