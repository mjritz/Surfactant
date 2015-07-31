#Reads in a DCD file and writes xyz format out for each frame.

#!/usr/bin/env python 

import struct
import numpy as np

def read_DCD(filename):
    #warning no input checking is performed
    #assumes dcd is well behaved and of expected form
    D=3
    #open file filename for reading DCD
    f=open(filename,"rb")
    #read header information in dcd file
    header=struct.unpack('i4s9if13i80s80s4i',f.read(struct.calcsize('i4s9if13i160c4i')))
    nAt=header[-2]
    nFrames=header[2]
    xyz_frames=np.zeros((nFrames,nAt,D))
    BoxD_frames=np.zeros((nFrames,3))
    BoxAngle_frames=np.zeros((nFrames,3))

    for i in range(nFrames):
        f.read(4)
        frame=struct.unpack('6d4x'+3*('4x'+str(nAt)+'f4x'),f.read(struct.calcsize('6d4x'+3*('4x'+str(nAt)+'f4x'))))
        BoxD_frames[i]=np.array((frame[0],frame[2],frame[5]),float)
        BoxAngle_frames[i]=np.degrees(np.array((np.arccos(frame[4]),np.arccos(frame[3]),np.arccos(frame[1]))))
        xyz_frames[i]=np.array(frame[6:],float).reshape(3,nAt).T
    f.close()

    return nAt,BoxD_frames,BoxAngle_frames,xyz_frames

    
