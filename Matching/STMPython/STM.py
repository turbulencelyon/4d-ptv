#! /usr/bin/env python
import time
import os
import socket
import sys

print(socket.gethostname())
print("Python", sys.version)
print("curdir =", os.path.realpath(os.curdir))
sys.path.insert(0, "/home/eberna07/Stage_EB_2020/4d-ptv/Matching/STMPython")
print(sys.path)
import STMFunctions as stmf
import math
import numpy as np
import itertools as it
import copy
import struct
from datetime import datetime 


def STM(filename,minframes,maxframes,cammatch,maxdistance,nx,ny,nz,maxmatchesperray, boundingbox=[[-140, 140], [-150, 150], [5,170]], neighbours=6):
    """
    Compute matches from rays projecting them into voxels.
    
    # Parameters
    #   filename                           # name of the file containing rays
    #   minframes                          # number of the first frame
    #   maxframes                          # number of the last frame
    #   cammatch                           # minimum number of rays crossing to get a match. We require a match to satisfy cammatchfunc = lambda x: len(x)>=cammatch for the number of cameras (and thus the number of rays)
    #   maxdistance                        # max distance allowed for a match.
    #   nx,ny,nz                           # number of voxels in each direction 
    #   maxmatchesperray                   # number of matches/ray 
    #   boundingbox                        # correspond to the volume visualized [[minX,maxX],[minY,maxY],[minZ,maxZ]] ATTENTION Does not work currently -> To DO !!
    #   neighbours                         # number of illuminated voxels: due to noise, when a ray crosses a voxel, it is possible that in reality, the ray crosses a close voxel. neighbours indicates how many neighbours we consider in reality when a ray crosses a voxel. =6 by defaut.
    """
    #############################################################################################################
    # Parameters to adjust
    tstart = datetime.now().timestamp()
    cammatchfunc = lambda x: len(x)>=cammatch       # We require a match to satisfy this requirement for the number of cameras (and thus the number of rays)
    #############################################################################################################
    
    fileout = copy.copy(filename).split(".")
    fileout = ".".join(fileout[0:len(fileout)-1])
    filelog = fileout + ".log"
    print(filelog)
    fileout = fileout.replace("rays","matched")+"cam{}_{}-{}.dat".format(cammatch,minframes,maxframes)

    fout = open(fileout, 'wb')
    fin = open(filename, 'rb')
    frameid = minframes
    numpts = fin.read(4)                                                 # Read 4 bytes header
    while(len(numpts)>0 and frameid < maxframes):                        # If something is read
        numpts = struct.unpack('I', numpts)[0]                           # Interpret header as 4 byte uint
        flog = open(filelog, 'a')
        flog.write("#######\n")
        flog.write("Frame: " + str(frameid) + "\nNumber of rays: " + str(numpts) + "\n")
        flog.close()
        
        print("Frame:",frameid,". # of rays:", numpts)
        
        # Read rays
        raydata = fin.read(numpts*27)                                     # 27 bytes per line 2+1+6*4
        raydata = struct.unpack('='+('BH6f'*numpts),raydata)              # Create string '=BHFFFFFFBHFFFFFFBHFFFFFF...BHFFFFFF'
        raydata = list(map(lambda i: list(raydata[8*i:8*i+8]) ,range(len(raydata)//8)))     # Reshape to 8*N np.arreyreshape converts everything to floats...
        # The actual call
        output = stmf.SpaceTraversalMatching(list(raydata),boundingbox,nx=nx,nz=nz,ny=ny,cammatchfunc=cammatchfunc,neighbours=neighbours,logfile=filelog,maxdistance=maxdistance)
        
        # Prepare output
        print("Matches found:",len(output))
        coutput = []
        for o in output:
            tmp = [len(o[0])]
            tmp.extend(o[1])
            tmp.append(o[2])
            for j in o[0]:
                tmp.extend(j)
            coutput.append(tmp)
        output = coutput
        del coutput
        
        # Output the output
        buf = struct.pack('I', len(output))                               # Write the number of matches
        fout.write(buf)
        
        for out in output:
            buf = struct.pack('=B4f' + out[0]*"BH", *out)                 # Write each to the output file
            fout.write(buf)

        numpts = fin.read(4)                                              # Read next header
        #print("numpts:",numpts)
        #print("type:",type(numpts))
        frameid += 1
    fout.close()
    fin.close()
    print("Finished")

    elapsed = datetime.now().timestamp() - tstart
    print("Elapsed time:",elapsed)
    print("Elapsed time/frame:",elapsed/frameid)

if(len(sys.argv)==10):
    STM(str(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]))
elif (len(sys.argv)==11):
    STM(str(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]),np.array(sys.argv[10]))
else:
    print("Only {} arguments".format(len(sys.argv)-1))
    print("There should be an argument with the filename, minframe, maxframe, cammatch, maxdistance, nx, ny, nz, maxmatchesperray, boundingbox(optional)!")
    
