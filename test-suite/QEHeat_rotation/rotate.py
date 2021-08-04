import numpy as np
import sys
from math import *
a=1.37*pi/4.0
b=-0.77*pi
c=0.23*pi
r1=np.array([[1,0,0],[0,cos(a),-sin(a)],[0,sin(a),cos(a)]])
r2=np.array([[cos(b),0,sin(b)],[0,1,0],[-sin(b),0,cos(b)]])
r3=np.array([[cos(c),-sin(c),0],[sin(c),cos(c),0],[0,0,1]])
rm=r1@r2@r3 #rotation matrix

if len(sys.argv)==2:
    np.savetxt(sys.argv[1],rm)
    sys.exit()
for infile, rotfile in zip(sys.argv[1::2], sys.argv[2::2]):
    pos=np.loadtxt(infile)
    rp=np.matmul(pos,rm)
    np.savetxt(rotfile,rp)
