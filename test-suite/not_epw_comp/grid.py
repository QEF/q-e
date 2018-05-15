import numpy as np

tot = 0  

for ii in np.arange(0,1,0.25):
  for jj in np.arange(0,1,0.25):
    for kk in np.arange(0,1,0.25):
      print ii,' ',jj,' ',kk,' ',1.0/64
      tot += 1

print tot
