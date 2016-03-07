#
# 14/07/2015  Samuel Ponce
#
import numpy as np

i=0
for ii in np.arange(0.5,0.0,-1.0/200):
  print str(ii)+' 0.0 0.0 '+str(1.0/201)
  i +=1

for ii in np.arange(0.0,0.5+1.0/200,1.0/200):
  print str(ii)+' '+str(ii)+' 0.0 '+str(1.0/201)
  i +=1

print i
