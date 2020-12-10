#
# Post-processing script QE --> EPW
# 14/07/2015 - Samuel Ponce
#

from builtins import input
import numpy as np
import os

# Enter the number of irr. q-points
prefix = input('Enter the prefix used for PH calculations (e.g. diam)\n')

# Enter the number of irr. q-points
nqpt = input('Enter the number of irreducible q-points\n')
try:
  nqpt = int(nqpt)
except ValueError:
  raise Exception('The value you enter is not an integer!')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
  if (iqpt == 1):
    os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r _ph0/'+prefix+'.phsave save/')
  else:
    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('rm -f _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
