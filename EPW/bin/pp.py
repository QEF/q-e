#
# Post-processing script from of PH data in format used by EPW
# 14/07/2015 - Creation of the script - Samuel Ponce
# 14/03/2018 - Automatically reads the number of q-points - Michael Waters
# 14/03/2018 - Detect if SOC is included in the calculation - Samuel Ponce 
# 
import numpy as np
import os
from xml.dom import minidom

# Return the number of q-points in the IBZ
def get_nqpt(prefix):
  fname = '_ph0/' +prefix+'.phsave/control_ph.xml'

  fid = open(fname,'r')
  lines = fid.readlines() # these files are relatively small so reading the whole thing shouldn't be an issue
  fid.close()

  line_number_of_nqpt = 0
  while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]: # increment to line of interest
    line_number_of_nqpt +=1
  line_number_of_nqpt +=1 # its on the next line after that text

  nqpt = int(lines[line_number_of_nqpt])

  return nqpt

# Check if the calculation include SOC
def hasSOC(prefix):
  fname = prefix+'.save/data-file-schema.xml'

  xmldoc = minidom.parse(fname)
  item = xmldoc.getElementsByTagName('spinorbit')[0]
  lSOC = item.childNodes[0].data
  
  return lSOC

# Check if the calculation was done in sequential
def isSEQ(prefix):
  fname = '_ph0/'+str(prefix)+'.dvscf'
  if (os.path.isfile(fname)):
    lseq = True
  else:
    lseq = False
 
  return lseq
    
# Enter the number of irr. q-points
user_input = raw_input('Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = str(user_input)

# Test if SOC
SOC = hasSOC(prefix)

# Test if seq. or parallel run
SEQ = isSEQ(prefix)

if True: # this gets the nqpt from the outputfiles
  nqpt =  get_nqpt(prefix)

else:
  # Enter the number of irr. q-points
  user_input = raw_input('Enter the number of irreducible q-points\n')
  nqpt = user_input
  try:
    nqpt = int(user_input)
  except ValueError:
    raise Exception('The value you enter is not an integer!')

os.system('mkdir save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  # Case calculation in seq.
  if SEQ:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
 
  else:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
