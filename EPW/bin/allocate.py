# 
# Script to automatically add status and error message to ALLOCATE and DEALLOCATE in 
# Fortran files. 
#
# Author: S. Ponce
# Date: Sept. 2019
#

import numpy as np

# File name
file_name = 'elphel2_shuffle.f90'

with open(file_name,'r') as F:
  for lines in F:
    tmp = lines
    tmp_split = lines.split()
    #tmp = lines.split()
    if len(tmp_split) < 1:
      continue
    # 
    #print tmp_split[0][0:9]
    if tmp_split[0] == 'SUBROUTINE':
      for ii in np.arange(len(tmp)):
        if tmp[ii] == 'E' and tmp[ii+1] == ' ':
          start = ii
        if tmp[ii] == '(':
          end = ii

      name_sub = str(tmp[start+2:end])
    
    if tmp_split[0] == 'FUNCTION':
      for ii in np.arange(len(tmp)):
        if tmp[ii] == 'N' and tmp[ii+1] == ' ':
          start = ii
        if tmp[ii] == '(':
          end = ii

      name_sub = str(tmp[start+2:end])

    if tmp_split[0][0:9] == 'ALLOCATE(':
      for ii in np.arange(len(tmp)):
        if tmp[ii] == 'E':
          start = ii + 2
        if tmp[ii] == '(':
          end = ii
        if tmp[ii] == ')':
          final = ii
      #print 'start ',start
      print str(tmp[0:start-9])+'ALLOCATE('+str(tmp[start:final])+', STAT = ierr)'
      print str(tmp[0:start-9])+'IF (ierr /= 0) CALL errore(\''+str(name_sub)+'\', \'Error allocating '+str(tmp[start:end])+'\', 1)'     
    elif tmp_split[0][0:11] == 'DEALLOCATE(':
      for ii in np.arange(len(tmp)):
        if tmp[ii] == '(':
          start = ii
        if tmp[ii] == ')':
          end = ii
      print str(tmp[0:start-10])+'DEALLOCATE('+str(tmp[start+1:end])+', STAT = ierr)' 
      print str(tmp[0:start-10])+'IF (ierr /= 0) CALL errore(\''+str(name_sub)+'\', \'Error deallocating '+str(tmp[start+1:end])+'\', 1)' 
    else:
      print str(lines),




