#
# Script to automatically add status and error message to ALLOCATE and DEALLOCATE in
# Fortran files.
#
# Author: S. Ponce
# Date: Sept. 2019
#

from __future__ import print_function
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
      print('{}ALLOCATE({}, STAT = ierr)'.format(tmp[0:start-9], tmp[start:final]))
      print("{}IF (ierr /= 0) CALL errore('{}', 'Error allocating {}', 1)".format(
          tmp[0:start-9], name_sub, tmp[start:end]
      ))
    elif tmp_split[0][0:11] == 'DEALLOCATE(':
      for ii in np.arange(len(tmp)):
        if tmp[ii] == '(':
          start = ii
        if tmp[ii] == ')':
          end = ii
      print('{}DEALLOCATE({}, STAT = ierr)'.format(tmp[0:start-10], tmp[start+1:end]))
      print("{}IF (ierr /= 0) CALL errore('{}', 'Error deallocating {}', 1)".format(
          tmp[0:start-10], name_sub, tmp[start+1:end]))
    else:
      print(lines, end='')
