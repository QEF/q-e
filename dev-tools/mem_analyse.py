from __future__ import print_function
import numpy as np

#1) Run dev-tools/mem_counter inside EPW/src
#2) compile UtilXlib/mem_counter.f90  with -D__DEBUG flag
#3) Run EPW
#4) grep ' allocating' epw1.out > alloc.txt
#5) grep 'deallocating' epw1.out > dealloc.txt
#6) Run this script after having changed the correct allocation lengths

alloc_len = 38817
dealloc_len = 38769

ii = 0
alloc_name = [None] * alloc_len
alloc_size = np.zeros((alloc_len))
alloc_sub = [None] * alloc_len

with open('alloc.txt','r') as R:
  for lines in R:
    tmp = lines.split()
    alloc_name[ii] = str(tmp[4])
    alloc_sub[ii] = str(tmp[5])
    alloc_size[ii] = np.float(tmp[1])
    ii+=1

ii = 0
dealloc_name = [None] * dealloc_len
dealloc_size = np.zeros((dealloc_len))
with open('dealloc.txt','r') as R:
  for lines in R:
    tmp = lines.split()
    dealloc_name[ii] = str(tmp[4])
    dealloc_size[ii] = np.float(tmp[1])
    ii+=1


deall_found = [ False ] * dealloc_len

for ii in np.arange(alloc_len):
  print(ii, ' / ', alloc_len)
  name = alloc_name[ii]

  found = False
  for jj in np.arange(dealloc_len):
    if name == dealloc_name[jj]:
      if alloc_size[ii] == dealloc_size[jj] and not deall_found[jj]:
        # We found the corresponding all/deall pair
        deall_found[jj] = True
        found = True
        break
  if not found:
    with open('mem_analyse.out','a') as O:
      O.write('We did not find a maching pair in '+str(alloc_sub[ii])+'\n')
      O.write('Allocate:   '+str(name)+' '+str(alloc_size[ii])+'\n')
      O.write('Deallocate: '+str(dealloc_name[jj])+' '+str(dealloc_size[jj])+ '\n')
#    print 'We did not find a maching pair in ', alloc_sub[ii]
#    print 'Allocate:   ',name,' ',alloc_size[ii]
#    print 'Deallocate: ',dealloc_name[jj],' ',dealloc_size[jj]
