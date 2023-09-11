#
# 14/07/2015  Samuel Ponce
#
from __future__ import print_function
import numpy as np

i = 0
for ii in np.arange(0.5, 0.0, -1.0 / 200):
    print('{0} 0.0 0.0 {1}'.format(ii, 1.0 / 201))
    i += 1

for ii in np.arange(0.0, 0.5 + 1.0 / 200, 1.0 / 200):
    print('{0} {0} 0.0 {1}'.format(ii, 1.0 / 201))
    i += 1

print(i)
