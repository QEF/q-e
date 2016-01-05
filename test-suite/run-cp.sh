#!/bin/bash
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainer: Filippo Spiga (filippo.spiga@quantum-espresso.org)

include ${ESPRESSO_ROOT}/test-suite/ENVIRONMENT

# -- Uncomment to run in parallel
#export PARA_PREFIX="mpirun -np 4"
unset PARA_PREFIX

${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/cp.x "$@"

rm -f input_tmp.in
