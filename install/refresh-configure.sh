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
# Dependency: GNU AUTOCONF (minimum 2.60, 2.69 suggested)
#
# Tested using:
# - autoconf (GNU Autoconf) 2.63
# - aclocal (GNU automake) 1.11.1
# - m4 (GNU M4) 1.4.13

aclocal -I m4 --install

autoconf -f -v --output=configure configure.ac
