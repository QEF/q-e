#!/bin/bash
aclocal -I m4 --install
autoconf -f -v --output=configure-new configure.ac.new
