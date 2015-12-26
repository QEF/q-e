#!/bin/bash
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

source ${ESPRESSO_ROOT}/test-suite/ENVIRONMENT

if test "`which curl`" = "" ; then
   if test "`which wget`" = "" ; then
      echo "### wget or curl not found: will not be able to download missing PP ###"
   else
      DOWNLOADER="wget -O"
      # echo "wget found"
   fi
else
   DOWNLOADER="curl -o"
   # echo "curl found"
fi

inputs=`find $1* -type f -name "*.in" -not -name "test.*" -not -name "benchmark.*"`
pp_files=`for x in ${inputs}; do grep UPF ${x} | awk '{print $3}'; done`

for pp_file in ${pp_files} ; do
if ! test -f ${ESPRESSO_PSEUDO}/${pp_file} ; then 
	#echo -n "Downloading ${pp_file} to ${ESPRESSO_PSEUDO} ... "
	${DOWNLOADER} ${ESPRESSO_PSEUDO}/${pp_file} ${NETWORK_PSEUDO}/${pp_file} 2> /dev/null
	if test $? != 0 ; then 
		echo "FAILED, do it manually -- Testing aborted!"
		exit -1
	#else
		#echo "done."
	fi
#else
	#echo "No need to download ${pp_file}."
fi
done
