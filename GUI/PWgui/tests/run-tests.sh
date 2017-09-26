#!/bin/bash

module_message() {
    # Usage: $0 module
    echo "
 ==============================================================================
 
             * * *  Running GUI($1) tests. Please wait * * *

 ==============================================================================
"
}
file_message() {
    # usage: $0 file
    echo "
--
  Testing file: $1
--
"
}
run() {
    # usage: $0 dir module
    module_message $2
    here=`pwd`
    pwgui_dir=$here/..
    cd $1
    files=`ls *.{in,inp} neb.dat 2> /dev/null`
    for file in $files; do
	file_message $file
	$pwgui_dir/pwgui_reformat $2 $file
	if test $? != 0 ; then
	    echo "an error occurred for module: $2 , file: $file"
	    exit 1
	fi
    done
    cd $here
}

for dir in atomic  neb.dat  ph  pp  projwfc  pw
do
    module=${dir%.dat}; # neb.dat --> neb
    run ../examples/$dir $module
done
