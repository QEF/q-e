#!/bin/sh

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
    cd $1
    for file in `ls *.inp`; do
	file_message $file
	$PWGUI/pwgui_reformat $2 $file
	if test $? != 0 ; then
	    echo "an error occurred for module: $2 , file: $file"
	    exit 1
	fi
    done
    cd $here
}

#for module in \
#    pw \
#    ph \
#    pp \
#    projwfc \
#    chdens \
#    d3
#  do
#  run ../examples/$module $module
#done

run ../examples/pw pw
run ../examples/ph ph
run ../examples/pp pp
run ../examples/projwfc projwfc
run ../examples/chdens chdens
run ../examples/d3 d3