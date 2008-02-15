#!/bin/sh

MakeDir() {
    # usage:   MakeDir directory
    # purpose: recursively make a directory only if it does not exists !!!
    # example: MakeDir /a/b/c/d

    dir=${1%/*}
    if [ \( "$dir" != "" \) -a \( "$dir" != "$1" \) ]; then
	MakeDir $dir
    fi
	
    if [ ! -d $1 ]; then
	mkdir $1
    fi
}

if [ $# -ne 1 ]; then
    echo "Usage: MakeDir directory"
    exit 1
fi
MakeDir $1