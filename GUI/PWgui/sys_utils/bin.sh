#!/bin/sh

cat $1 | awk '
/^exec/ { print "exec $PWGUI/bin/itkwish \"$0\""; next; }
/a*/    { print; }' > $1.bin