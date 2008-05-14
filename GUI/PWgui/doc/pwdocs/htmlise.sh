#!/bin/sh

if [ $# -ne 3 ]; then
    echo "Usage: $0 NAME program_name txt_docfile"
    exit 1
fi

docfile=$3

#version=$(ls $docfile-*)
#version=${version#$docfile-}

#version=`grep def users-guide.tex | grep '\\\def\\\version' | awk '{split($1,v,"[{}]"); print v[2]}'`
version=4.0

echo "
<html>
  <body>
    <pre>" > head.$$

echo "
    </pre>
  </body>
</html>" > tail.$$

# filter out line: "*** FILE AUTOMATICALLY CREATED ..."
awk '$0 !~ /FILE AUTOMATICALLY CREATED: DO NOT EDIT/ { print; }' $docfile > body.$$
sed 's/</\&lt\;/g' body.$$ | sed "s/>/\&gt\;/g" - > body2.$$

cat head.$$ body2.$$ tail.$$
rm -f head.$$ body.$$ body2.$$ tail.$$



