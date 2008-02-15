#!/bin/sh

if [ $# -ne 3 ]; then
    echo "Usage: $0 NAME program_name txt_docfile"
    exit 1
fi

docfile=$3

#version=$(ls $docfile-*)
#version=${version#$docfile-}

version=`grep def users-guide.tex | grep '\\\def\\\version' | awk '{split($1,v,"[{}]"); print v[2]}'`

echo "
<html>
<body>
  <H1> Description of $1 ($2) Input File</H1>

    <h2> PWscf version:</b> $version</h2> 
    <pre>" > head.$$

echo "
    </pre>
  </body>
</html>" > tail.$$

sed 's/</\&lt\;/g' $docfile | sed "s/>/\&gt\;/g" - > body.$$

cat head.$$ body.$$ tail.$$
rm -f head.$$ body.$$ tail.$$



