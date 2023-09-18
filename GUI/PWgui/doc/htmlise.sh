#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Usage: $0 txt_docfile"
    exit 1
fi

docfile=$1

echo "
<html>
  <body>
    <pre>" > head.$$

echo "
    </pre>
  </body>
</html>" > tail.$$

# filter out line: "*** FILE AUTOMATICALLY CREATED ..."
awk '$0 !~ /FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST/ { print; }' $docfile > body.$$
sed 's/</\&lt\;/g' body.$$ | sed "s/>/\&gt\;/g" - > body2.$$

cat head.$$ body2.$$ tail.$$
rm -f head.$$ body.$$ body2.$$ tail.$$



