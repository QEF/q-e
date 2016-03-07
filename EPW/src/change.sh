# Copyright (C) 2015 Samuel Ponce, Feliciano Giustino
#
# Small script to mass replace argument in lots of files
#

for file in `find . -maxdepth 1 -type f`
do
  sed -i "s,io_global\,,io_epw\,   ,"  $file
done
