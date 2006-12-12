#!/bin/sh
# namedep.sh -- script that computes dependencies on Fortran 90 modules

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# first argument is mandatory
if test $# = 0
then
    echo usage: $0 name [files]
    exit 1
fi

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname

# module, function or subroutine whose dependencies must be computed
name=$1
shift

# list of files to be searched
sources_all=`ls */*.f90`
if test $# = 0 ; then sources="$sources_all"
else sources="$* /dev/null" ; fi

# search for declaration of $name
# caution: must not select names that _contain_ $name
decls=`egrep -ni -e "^ *subroutine  *$name *(\(.*)?$" \
                 -e "^ *function  *$name *(\(.*)?$" \
                 -e "^ *module  *$name *$" \
       $sources | sed 's/[:(]/ /g' | awk '{print $1 "@" $2 "@" $4}'`

num=`echo $decls | wc | awk '{print $2}'`
if test $num = 0
then
    echo error: $name not found
    exit 1
elif test $num -gt 1
then
    # $name is defined in more than one place, must choose one
    echo error: there are multiple declarations:
    for decl in $decls
    do
      file=`echo $decl | sed 's/@/ /g' | awk '{print $1}'`
      echo "    $name [$file]"
    done
    echo please specify file
    exit 1
fi

# build list of all module declarations
# list format is: file_name starting_line module_name
egrep -ni "^ *module  *[a-zA-Z_][a-zA-Z_]*" $sources_all |
grep -iv procedure | # exclude "module procedure" declarations
sed 's/:/ /g' | awk '{print $1, $2, $4}' > namedep.sh.tmp1

decl=`echo $decls | sed 's/@/ /g'`
file=`echo $decl | awk '{print $1}'`
echo $name [$file]:

# find starting and ending line
start=`echo $decl | awk '{print $2}'`
end=`egrep -ni -e "^ *end  *subroutine  *$name *$" \
               -e "^ *end  *function  *$name *$" \
               -e "^ *end  *module  *$name *$" \
     $file | sed 's/:.*//'`

# look for use declarations
modules=`sed -n "$start,${end}p" $file | egrep -i "^ *use " |
         sed 's/,.*//' |             # remove ", only: ..."
         awk '{print tolower($2)}' | # cast module name to lowercase
         sort | uniq                 # remove duplicates`

# look for recursive dependencies
modules_prev=""
until test "$modules_prev" = "$modules"
do
  modules_tested="$modules_prev"
  modules_prev="$modules"
  for module in $modules
  do
    # skip module if already tested
    if test "`echo $modules_tested | tr ' ' '\n' | grep ^$module\$`" = ""
    then
	mdecl=`egrep -i " $module *$" namedep.sh.tmp1`
	file=`echo $mdecl | awk '{print $1}'`

	# find starting and ending line
	start=`echo $mdecl | awk '{print $2}'`
	end=`egrep -ni "^ *end  *module  *$module *$" $file | sed 's/:.*//'`

	# look for use declarations
	recur=`sed -n "$start,${end}p" $file | egrep -i "^ *use " |
               sed 's/,.*//' |             # remove ", only: ..."
               awk '{print tolower($2)}' | # cast module name to lowercase
               sort | uniq                 # remove duplicates`
	modules="$modules $recur"
    fi
  done
  # remove duplicates
  modules=`echo $modules | tr " " "\n" | sort | uniq`
done

# print final list of modules
for module in $modules
do
  mdecl=`grep -i $module namedep.sh.tmp1`
  file=`echo $mdecl | awk '{print $1}'`
  echo "    $module [$file]"
done

# remove temporary file
rm -f namedep.sh.tmp1
