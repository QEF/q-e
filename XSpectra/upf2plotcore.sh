#!/bin/sh

LANG=C

cat $1 | awk '

/Number of points in mesh/ {
  nris = $1
}
/PP/ {
  val = 0
}

/<PP_R>/,/<\/PP_R>/ {
  if ( $1 != "<PP_R>" && $1 != "</PP_R>" ) {
    for ( i = 1; i <= NF; i++ ) {
      r[++val] = 1*$i
    }
  }
}

/<PP_GIPAW_CORE_ORBITAL>/ {
  proj++
}

/<PP_GIPAW_CORE_ORBITAL>/,/<\/PP_GIPAW_CORE_ORBITAL>/ {
  if ( $1 == "<PP_GIPAW_CORE_ORBITAL>" ) {
    getline
    n[proj] = $1
    l[proj] = $2
    label[proj] = $5
    next
  }
  if ( $1 != "<PP_GIPAW_CORE_ORBITAL>" && $1 != "</PP_GIPAW_CORE_ORBITAL>" ) {
    for ( i = 1; i <= NF; i++ ) {
      f[++val,proj] = $i
    }
  }
}


END {

  nprojs = proj
  printf "# number of core states %d = ", nprojs, " n,l = "
  for ( proj = 1; proj <= nprojs; proj++ ) {
    printf " %d %d; ", n[proj], l[proj]
  }
  printf "\n"
  for ( proj = 1; proj <= nprojs; proj++ ) {
    for ( ri = 1; ri <= nris; ri++ ) {
      print r[ri], f[ri,proj]
    }
    print ""
  }
}

'
