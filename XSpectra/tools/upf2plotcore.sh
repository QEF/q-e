#!/bin/sh

LANG=C

cat $1 | awk '


#<UPF v1>

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

#</UPF v1>


#<UPF v2>

/mesh_size=/ {
  nris = $0
  gsub ( "[ ]*mesh_size=\"", "", nris )
  gsub ( "\"[ ]*$", "", nris )
}

/<PP_R /,/<\/PP_R>/ {
  if ( $1 != "<PP_R" && $1 != "</PP_R>" ) {
    for ( i = 1; i <= NF; i++ ) {
      r[++val] = 1*$i
    }
  }
}

/<PP_GIPAW_CORE_ORBITAL.[0-9] / {
  proj++
  val_wfc = 0
}
/<PP_GIPAW_CORE_ORBITAL.[0-9]/,/<\/PP_GIPAW_CORE_ORBITAL.[0-9]>/ {
  st = $1; gsub ( "[0-9]*", "", st )
#  if ( index( $0, "label=" ) != 0 && index( $0, "label=\"1S\"" ) == 0 ) { next }

  if ( index( $0, " n=\"" ) != 0 || index( $0, "^n=\"" ) != 0 ) {
    n[proj] = $0
    gsub( "[a-z0-9A-Z=_.<> \"]*n=\"", "", n[proj] )
    gsub( "[\">]*$", "", n[proj] )
    n[proj] = int ( n[proj] + 1e-5 )
  }

  if ( index( $0, " l=\"" ) != 0 || index( $0, "l=\"" ) == 1 ) {
    l[proj] = $0
    gsub( "[a-z0-9A-Z=_.<> \"]*l=\"", "", l[proj] )
    gsub( "[\">]*$", "", l[proj] )
    l[proj] = int ( l[proj] + 1e-5 )
  }

  if ( index( $0, ">" ) != 0 || index( $0, "<" ) != 0 ) { next }
  for ( i = 1; i <= NF; i++ ) {
    f[++val_wfc,proj] = $i
  }

}

#</UPF v2>


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
