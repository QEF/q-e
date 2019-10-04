# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_AR], [

  # default from the environment (shouldn't be needed)
  ar=$AR    
  arflags=$ARFLAGS
  
  try_ar="ar"
  AC_MSG_CHECKING([setting AR... ])
  if test "$arch" = "necsx"; then
    try_ar="sxar"
  fi  
  if test "$ar" = "" ; then ar="$try_ar" ; fi
  AC_MSG_RESULT(${ar})
  AC_SUBST(ar)
  
  try_arflags="ruv"
  AC_MSG_CHECKING([setting ARFLAGS... ])
  if test "$arch" = "necsx"; then
        try_arflags="rv"
  fi
  if test "$arflags" = "" ; then arflags="$try_arflags" ; fi
  AC_MSG_RESULT(${arflags})  
  AC_SUBST(arflags)
  
  ]
)
