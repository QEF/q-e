# Copyright (C) 2001-2022 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_GIT], [

  AC_CHECK_PROG(git, git, git)
  if test "$git" = ""; then
    AC_MSG_ERROR([git needed])
  fi

  if test -d $topdir/.git ; then
    echo Source files are cloned from a git repository.
    gitbranch=`(cd $topdir && $git rev-parse --abbrev-ref HEAD)`
    echo On git branch $gitbranch
    githash=`(cd $topdir && $git describe --always --dirty --abbrev=40 --match="NoTagWithThisName")` 
    echo On git commit hash $githash
    (cd $topdir && $git submodule init)
  fi

  AC_SUBST(git)

  ]
)
