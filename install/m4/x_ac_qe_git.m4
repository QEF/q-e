# Copyright (C) 2001-2022 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_GIT], [

  AC_CHECK_PROG(git, git, git)
  if test "$git" = ""; then
    AC_MSG_ERROR([git needed])
  fi

  if test -d $topdir/.git ; then
    echo Source files are cloned from a git repository.
    echo On git branch `$git rev-parse --abbrev-ref HEAD`
    echo On git commit hash `$git describe --always --dirty --abbrev=40 --match="NoTagWithThisName"`
    git submodule init
  fi

  AC_SUBST(git)

  ]
)
