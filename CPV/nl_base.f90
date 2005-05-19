!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!  ----------------------------------------------
!  BEGIN manual

   MODULE nl_base

!  this module handles nonlocal pseudopotential calculations
!  ----------------------------------------------
!  routines in this module:
!  REAL(dbl) FUNCTION ene_nl(fnl,wsg,occ,ngh,nspnl,na)
!  REAL(dbl) FUNCTION ene_nl_kp(fnlk,wk,wsg,occ,nk,ngh,nspnl,na)
!  ----------------------------------------------
!  END manual

!  end of module-scope declarations
!  ----------------------------------------------

! ...   declare modules
        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: ene_nl

      CONTAINS

!  subroutines
!
!  ----------------------------------------------
!  BEGIN manual

        REAL(dbl) FUNCTION ene_nl(fnl, wsg, fi, nspnl, na)

!  this function computes and returns the nonlocal contribution to total
!  energy, for both Gamma point and Kpoints calculations
!
!    ene_nl =  (sum over ia,ib,igh,is) f(ib) wsg(igh,is) fnl(ia,igh,ib)**2
!
!    fi(ib,1)       = occupation numbers
!    fnl(ia,igh,ib) = Kleinman-Bylander factor (see nlsm1)
!    wsg(igh,is)    = inverse denominator in KB's formula
!
!    ia = index of ion
!    ib = index of band
!    igh = index of orbital
!    is = index of atomic species
!  ----------------------------------------------
!  END manual

          USE pseudo_projector, ONLY: projector

! ...     declare function arguments
          TYPE (projector), INTENT(IN) :: fnl
          REAL(dbl), INTENT(IN) :: wsg(:,:), fi(:)
          INTEGER, INTENT(IN) :: nspnl
          INTEGER, INTENT(IN) :: na(:)

! ...     declare other variables
          INTEGER   :: igh, isa, is, ia, ib, nb, ngh
          REAL(dbl) :: enl, fsum
          COMPLEX(dbl) :: tt

! ...     end of declarations
!  ----------------------------------------------

          IF( fnl%gamma_symmetry ) THEN
            nb  = MIN( SIZE( fnl%r, 3), SIZE(fi) )
            ngh = SIZE( fnl%r, 2)
          ELSE
            nb  = MIN( SIZE( fnl%c, 3), SIZE(fi) )
            ngh = SIZE( fnl%c, 2)
          END IF

          enl=0.d0
          DO igh = 1, ngh
            DO ib = 1, nb
              isa = 0
              DO is = 1, nspnl
                fsum = 0.0d0
                IF( fnl%gamma_symmetry ) THEN
                  DO ia = 1, na(is)
                    fsum = fsum + fnl%r(isa+ia,igh,ib)**2
                  END DO
                ELSE
                  DO ia = 1, na(is)
                    tt = fnl%c(isa+ia,igh,ib)
                    fsum = fsum + REAL( CONJG(tt) * tt )
                  END DO
                END IF
                enl = enl + fi(ib) * wsg(igh, is) * fsum
                isa = isa + na(is)
              END DO
            END DO
          END DO

          ene_nl = enl

          RETURN
        END FUNCTION ene_nl


   END MODULE nl_base
