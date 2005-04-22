!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 14 08:42:17 MET 1999
!  ----------------------------------------------
!  BEGIN manual

      MODULE nl_base

!  this module handles nonlocal pseudopotential calculations
!  ----------------------------------------------
!  routines in this module:
!  REAL(dbl) FUNCTION ene_nl(fnl,wsg,occ,ngh,nspnl,na)
!  REAL(dbl) FUNCTION ene_nl_kp(fnlk,wk,wsg,occ,nk,ngh,nspnl,na)
!  SUBROUTINE dedh_nls(wnl,wnla,eigr,ll,auxc,gwork,na)
!  SUBROUTINE dedh_nlp(gagx,wnl,wnla,eigr,ll,kk,auxc,gwork,na)
!  SUBROUTINE dedh_nld(gagx,wnl,wnla,eigr,ll,kk,dm,auxc,gwork,na)
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
!  ----------------------------------------------

!  ----------------------------------------------
        SUBROUTINE dedh_nls( wnl, wnla, eigr, ll, auxc, gwork, na )

!  this routine computes the nonlocal contribution of s-orbitals to
!  cell derivatives of total energy
!  version for Gamma-point calculations
!
!  orbitals: s(r) = 1
!  ----------------------------------------------

        USE gvecw,              ONLY: ngw
        USE reciprocal_vectors, ONLY: gstart

! ...     declare subroutine arguments
          REAL(dbl), INTENT(IN) :: wnl(:,:)
          REAL(dbl), INTENT(IN) :: wnla(:,:)
          COMPLEX(dbl), INTENT(IN) :: eigr(:,:)
          INTEGER, INTENT(IN) :: ll,na
          COMPLEX(dbl), INTENT(OUT) :: auxc(:,:)
          REAL(dbl), INTENT(IN)  :: gwork(:)

! ...     declare other variables
          INTEGER   :: ia, ig
          REAL(dbl) :: arg
          REAL(dbl), ALLOCATABLE :: gwtmp(:)

! ...     end of declarations
!  ----------------------------------------------

          IF(ll.GT.0) THEN
            ALLOCATE(gwtmp(ngw))
            gwtmp(1) = 0.0d0
            DO ig = gstart, ngw
              gwtmp(ig) = gwork(ig) * (wnl(ig,ll)-wnla(ig,ll))
            END DO
            DO ia = 1, na
              auxc(1,ia) = CMPLX( 0.0d0 )
              DO ig = gstart, ngw
                auxc(ig,ia) = gwtmp(ig) * eigr(ig,ia)
              END DO
            END DO
            DEALLOCATE(gwtmp)
          ELSE
            auxc = CMPLX( 0.0d0 )
          END IF

          RETURN
        END SUBROUTINE dedh_nls

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE dedh_nlp( gagx, wnl, wnla, eigr, ll, kk, auxc, gwork, na )

!  this routine computes the nonlocal contribution of p-orbitals to
!  cell derivatives of total energy
!  version for Gamma-point calculations
!
!  orbitals: p_x(r) = x/r
!            p_y(r) = y/r
!            p_z(r) = z/r
!  ----------------------------------------------

        USE reciprocal_vectors, ONLY: gstart
        USE gvecw,              ONLY: ngw

! ...     declare subroutine arguments
          REAL(dbl), INTENT(IN)  :: gagx(:,:)
          REAL(dbl), INTENT(IN) :: wnl(:,:)
          REAL(dbl), INTENT(IN) :: wnla(:,:)
          COMPLEX(dbl), INTENT(IN) :: eigr(:,:)
          INTEGER, INTENT(IN) :: ll,na,kk
          COMPLEX(dbl), INTENT(OUT) :: auxc(:,:)
          REAL(dbl), INTENT(IN)  :: gwork(:)

! ...     declare other variables
          INTEGER  :: ia, ig
          REAL(dbl), ALLOCATABLE :: gwtmp(:)
          COMPLEX(dbl), PARAMETER :: uimag = (0.0d0,1.0d0)

! ...     end of declarations
!  ----------------------------------------------

          IF(ll.GT.0) THEN
            ALLOCATE(gwtmp(ngw))
            gwtmp(1) = 0.0d0
            DO ig = gstart, ngw
              gwtmp(ig) = gagx(kk,ig) * gwork(ig)* &
              ( 3.d0 * wnl(ig,ll) - wnla(ig,ll) )
            END DO
            DO ia = 1, na
              auxc(1,ia) = CMPLX( 0.0d0 )
              DO ig = gstart, ngw
                auxc(ig,ia) = uimag * gwtmp(ig) * eigr(ig,ia)
              END DO
            END DO
            DEALLOCATE(gwtmp)
          ELSE
            auxc = CMPLX( 0.0d0 )
          END IF

          RETURN
        END SUBROUTINE dedh_nlp

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE dedh_nld( gagx, wnl, wnla, eigr, ll, kk, dm, auxc, gwork, na )

!  this routine computes the nonlocal contribution of d-orbitals to
!  cell derivatives of total energy
!  version for Gamma-point calculations
!
!  orbitals: d_m(r) = gkl(x,y,z,r,m)
!  ----------------------------------------------

        USE reciprocal_vectors, ONLY: gstart
        USE gvecw,              ONLY: ngw

! ...     declare subroutine arguments
          REAL(dbl), INTENT(IN)  :: gagx(:,:)
          REAL(dbl), INTENT(IN) :: wnl(:,:)
          REAL(dbl), INTENT(IN) :: wnla(:,:)
          COMPLEX(dbl), INTENT(IN) :: eigr(:,:)
          INTEGER, INTENT(IN) :: ll,na,kk
          COMPLEX(dbl), INTENT(OUT) :: auxc(:,:)
          REAL(dbl), INTENT(IN)  :: dm

! ...     declare other variables
          REAL(dbl), INTENT(IN)  :: gwork(:)
          INTEGER :: ia, ig
          REAL(dbl), ALLOCATABLE :: gwtmp(:)

! ...     end of declarations
!  ----------------------------------------------

          IF(ll.GT.0) THEN
            ALLOCATE(gwtmp(ngw))
            gwtmp(1) = 0.0d0
            DO ig = gstart, ngw
               gwtmp(ig) = gwork(ig) * gagx(kk,ig) *  &
               ( 5.d0 * wnl(ig,ll) - wnla(ig,ll) )
               gwtmp(ig) = gwtmp(ig) - 2.0d0/3.0d0 * dm * wnl(ig,ll)
            END DO
            DO ia= 1 , na
              auxc(1,ia) = CMPLX( 0.0d0 )
              DO ig = gstart, ngw
                auxc(ig,ia) = - gwtmp(ig) * eigr(ig,ia)
              END DO
            END DO
            DEALLOCATE(gwtmp)
          ELSE
            auxc = CMPLX( 0.0d0 )
          END IF

          RETURN
        END SUBROUTINE dedh_nld

!  ----------------------------------------------
!  ----------------------------------------------

      END MODULE nl_base

