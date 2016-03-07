  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine reset_sym (xq, nsym, s, invs, irt, rtau)
  !-----------------------------------------------------------------------
  !
  ! output values of symmetry arrays (nsym, s, rtau, irt) are those
  ! appropriate to the small-qroup of q.
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,         only : DP
  USE phcom,         only : t, tmq, npertx, u, npert, nirr
  USE lr_symm_base,  ONLY : irotmq, minus_q, nsymq, gi, gimq, irgq
  USE pwcom,         only : at, bg
  USE symm_base,     ONLY : sname, ftau, invsym, time_reversal, t_rev, &
                            find_sym, set_sym_bl
  USE ions_base,     only : nat, tau, ityp
  USE control_flags, only : iverbosity, noinv, modenum
  USE cryst_ph,      ONLY : magnetic_sym
  !
  implicit none
  real(kind=DP) :: xq (3)
  ! input: q vector

  !-output variables
  integer :: nsym, s (3, 3, 48), invs (48), irt (48, nat)
  ! output: number of symmetry operations
  ! output: the first nq matrices are those that generate the star of q
  !         starting from it
  ! output: list of inverse operation indices
  ! output: for each atom gives the rotated atom

  real(kind=DP) :: rtau (3, 48, nat)
  ! output: for each atom and rotation gives the R vector involved
  ! output: list of vectors in the star of q
  !
  ! Local variables
  !
  integer :: nrot, isym, jsym, table (48, 48), &
       i, j, nks0, npk0
  ! number of symmetry ops. of bravais lattice.
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! group table
  ! counter on q-vectors
  ! generic counter
  ! number of dummy k-points
  ! maximum allowed number of dummy k-points
  ! dummy (zero) value of iswitch passed to sgama

  integer ::  irot, jrot, ipol, jpol, temp, na
  ! counter over the rotations
  ! counter over the rotations
  ! counter over the polarizations
  ! counter over the polarizations

  real(kind=DP) :: xk0 (3), wk0(1), zero (3), mdum(3,nat)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! dummy k-points list
  ! a zero vector: used in eqvect and as dummy q-vector in sgama

  logical :: sym (48)
  ! .t. if the crystal has inversion
  ! dummy output from sgama
  ! input for sgama
  !
  !  initialize dummy k-point list and zero vector
  !
  npk0 = 1
  nks0 = 1
  wk0(:) = 1.d0
  xk0(:)= 0.d0
  zero(:) = 0.d0
  !
  !  generate transformation matrices for the bravais lattice
  !
  CALL set_sym_bl ( )

!  if (noinv) then
  IF (noinv) then
     jsym = 0
     DO isym = 1, nrot
        IF ( s (1, 3, isym) == 0 .and. s (3, 1, isym) == 0 .and. &
             s (2, 3, isym) == 0 .and. s (3, 2, isym) == 0 .and. &
             s (3, 3, isym) == 1) then
           jsym = jsym + 1
           DO i = 1, 3
              DO j = 1, 3
                 s (i, j, jsym) = s (i, j, isym)
              ENDDO
           ENDDO
           sname (jsym) = sname (isym)
        ENDIF
     ENDDO
     nrot = jsym
  ENDIF
  !
  ! extract from it the crystal symmetry group by calling sgama
  !
!  CALL sgama (nrot, nat, s, sname, at, bg, tau, ityp, nsym, dfftp%nr1, &
!       dfftp%nr2, dfftp%nr3, irt, ftau, npk0, nks0, xk0, wk0, invsym, minus_q, xq, &
!       iswitch, modenum, .false., mdum)

  CALL find_sym ( nat, tau, ityp, 6, 6, 6, .not.time_reversal, mdum )
  call smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)

  IF ( .not. time_reversal ) THEN
     minus_q=.false.
  ENDIF
  !
  if (modenum .ne. 0) then
     call sgam_ph_new (at, bg, nrot, s, irt, tau, rtau, nat)
     call mode_group (modenum, xq, at, bg, nat, nrot, s, irt, rtau, &
          sym, minus_q)
  endif
  !
  !nqx = nk1*nk2*nk3
  CALL irreducible_BZ (nrot, s, nsym, time_reversal, magnetic_sym, &
                       at, bg, npk0, nks0, xk0, wk0, t_rev)

  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  !
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        do ipol = 1, 3
           do jpol = 1, 3
              temp = s (ipol, jpol, jrot)
              s (ipol, jpol, jrot) = s (ipol, jpol, irot)
              s (ipol, jpol, irot) = temp
           enddo
           ftau (ipol, jrot) = ftau (ipol, irot)
        enddo
        do na = 1, nat
           irt (jrot, na) = irt (irot, na)
        enddo
        sname (jrot) = sname (irot)
        t_rev(jrot) = t_rev(irot)
     endif
  enddo
  if (jrot.ne.nsym) call errore ('sgama', 'unexpected', 1)
  !
  ! Sets to zero the first matrix that is not a symmetry of the crystal.
  ! This will be used by d3toten program.
  !
  if (nrot.lt.48) then
     do ipol = 1, 3
        do jpol = 1, 3
           s (ipol, jpol, nrot + 1) = 0
        enddo
     enddo
  endif
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = .true.
  irot = nsym / 2 + 1
  do ipol = 1, 3
     do jpol = 1, 3
        invsym = invsym.and.s (ipol, jpol, irot) .eq. - s (ipol, jpol, 1)
     enddo
  enddo  

  DO isym = 1, nsym
     sym (isym) = .true.
  ENDDO
  CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
  !
  ! computes the inverse of each matrix
  !
  CALL multable (nsym, s, table)
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (table (isym, jsym) .eq.1) invs (isym) = jsym
     ENDDO
  ENDDO
  !
  IF (nsym.gt.1) then
     CALL set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, npert, &
          nirr, gi, gimq, iverbosity)
  ELSE
     CALL set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, npert, &
          nirr, gi, gimq, iverbosity)
  ENDIF
  !
  end subroutine reset_sym
