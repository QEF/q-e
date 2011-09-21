!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE writed3dyn_5 (d3dyn_x, filename, isw)
  !-----------------------------------------------------------------------
  !
  !     writes in a file the third derivative of dynamical matrix
  !     isw = +1  :  d3dyn_x is in cartesian axis
  !     isw = -1  :  rotates d3dyn_x from the basis of pattern to
  !                      cartesian axis
  !
  USE ions_base,  ONLY : nat
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode
  USE pwcom
  USE phcom
  USE d3com
  !
  IMPLICIT NONE
  !
  INTEGER :: isw, iud3dyn, n_d3, na, nb, icart, jcart, kcart, na_i, &
       na_j, na_k
  ! input: switch
  ! index on cartesian coordinates
  ! index on cartesian coordinates
  ! index on cartesian coordinates
  ! index on modes
  ! index on modes
  ! index on modes

  COMPLEX (DP) :: d3dyn_x (3 * nat, 3 * nat, 3 * nat), work
  ! input: the third derivative of the dynamical matrix
  COMPLEX (DP), ALLOCATABLE :: aux (:,:,:)
  ! auxiliary space

  CHARACTER (len=*) :: filename
  ! input: the name of the file

  IF ( .NOT. ionode ) RETURN

  ALLOCATE  (aux( 3 * nat, 3 * nat, 3 * nat))
  IF (isw.EQ. + 1) THEN
     CALL zcopy (27 * nat * nat * nat, d3dyn_x, 1, aux, 1)
  ELSEIF (isw.EQ. - 1) THEN
     !
     !   Rotates third derivative of the dynamical basis from the basis
     !   of modes to cartesisn axis
     !
     DO kcart = 1, 3 * nat
        DO icart = 1, 3 * nat
           DO jcart = 1, 3 * nat
              work = (0.d0, 0.d0)
              DO na_k = 1, 3 * nat
                 DO na_i = 1, 3 * nat
                    DO na_j = 1, 3 * nat
                       work = work + CONJG (ug0 (kcart, na_k) ) * u (icart, na_i) &
                            * d3dyn_x (na_k, na_i, na_j) * CONJG (u (jcart, na_j) )
                    ENDDO
                 ENDDO
              ENDDO
              aux (kcart, icart, jcart) = work
           ENDDO
        ENDDO

     ENDDO

  ENDIF
  iud3dyn = 57

  OPEN (unit = iud3dyn, file = TRIM(filename), status = 'unknown')
  DO n_d3 = 1, 3 * nat
     WRITE (iud3dyn, * )
     WRITE (iud3dyn,  * ) '               modo:', n_d3
     WRITE (iud3dyn, * )
     DO na = 1, nat
        DO nb = 1, nat
           WRITE (iud3dyn, '(2i3)') na, nb
           DO icart = 1, 3
              WRITE (iud3dyn, '(3E24.12)') (aux (n_d3, icart + 3 * (na - 1) , &
                   jcart + 3 * (nb - 1) ) , jcart = 1, 3)
           ENDDO
        ENDDO
     ENDDO

  ENDDO

  CLOSE (iud3dyn)

  DEALLOCATE (aux)
  RETURN
END SUBROUTINE writed3dyn_5
