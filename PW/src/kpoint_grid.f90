!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                         bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nrot, npk, k1, k2, k3, nk1, nk2, nk3, &
                        t_rev(48), s(3,3,48)
  LOGICAL, INTENT(in):: time_reversal, skip_equivalence
  real(DP), INTENT(in):: bg(3,3)
  !
  INTEGER, INTENT(out) :: nks
  real(DP), INTENT(out):: xk(3,npk)
  real(DP), INTENT(out):: wk(npk)
  ! LOCAL:
  real(DP), PARAMETER :: eps=1.0d-5
  real(DP) :: xkr(3), fact, xx, yy, zz
  real(DP), ALLOCATABLE:: xkg(:,:), wkk(:)
  INTEGER :: nkr, i,j,k, ns, n, nk
  INTEGER, ALLOCATABLE :: equiv(:)
  LOGICAL :: in_the_list
  !
  nkr=nk1*nk2*nk3
  ALLOCATE (xkg( 3,nkr),wkk(nkr))
  ALLOCATE (equiv( nkr))
  !
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  DO nk=1,nkr
     equiv(nk)=nk
  ENDDO

  IF ( skip_equivalence ) THEN
    CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
    wkk = 1.d0
  ELSE
    DO nk=1,nkr
    !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk)   = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns=1,nrot
           DO i=1,3
              xkr(i) = s(i,1,ns) * xkg(1,nk) &
                     + s(i,2,ns) * xkg(2,nk) &
                     + s(i,3,ns) * xkg(3,nk)
              xkr(i) = xkr(i) - nint( xkr(i) )
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_the_list = abs(xx-nint(xx))<=eps .and. &
                         abs(yy-nint(yy))<=eps .and. &
                         abs(zz-nint(zz))<=eps
           IF (in_the_list) THEN
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              IF (n>nk .and. equiv(n)==n) THEN
                 equiv(n) = nk
                 wkk(nk)=wkk(nk)+1.0d0
              ELSE
                 IF (equiv(n)/=nk .or. n<nk ) CALL errore('kpoint_grid', &
                    'something wrong in the checking algorithm',1)
              ENDIF
           ENDIF
           IF ( time_reversal ) THEN
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_the_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                                 .and. abs(zz-nint(zz))<=eps
              IF (in_the_list) THEN
                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 IF (n>nk .and. equiv(n)==n) THEN
                    equiv(n) = nk
                    wkk(nk)=wkk(nk)+1.0d0
                 ELSE
                    IF (equiv(n)/=nk.or.n<nk) CALL errore('kpoint_grid', &
                    'something wrong in the checking algorithm',2)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !  count irreducible points and order them

  nks=0
  fact=0.0d0
  DO nk=1,nkr
     IF (equiv(nk)==nk) THEN
        nks=nks+1
        IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
        wk(nks) = wkk(nk)
        fact    = fact+wk(nks)
        !  bring back into to the first BZ
        DO i=1,3
           xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
        ENDDO
     ENDIF
  ENDDO
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk,bg,1)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = wk(nk)/fact
  ENDDO

  DEALLOCATE(equiv)
  DEALLOCATE(xkg,wkk)

  RETURN
END SUBROUTINE kpoint_grid
!-----------------------------------------------------------------------
SUBROUTINE kpoint_grid_efield (at, bg, npk, &
                         k1,k2,k3, nk1,nk2,nk3, nks, xk, wk, nspin)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
! for Berry's phase electric field
!
  USE kinds, ONLY : DP
  USE bp,    ONLY : nppstr_3d, nx_el, l3dstring, efield_cart, efield_cry,&
                    transform_el
  USE io_global,  ONLY : stdout
  USE noncollin_module,   ONLY : noncolin
  USE matrix_inversion

  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: npk, k1, k2, k3, nk1, nk2, nk3,nspin
  real(DP), INTENT(in):: bg(3,3), at(3,3)
  !
  INTEGER, INTENT(out) :: nks
  real(DP), INTENT(out):: xk(3,npk)
  real(DP), INTENT(out):: wk(npk)

  INTEGER :: i,j,k,n,nk,m
  INTEGER :: nppstr_max
  real(DP) :: fact, sca
  real(DP) :: cry_to_cart(3,3)
  real(DP) :: bg_n(3,3)
  !
  !
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xk(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xk(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xk(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  nks=nk1*nk2*nk3
  !  go to cartesian axis (in units 2pi/a0)
  CALL cryst_to_cart(nks,xk,bg,1)
  fact=1.d0/dble(nks)
  !  normalize weights to one
  DO nk=1,nks
     wk(nk) = fact
  ENDDO

!setup nppstr_3d
  nppstr_3d(1)=nk1
  nppstr_3d(2)=nk2
  nppstr_3d(3)=nk3

!allocate and set up correspondence
  nppstr_max=nk1*nk2*nk3

  IF(noncolin) THEN
     ALLOCATE(nx_el(nppstr_max,3))
  ELSE
     ALLOCATE(nx_el(nppstr_max*nspin,3))
  END IF

 !establih correspondence

   DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           nx_el(n,3)=n
           m = (i-1) + (k-1)*nk1 + (j-1)*nk3*nk1 + 1
           nx_el(m,1)=n
           m = (j-1) + (i-1)*nk2 + (k-1)*nk1*nk2 + 1
           nx_el(m,2)=n
      ENDDO
     ENDDO
  ENDDO

  IF(nspin==2) THEN
     DO i=1,nks
        nx_el(i+nks,:)=nx_el(i,:)+nks
     ENDDO
  ENDIF
  l3dstring=.true.

  DO i=1,3
     sca=at(1,i)**2.d0+at(2,i)**2.d0+at(3,i)**2.d0
     sca=dsqrt(sca)
     bg_n(1:3,i)=(1.d0/sca)*at(1:3,i)
  ENDDO

  DO i=1,3
     DO j=1,3
        cry_to_cart(j,i) = bg_n(1,j)*bg_n(1,i) + &
                           bg_n(2,j)*bg_n(2,i) + &
                           bg_n(3,j)*bg_n(3,i)
     ENDDO
 ENDDO
 CALL  invmat (3, cry_to_cart, transform_el)

! calculate EFFECTIVE electric field on crystal axis

  efield_cry(:)=0.d0
  DO i=1,3
        efield_cry(i) = efield_cry(i) + efield_cart(1)*bg_n(1,i) + &
                                        efield_cart(2)*bg_n(2,i) + &
                                        efield_cart(3)*bg_n(3,i)
  ENDDO

  RETURN

 END SUBROUTINE kpoint_grid_efield
