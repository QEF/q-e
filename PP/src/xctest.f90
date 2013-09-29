!
! Copyright (C) 2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM xctest
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode
  USE kinds, ONLY: DP
  USE funct, ONLY: set_dft_from_indices
  IMPLICIT NONE
  INTEGER :: nnr = 1000
  INTEGER :: nspin = 2
  real(DP), ALLOCATABLE :: rhor( :, : )
  real(DP), ALLOCATABLE :: grhor( :, :, : )
  INTEGER :: iexch,icorr,igcx,igcc,inlc
  INTEGER :: nproc, mpime

  CALL mp_startup( )

  if ( ionode ) then
  iexch=1
  icorr=3
  igcx=1
  igcc=3
  inlc=0
  CALL set_dft_from_indices(iexch,icorr,igcx,igcc,inlc)

  OPEN(unit=17,form='unformatted',status='old')
  READ(17) nnr, nspin
  ALLOCATE(rhor( nnr, nspin ))
  ALLOCATE(grhor( nnr, 3, nspin ))
  READ(17) rhor
  READ(17) grhor
  CLOSE(17)

  !CALL test_gcxc( nnr, nspin, rhor, grhor )
  CALL test_xc( nnr, nspin, rhor, grhor )

  end if
  CALL mp_global_end()
END PROGRAM xctest

SUBROUTINE test_gcxc( nnr, nspin, rhor, grhor )
  USE kinds, ONLY: DP
!  use funct, only: gcxc
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nnr, nspin
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: grhor( nnr, 3, nspin )
  !
  real(DP), PARAMETER :: epsr = 1.0d-10, epsg = 1.0d-10
  real(DP), PARAMETER :: e2   = 1.0d0
  real(DP) :: grho2( nspin )
  real(DP) :: arho, segno
  real(DP) :: sx_w, sc_w, v1x_w, v2x_w, v1c_w, v2c_w
  real(DP) :: sx, sc, v1x, v2x, v1c, v2c
  real(DP) :: sx_m, sc_m, v1x_m, v2x_m, v1c_m, v2c_m
  real(DP) :: sx_d, sc_d, v1x_d, v2x_d, v1c_d, v2c_d
  INTEGER :: k, is, ipol

    DO k = 1, nnr
       !
       !
       DO is = 1, nspin
          grho2 (is) = grhor(k, 1, is)**2 + grhor(k, 2, is)**2 + grhor(k, 3, is)**2
       ENDDO
       !
       !
       IF (nspin == 1) THEN
          !
          !    This is the spin-unpolarised case
          !
          arho = abs (rhor (k, 1) )
          segno = sign (1.d0, rhor (k, 1) )
          IF (arho > epsr .and. grho2 (1) > epsg) THEN

             ! call gcxc (arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c)

             CALL becke88 (arho, grho2(1), sx, v1x, v2x)
             CALL wrap_b88 (arho, grho2(1), sx_w, v1x_w, v2x_w)  ! DEBUG
             CALL glyp (arho, grho2(1), sc, v1c, v2c)
             CALL wrap_glyp (arho, grho2(1), sc_w, v1c_w, v2c_w)  ! DEBUG

             sx_d = (sx_w - sx) / (abs(sx) + abs(sx_w))
             sc_d = (sc_w - sc) / (abs(sc) + abs(sc_w))
             v1x_d = (v1x_w - v1x) / (abs(v1x) + abs(v1x_w))
             v1c_d = (v1c_w - v1c) / (abs(v1c) + abs(v1c_w))
             v2x_d = (v2x_w - v2x) / (abs(v2x) + abs(v2x_w))
             v2c_d = (v2c_w - v2c) / (abs(v2c) + abs(v2c_w))


             WRITE(18,*) arho,grho2(1), sx_d, sc_d
             WRITE(19,*) arho,grho2(1), v1x_d, v1c_d
             WRITE(20,*) arho,grho2(1), v2x_w, v2x, v2x_d
             WRITE(21,*) arho,grho2(1), v2c_w, v2c, v2c_d

             !
             ! first term of the gradient correction : D(rho*Exc)/D(rho)

             ! v (k, 1) = v (k, 1) + e2 * (v1x + v1c)

             ! HERE h contains D(rho*Exc)/D(|grad rho|) / |grad rho|
             !
             ! h (k, 1, 1) = e2 * (v2x + v2c)
             ! etxc = etxc + e2 * (sx + sc) * segno

          ELSE
             ! h (k, 1, 1) = 0.d0
             sx = 0.0d0
             sc = 0.0d0
          ENDIF
          !
       ENDIF
       !
    ENDDO

    RETURN
END SUBROUTINE test_gcxc

!
!
!

SUBROUTINE test_xc( nnr, nspin, rhor, grhor )
  USE kinds, ONLY: DP
  USE funct, ONLY: get_iexch, get_icorr, get_igcx, get_igcc

  IMPLICIT NONE
  INTEGER, INTENT(in) :: nnr, nspin
  real(DP) :: rhor( nnr, nspin )
  real(DP) :: grhor( nnr, 3, nspin )
  !
  real(DP) :: rhon( nnr, nspin )
  real(DP) :: grhon( nnr, 3, nspin )
  real(DP) :: exc, excn, rhod, grhod
  INTEGER :: ir, is, ipol
  INTEGER iexch,icorr,igcx,igcc


  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  rhon  = rhor
  grhon = grhor
  !
  ! original CP xc selection
  !
      IF (iexch==1.and.icorr==1.and.igcx==0.and.igcc==0) THEN
         ! LDA (Perdew-Zunger)
         CALL expxc(nnr,nspin,rhor,exc)
      ELSEIF (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) THEN
         ! PW91
         CALL ggapwold(nnr,nspin,grhor,rhor,exc)
      ELSEIF (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) THEN
         ! BLYP
         CALL ggablyp4(nnr,nspin,grhor,rhor,exc)
      ELSEIF (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) THEN
         ! PBE
         CALL ggapbe(nnr,nspin,grhor,rhor,exc)
      ELSE
         CALL errore('exc-cor','no such exch-corr',1)
      ENDIF
  !
  ! Wrapper to PW xc selection
  !
  CALL exch_corr_cp(nnr,nspin,grhon,rhon,excn)
  !
  WRITE(6,*) 'EXC = ', exc, excn
  DO is = 1, nspin
    DO ir = 1, nnr
      rhod = abs( rhor( ir, is ) - rhon( ir, is ) ) / ( abs( rhor( ir, is ) ) + abs( rhon( ir, is ) ) )
      WRITE(18,100) ir,is,rhod
    ENDDO
  ENDDO
  DO is = 1, nspin
    DO ir = 1, nnr
      DO ipol = 1, 3
      grhod = abs( grhor( ir, ipol, is ) - grhon( ir, ipol, is ) ) / &
             ( abs( grhor( ir, ipol, is ) ) + abs( grhon( ir, ipol, is ) ) )
      WRITE(19,100) ir,is,grhod
      ENDDO
    ENDDO
  ENDDO
100 FORMAT( I5, I2, 1X, E15.8, 1X, E15.8 )
END SUBROUTINE test_xc
