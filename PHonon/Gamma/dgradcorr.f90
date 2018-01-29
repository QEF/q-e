!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE dgradcor1 (dfft, rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, &
     drho, drhoc, nspin, g, dvxc)
  !     ===================
  !--------------------------------------------------------------------
  !  ADD Gradient Correction contibution to screening potential
  !  phonon calculation, half G-vectors
  USE kinds,      ONLY : DP
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  INTEGER, INTENT(IN) :: nspin

  REAL(DP), INTENT(IN) :: rho (dfft%nnr, nspin), grho (3, dfft%nnr, nspin), &
       g (3, dfft%ngm)
  REAL(DP), INTENT(OUT):: drho (dfft%nnr,nspin),&
       dvxc_rr(dfft%nnr, nspin, nspin), dvxc_sr (dfft%nnr, nspin, nspin), &
       dvxc_ss (dfft%nnr,nspin, nspin), dvxc_s (dfft%nnr, nspin, nspin)
  REAL(DP), INTENT(INOUT) ::  dvxc (dfft%nnr, nspin)
  COMPLEX(DP) :: drhoc(dfft%nnr, nspin)
  !
  INTEGER :: k, ipol, is, js, ks, ls
  real(DP) :: epsr, epsg, grho2
  COMPLEX(DP) :: s1
  COMPLEX(DP) :: a (2, 2, 2), b (2, 2, 2, 2), c (2, 2, 2), &
                      ps (2, 2), ps1 (3, 2, 2), ps2 (3, 2, 2, 2)
  REAL(DP), ALLOCATABLE  :: gdrho (:,:,:)
  REAL(DP), ALLOCATABLE :: h (:,:,:), dh (:)
  PARAMETER (epsr = 1.0d-6, epsg = 1.0d-10)

  ALLOCATE (gdrho( 3, dfft%nnr , nspin))
  ALLOCATE (h(  3, dfft%nnr , nspin))
  ALLOCATE (dh( dfft%nnr))

  h (:,:,:) = 0.d0
  DO is = 1, nspin
     CALL fft_gradient_g2r (dfft, drhoc(1, is), g, gdrho (1,1,is) )
  ENDDO
  DO k = 1, dfft%nnr
     grho2 = grho(1, k, 1)**2 + grho(2, k, 1)**2 + grho(3, k, 1)**2
     IF (nspin==1) THEN
        !
        !    LDA case
        !
        IF (abs (rho (k, 1) ) >epsr.and.grho2>epsg) THEN
           s1 = grho (1, k, 1) * gdrho (1, k, 1) + &
                grho (2, k, 1) * gdrho (2, k, 1) + &
                grho (3, k, 1) * gdrho (3, k, 1)
           !
           ! linear variation of the first term
           !
           dvxc (k, 1) = dvxc (k, 1) + dvxc_rr (k, 1, 1) * drho (k, 1) &
                + dvxc_sr (k, 1, 1) * s1
           DO ipol = 1, 3
              h (ipol, k, 1) = (dvxc_sr(k, 1, 1) * drho(k, 1) + &
                                dvxc_ss(k, 1, 1) * s1 )*grho(ipol, k, 1) + &
                                dvxc_s (k, 1, 1) * gdrho (ipol, k, 1)
           ENDDO
        ELSE
           DO ipol = 1, 3
              h (ipol, k, 1) = (0.d0, 0.d0)
           ENDDO
        ENDIF
     ELSE
        !
        !    LSDA case
        !
        ps (:,:) = (0.d0, 0.d0)
        DO is = 1, nspin
           DO js = 1, nspin
              DO ipol = 1, 3
                 ps1(ipol, is, js) = drho (k, is) * grho (ipol, k, js)
                 ps(is, js) = ps(is, js) + grho(ipol,k,is)*gdrho(ipol,k,js)
              ENDDO
              DO ks = 1, nspin
                 IF (is==js.and.js==ks) THEN
                    a (is, js, ks) = dvxc_sr (k, is, is)
                    c (is, js, ks) = dvxc_sr (k, is, is)
                 ELSE
                    IF (is==1) THEN
                       a (is, js, ks) = dvxc_sr (k, 1, 2)
                    ELSE
                       a (is, js, ks) = dvxc_sr (k, 2, 1)
                    ENDIF
                    IF (js==1) THEN
                       c (is, js, ks) = dvxc_sr (k, 1, 2)
                    ELSE
                       c (is, js, ks) = dvxc_sr (k, 2, 1)
                    ENDIF
                 ENDIF
                 DO ipol = 1, 3
                    ps2 (ipol, is, js, ks) = ps (is, js) * grho (ipol, k, ks)
                 ENDDO
                 DO ls = 1, nspin
                    IF (is==js.and.js==ks.and.ks==ls) THEN
                       b (is, js, ks, ls) = dvxc_ss (k, is, is)
                    ELSE
                       IF (is==1) THEN
                          b (is, js, ks, ls) = dvxc_ss (k, 1, 2)
                       ELSE
                          b (is, js, ks, ls) = dvxc_ss (k, 2, 1)
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        DO is = 1, nspin
           DO js = 1, nspin
              dvxc (k, is) = dvxc (k, is) + dvxc_rr (k, is, js) * drho (k, &
                   js)
              DO ipol = 1, 3
                 h (ipol, k, is) = h (ipol, k, is) + &
                      dvxc_s (k, is, js) * gdrho(ipol, k, js)
              ENDDO
              DO ks = 1, nspin
                 dvxc (k, is) = dvxc (k, is) + a (is, js, ks) * ps (js, ks)
                 DO ipol = 1, 3
                    h (ipol, k, is) = h (ipol, k, is) + &
                         c (is, js, ks) * ps1 (ipol, js, ks)
                 ENDDO
                 DO ls = 1, nspin
                    DO ipol = 1, 3
                       h (ipol, k, is) = h (ipol, k, is) + &
                            b (is, js, ks, ls) * ps2 (ipol, js, ks, ls)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  ! linear variation of the second term
  DO is = 1, nspin
     CALL fft_graddot (dfft, h (1, 1, is), g, dh)
     DO k = 1, dfft%nnr
        dvxc (k, is) = dvxc (k, is) - dh (k)
     ENDDO
  ENDDO
  DEALLOCATE (dh)
  DEALLOCATE (h)
  DEALLOCATE (gdrho)
  RETURN
END SUBROUTINE dgradcor1
