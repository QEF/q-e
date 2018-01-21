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
     !CALL gradient1 (dfft, drhoc(1, is), g, gdrho (1,1,is) )
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
     !CALL grad_dot1 (dfft, h (1, 1, is), g, dh)
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
!
!--------------------------------------------------------------------
SUBROUTINE gradient1( dfft, a, g, ga)
  !--------------------------------------------------------------------
  ! Calculates ga = \grad a in R-space (a is G-space)
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE fft_interfaces, ONLY : fwfft, invfft
  USE fft_types,      ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  COMPLEX(DP) :: a (dfft%nnr)
  real(DP) :: ga (3, dfft%nnr), g (3, dfft%ngm)
  !
  INTEGER :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: gaux (:)

  ALLOCATE (gaux(dfft%nnr))

  ! a(G) multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  !  do ipol = 1, 3
  ! x, y
     ipol=1
     DO n = 1, dfft%nnr
        gaux (n) = (0.d0, 0.d0)
     ENDDO
     DO n = 1, dfft%ngm
        gaux(dfft%nl (n)) = cmplx(0.d0, g(ipol, n),kind=DP)* a (dfft%nl(n)) - &
                                        g(ipol+1, n) * a (dfft%nl(n))
        gaux(dfft%nlm(n)) = cmplx(0.d0,-g(ipol,n),kind=DP)* conjg(a (dfft%nl(n))) + &
                                     g(ipol+1, n) * conjg(a (dfft%nl(n)))
     ENDDO
     ! bring back to R-space, (\grad_ipol a)(r) ...

     CALL invfft ('Rho', gaux, dfft )
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     DO n = 1, dfft%nnr
        ga (ipol  , n) =  dble(gaux (n)) * tpiba
        ga (ipol+1, n) = aimag(gaux (n)) * tpiba
     ENDDO
  ! z
     ipol=3
     DO n = 1, dfft%nnr
        gaux (n) = (0.d0, 0.d0)
     ENDDO
     DO n = 1, dfft%ngm
        gaux(dfft%nl (n)) = cmplx(0.d0, g(ipol, n),kind=DP) * a (dfft%nl(n))
        gaux(dfft%nlm(n)) = conjg(gaux(dfft%nl(n)))
     ENDDO
     ! bring back to R-space, (\grad_ipol a)(r) ...
     CALL invfft ('Rho', gaux, dfft )
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     DO n = 1, dfft%nnr
        ga (ipol, n) =  dble(gaux (n)) * tpiba
     ENDDO
!  enddo
  DEALLOCATE (gaux)
  RETURN

END SUBROUTINE gradient1
!--------------------------------------------------------------------
SUBROUTINE grad_dot1 (dfft, a, g, da)
  !--------------------------------------------------------------------
  ! Calculates da = \sum_i \grad_i a_i in R-space
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE fft_interfaces, ONLY : fwfft, invfft
  USE fft_types,      ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  TYPE(fft_type_descriptor),INTENT(IN) :: dfft
  REAL(DP) :: a (3,dfft%nnr), da(dfft%nnr)
  real(DP) :: g (3, dfft%ngm)
  !
  INTEGER :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux (:), daux(:)
  COMPLEX(DP) :: fp, fm, aux1, aux2

  ALLOCATE (aux (dfft%nnr))
  ALLOCATE (daux (dfft%nnr))

  DO n = 1, dfft%nnr
     daux(n) = (0.d0, 0.d0)
  ENDDO
!!!  do ipol = 1, 3
     ! x, y
     ipol=1
     ! copy a(ipol,r) to a complex array...
     DO n = 1, dfft%nnr
        aux (n) = cmplx( dble(a(ipol, n)), dble(a(ipol+1, n)),kind=DP)
     ENDDO
     ! bring a(ipol,r) to G-space, a(G) ...
     CALL fwfft ('Rho', aux, dfft)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     DO n = 1, dfft%ngm
        fp = (aux(dfft%nl (n)) + aux (dfft%nlm(n)))*0.5d0
        fm = (aux(dfft%nl (n)) - aux (dfft%nlm(n)))*0.5d0
        aux1 = cmplx( dble(fp), aimag(fm),kind=DP)
        aux2 = cmplx(aimag(fp),- dble(fm),kind=DP)
        daux (dfft%nl(n)) = daux (dfft%nl(n)) + cmplx(0.d0,g(ipol,n),kind=DP)*aux1 + &
                                  cmplx(0.d0, g(ipol+1, n),kind=DP) * aux2
     ENDDO
     ! z
     ipol=3
     ! copy a(ipol,r) to a complex array...
     DO n = 1, dfft%nnr
        aux (n) = a(ipol, n)
     ENDDO
     ! bring a(ipol,r) to G-space, a(G) ...
     CALL fwfft ('Rho', aux, dfft)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     DO n = 1, dfft%ngm
        daux (dfft%nl(n)) = daux (dfft%nl(n)) + cmplx(0.d0, g(ipol,n),kind=DP) * &
             aux(dfft%nl(n))
     ENDDO
!!!  enddo
  DO n = 1, dfft%ngm
     daux(dfft%nlm(n)) = conjg(daux(dfft%nl(n)))
  ENDDO
  !  bring back to R-space, (\grad_ipol a)(r) ...
  CALL invfft ('Rho', daux, dfft )
  ! ...add the factor 2\pi/a  missing in the definition of q+G and sum
  DO n = 1, dfft%nnr
     da (n) = DBLE(daux(n)) * tpiba
  ENDDO
  DEALLOCATE (daux)
  DEALLOCATE (aux)

  RETURN
END SUBROUTINE grad_dot1
