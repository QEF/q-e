!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
  MODULE exchange_correlation
!=----------------------------------------------------------------------------=!

        USE kinds, ONLY: dbl

        IMPLICIT NONE
        SAVE

        PRIVATE

! ... Gradient Correction & exchange and correlation

        LOGICAL :: tgc
        REAL(dbl), PARAMETER :: small_rho = 1.0d-10

        PUBLIC :: v2gc, exch_corr_energy
        PUBLIC :: exch_corr_setup, exch_corr_print_info
        PUBLIC :: tgc

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

      SUBROUTINE exch_corr_setup( xc_type )
        USE funct, ONLY: dft, which_dft, igcx, igcc
        CHARACTER(LEN = 20) :: xc_type
          dft = TRIM( xc_type )
          CALL which_dft( dft )
          tgc       = .FALSE.
          if( ( igcx > 0 ) .OR. ( igcc > 0 ) ) tgc = .TRUE.
        RETURN
      END SUBROUTINE exch_corr_setup

!=----------------------------------------------------------------------------=!
      SUBROUTINE exch_corr_print_info(iunit)
        USE funct, ONLY: dft, iexch, icorr, igcx, igcc
        INTEGER, INTENT(IN) :: iunit
        CHARACTER(LEN = 20) :: xcgc_type
        CHARACTER(LEN = 60) :: exch_info
        CHARACTER(LEN = 60) :: corr_info
        CHARACTER(LEN = 60) :: exgc_info
        CHARACTER(LEN = 60) :: cogc_info

        WRITE(iunit,800)

          ! ...     iexch => Exchange functional form 
          ! ...     icorr => Correlation functional form
          ! ...     igcx  => Gradient Correction to the Exchange potential
          ! ...     igcc  => Gradient Correction to the Correlation potential

          SELECT CASE ( iexch )
            CASE (0)
              exch_info = 'NONE'
            CASE (1)
              exch_info = 'SLATER'
            CASE (2)
              exch_info = 'SLATER (alpha=1)'
            CASE DEFAULT
              exch_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( icorr )
            CASE (0)
              corr_info = 'NONE'
            CASE (1)
              corr_info = 'PERDEW AND ZUNGER'
            CASE (2)
              corr_info = 'VOSKO, WILK AND NUSAIR'
            CASE (3)
              corr_info = 'LEE, YANG, AND PARR'
            CASE (4)
              corr_info = 'PERDEW AND WANG'
            CASE (9)
              corr_info = 'PADE APPROXIMATION'
            CASE DEFAULT
              corr_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( igcx )
            CASE (0)
              exgc_info = 'NONE'
            CASE (1)
              exgc_info = 'BECKE'
            CASE (2)
              exgc_info = 'PERDEW'
            CASE (3)
              exgc_info = 'PERDEW BURKE ERNZERHOF'
            CASE DEFAULT
              exgc_info = 'UNKNOWN'
          END SELECT
          SELECT CASE ( igcc )
            CASE (0)
              cogc_info = 'NONE'
            CASE (1)
              cogc_info = 'PERDEW'
            CASE (2)
              cogc_info = 'LEE, YANG AND PARR'
            CASE (3)
              cogc_info = 'PERDEW AND WANG'
            CASE (4)
              cogc_info = 'PERDEW BURKE ERNZERHOF'
            CASE DEFAULT
              cogc_info = 'UNKNOWN'
          END SELECT

          WRITE(iunit,910)
          WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exch_info
          WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') corr_info
          IF(tgc) THEN
            WRITE(iunit,810)
            WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exgc_info
            WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') cogc_info
          END IF

        WRITE( iunit, '(5x,"Exchange-correlation      = ",a, &
       &       " (",4i1,")")') trim(dft) , iexch, icorr, igcx, igcc

800 FORMAT(//,3X,'Exchange and correlations functionals',/ &
             ,3X,'-------------------------------------')
810 FORMAT(   3X,'Using Generalized Gradient Corrections with')
910 FORMAT(   3X,'Using Local Density Approximation with')

        RETURN
      END SUBROUTINE exch_corr_print_info

!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!

        SUBROUTINE v2gc(v2xc, grho, rhoer, vpot, gv)

          USE fft
          USE stick, ONLY: dfftp
          USE cell_base, ONLY: tpiba
          USE cp_types
          USE mp_global
!                                                                       
          implicit none
!                                                                       
          REAL(dbl) ::  vpot(:,:,:,:)
          REAL(dbl), intent(in)  ::  v2xc(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  grho(:,:,:,:,:)
          REAL(dbl), intent(in)  ::  rhoer(:,:,:,:)
          type (recvecs), intent(in)  ::  GV
!                                                                       
          integer ig, ipol, nxl, nyl, nzl, i, j, k, is, js, nspin
          COMPLEX(dbl), allocatable ::  psi(:,:,:)
          COMPLEX(dbl), allocatable ::  vtemp(:)
          COMPLEX(dbl), allocatable ::  vtemp_pol(:)
          REAL(dbl), ALLOCATABLE :: v(:,:,:)
          REAL(dbl) :: fac
          INTEGER :: kk(2,2), nn(2,2)
! ...                                                                   
          nxl   = dfftp%nr1
          nyl   = dfftp%nr2
          nzl   = dfftp%npl
          nspin = SIZE(rhoer,4)

          kk(1,1)=1
          kk(1,2)=3
          kk(2,1)=3
          kk(2,2)=2
          nn(1,1)=1
          nn(1,2)=2
          nn(2,1)=1
          nn(2,2)=2
          fac = REAL(nspin) / 2.0d0
          DO js = 1, nspin
            allocate(vtemp(gv%ng_l))
            allocate(vtemp_pol(gv%ng_l))
            allocate(psi(nxl,nyl,nzl))
            vtemp = CMPLX(0.0d0,0.0d0)

! 
!            DO is = 1, nspin
!              DO ipol = 1, 3
!                DO k = 1, nzl
!                  DO j = 1, nyl
!                    DO i = 1, nxl
!                      psi(i,j,k) = fac * v2xc(i,j,k,kk(js,is)) * grho(i,j,k,ipol,nn(js,is))
!                    END DO
!                  END DO
!                END DO
!                CALL pfwfft(vtemp_pol, psi)
!                DO ig = gv%gstart, gv%ng_l
!                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) * CMPLX( 0.d0, tpiba * gv%gx_l( ipol, ig ) )
!                END DO
!              END DO
!            END DO

            DO ipol = 1, 3
              DO is = 1, nspin
                DO k = 1, nzl
                  DO j = 1, nyl
                    DO i = 1, nxl
                      psi(i,j,k) = fac * v2xc(i,j,k,js,is) * grho(i,j,k,ipol,is)
                    END DO
                  END DO
                END DO
                CALL pfwfft(vtemp_pol, psi)
                DO ig = gv%gstart, gv%ng_l
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *  CMPLX( 0.d0, tpiba * gv%gx_l( ipol, ig ) )
                END DO
              END DO
            END DO
            DEALLOCATE(psi)
            DEALLOCATE(vtemp_pol)

            ALLOCATE(v(nxl,nyl,nzl))
            CALL pinvfft(v, vtemp)

            DO k = 1, nzl
              DO j = 1, nyl
                DO i = 1, nxl
                  vpot(i,j,k,js) = vpot(i,j,k,js) - v(i,j,k)
                END DO
              END DO
            END DO
            DEALLOCATE(vtemp)
            DEALLOCATE(v)
          END DO
          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!


        SUBROUTINE exch_corr_energy(rhoetr, rhoetg, grho, vpot, sxc, vxc, v2xc, gv)
          USE kinds, ONLY: dbl
          USE cp_types, ONLY: recvecs
          REAL (dbl) :: rhoetr(:,:,:,:)
          COMPLEX(dbl) :: rhoetg(:,:)
          REAL (dbl) :: grho(:,:,:,:,:)
          REAL (dbl) :: vpot(:,:,:,:)
          REAL(dbl) :: sxc              ! E_xc   energy
          REAL(dbl) :: vxc              ! SUM ( v(r) * rho(r) )
          REAL (dbl) :: v2xc(:,:,:,:,:)
          TYPE (recvecs), INTENT(IN) :: gv

          INTEGER    :: nxl, nyl, nzl, nspin, nnr

          ! ...     vpot = vxc(rhoetr); vpot(r) <-- u(r)

          nxl = SIZE(rhoetr,1)
          nyl = SIZE(rhoetr,2)
          nzl = SIZE(rhoetr,3)
          nspin = SIZE(rhoetr,4)
          nnr = nxl * nyl * nzl
          !
          CALL exch_corr_wrapper( nnr, nspin, grho, rhoetr, sxc, vpot, v2xc )
          !
          vxc = SUM( vpot * rhoetr )
          !
          IF(tgc) THEN
            ! ... vpot additional term for gradient correction
            CALL v2gc( v2xc, grho, rhoetr, vpot, gv )
          END If
          RETURN
        END SUBROUTINE

!=----------------------------------------------------------------------------=!
      END MODULE exchange_correlation
!=----------------------------------------------------------------------------=!
