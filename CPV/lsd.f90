!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      module local_spin_density

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

! ... Gradient Correction & exchange and correlation

        LOGICAL :: tgc
        CHARACTER(LEN = 20) :: xcgc_type
        CHARACTER(LEN = 60) :: exch_info
        CHARACTER(LEN = 60) :: corr_info
        CHARACTER(LEN = 60) :: exgc_info
        CHARACTER(LEN = 60) :: cogc_info
        INTEGER :: mfxcx, mfxcc, mgcx, mgcc
        REAL(dbl) :: salpha, bbeta, betapp
        REAL(dbl), PARAMETER :: small_rho = 1.0d-10

        INTERFACE vofxc_lsd
          MODULE PROCEDURE vofxc_lsd_gc, vofxc_lsd_nogc
        END INTERFACE

        PUBLIC :: tgc, v2gc, vofxc_lsd, lsd_setup, lsd_print_info

      contains


        SUBROUTINE lsd_setup(type_inp)
          CHARACTER(LEN = 20), intent(in) :: type_inp

          xcgc_type = type_inp
          tgc       = .FALSE.
          bbeta     = 0.0042d0
          betapp    = 0.0d0
          salpha    = 0.6666666666d0

! ...     mfxcx => Exchange functional form 
! ...     mfxcc => Correlation functional form
! ...     mgcx  => Gradient Correction to the Exchange potential
! ...     mgcc  => Gradient Correction to the Correlation potential

          IF(type_inp .EQ. 'BLYP') THEN
            mfxcx = 1
            mgcx  = 1
            mfxcc = 3
            mgcc  = 2
            tgc = .TRUE.
          ELSE IF(type_inp .EQ. 'BP') THEN
            mfxcx = 1
            mgcx  = 1
            mfxcc = 1
            mgcc  = 1
            tgc = .TRUE.
          ELSE IF(type_inp .EQ. 'PZ') THEN
            mfxcx = 1
            mfxcc = 1
            mgcx  = 0
            mgcc  = 0
          ELSE IF(type_inp .EQ. 'PW') THEN
            mfxcx = 1
            mfxcc = 4
            mgcx  = 0
            mgcc  = 0
          ELSE IF(type_inp .EQ. 'LDA') THEN
            mfxcx = 1
            mfxcc = 1
            mgcx  = 0
            mgcc  = 0
          ELSE IF(type_inp .EQ. 'PBE') THEN
            mfxcx = 1
            mfxcc = 1
            mgcx  = 3
            mgcc  = 4
            tgc = .TRUE.
          ELSE
            CALL errore(' module local spin density ', ' wrong xc type',0)
          END IF

          SELECT CASE ( mfxcx ) 
            CASE (1)
              exch_info = 'SLATER'
            CASE DEFAULT
              exch_info = 'UNKNOWN'
          END SELECT

          SELECT CASE ( mfxcc )
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

          SELECT CASE ( mgcx ) 
            CASE (1)
              exgc_info = 'BECKE'
            CASE (2)
              exgc_info = 'PERDEW'
            CASE (3)
              exgc_info = 'PERDEW BURKE ERNZERHOF'
            CASE DEFAULT
              exgc_info = 'UNKNOWN'
          END SELECT

          SELECT CASE ( mgcc ) 
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

          RETURN
        END SUBROUTINE lsd_setup

!=----------------------------------------------------------------------------=!

        SUBROUTINE lsd_print_info(iunit)
          INTEGER, INTENT(IN) :: iunit
          WRITE(iunit,910)
          WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exch_info
          WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') corr_info
          IF(tgc) THEN
            WRITE(iunit,810)
            WRITE(iunit,fmt='(5X,"Exchange functional: ",A)') exgc_info
            WRITE(iunit,fmt='(5X,"Correlation functional: ",A)') cogc_info
          END IF
810 FORMAT(   3X,'Using Generalized Gradient Corrections with')
910 FORMAT(   3X,'Using Local Density Approximation with')
          RETURN
        END SUBROUTINE lsd_print_info

!=----------------------------------------------------------------------------=!


        SUBROUTINE v2gc(v2xc, v2xc2, grho, rhoer, vpot, gv) 
                                                                        
          USE fft 
          USE stick, ONLY: dfftp
          USE cell_base, ONLY: tpiba
          USE cp_types 
          USE mp_global
!                                                                       
          implicit none 
!                                                                       
          REAL(dbl) ::  vpot(:,:,:,:) 
          REAL(dbl) ::  v2xc2(:,:,:) 
          REAL(dbl), intent(in)  ::  v2xc(:,:,:,:) 
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
          kk(2,2)=2
          kk(2,1)=3
          kk(1,2)=3
          nn(1,1)=1
          nn(2,2)=2
          nn(1,2)=2
          nn(2,1)=1                                                                
          fac = REAL(nspin) / 2.0_dbl


          DO js = 1, nspin
            allocate(vtemp(gv%ng_l)) 
            allocate(vtemp_pol(gv%ng_l)) 
            allocate(psi(nxl,nyl,nzl)) 
            vtemp = CMPLX(0.0d0,0.0d0)
            DO is = 1, nspin
              DO ipol = 1, 3 
                DO k = 1, nzl
                  DO j = 1, nyl
                    DO i = 1, nxl
                      psi(i,j,k) = fac * v2xc(i,j,k,kk(js,is)) * &
                        grho(i,j,k,ipol,nn(js,is)) 
                    END DO
                  END DO
                END DO
                CALL pfwfft(vtemp_pol, psi)                        
                DO ig = gv%gstart, gv%ng_l 
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *           &
     &                CMPLX( 0.d0, tpiba * gv%gx_l( ipol, ig ) )           
                END DO 
              END DO 
            END DO 
            DEALLOCATE(psi)
            DEALLOCATE(vtemp_pol)
                                                                        
            ALLOCATE(v(nxl,nyl,nzl)) 
            CALL pinvfft(v, vtemp) 
            DEALLOCATE(vtemp)

            DO k = 1, nzl
              DO j = 1, nyl
                DO i = 1, nxl
                  vpot(i,j,k,js) = vpot(i,j,k,js) - v(i,j,k)
                  v2xc2(i,j,k)   = v2xc2(i,j,k) + rhoer(i,j,k,js) * v(i,j,k)
                END DO
              END DO
            END DO
            DEALLOCATE(v)
          END DO
          RETURN 
        END SUBROUTINE                                         

        SUBROUTINE vofxc_lsd_gc(rhoe,vpot,grho,v2xc,v2xc2,sxc,vxc)
       
          REAL(dbl) :: rhoe(:,:,:,:)
          REAL(dbl) :: vpot(:,:,:,:)
          REAL(dbl) :: grho(:,:,:,:,:)
          REAL(dbl) :: v2xc(:,:,:,:)
          REAL(dbl) :: v2xc2(:,:,:)
          REAL(dbl) :: sxc, vxc

          INTEGER    ipol,i,j,k,nxl,nyl,nzl,is, nspin
          REAL(dbl)     roe,rhoa,rhob,v1x(2),v1c(2),v2c(3),v2x(3)
          REAL(dbl)     ex,ec,sx,sc,eta,grhoaa,grhobb,grhoab
          REAL(dbl)     VXA(2),VCA(2)

          sxc = 0.0d0
          vxc = 0.0d0
          nxl = SIZE(rhoe,1)
          nyl = SIZE(rhoe,2)
          nzl = SIZE(rhoe,3)
          nspin = SIZE(rhoe,4)

          IF(nspin .EQ. 1) THEN
            DO k = 1, nzl
              DO j = 1, nyl
                DO i = 1, nxl
                  ROE   = RHOE(i,j,k,1)
                  rhoa  = 0.5d0 * roe
                  rhob  = 0.5d0 * roe 
                  eta   = 0.0d0
  
                  IF ( roe .GT. small_rho ) THEN
                    CALL xc_lsd(roe,eta,ex,ec,vxa(1),vca(1),vxa(2),vca(2))
                  ELSE
                    ex = 0.0d0; ec = 0.0d0; vxa = 0.0d0; vca = 0.0d0;
                  END IF

                  SXC = SXC + (EX+EC)*ROE
                  VPOT(i,j,k,1)  =  VXA(1) + VCA(1)

                  IF( tgc ) THEN

                    v2xc2(i,j,k) =  ((EX+EC) - (VXA(1)+VCA(1)))*ROE

                    grhoaa =          grho(i,j,k,1,1)**2
                    grhoaa = grhoaa + grho(i,j,k,2,1)**2
                    grhoaa = grhoaa + grho(i,j,k,3,1)**2
                    grhoaa = grhoaa * 0.25d0
                    grhobb = grhoaa
                    grhoab = grhoaa

                    IF (roe .GT. small_rho) THEN
                      call GC_LSD(RHOA,RHOB,grhoaa,grhobb,grhoab,sx,sc, &
                      v1x(1),v2x(1),v1x(2),v2x(2),v1c(1),v2c(1),v1c(2), &
                      v2c(2),v2x(3),v2c(3))
                    ELSE
                      sx = 0.0d0
                      sc = 0.0d0
                      v1x = 0.0d0
                      v2x = 0.0d0
                      v1c = 0.0d0
                      v2c = 0.0d0
                    END IF

                    v2xc(i,j,k,1) = v2x(1) + &
                      2.0d0*(v2c(1)+v2c(2)+v2c(3)*2.d0)*0.25d0

                    VPOT(i,j,k,1) = VPOT(i,j,k,1) + v1x(1) + v1c(1)
                    sxc         = sxc + sx + sc

                    v2xc2(i,j,k)= v2xc2(i,j,k) + sx + sc &
                      - (v1x(1)+v1c(1))*rhoa - (v1x(2)+v1c(2))*rhob
                  END IF

                END DO
              END DO
            END DO

          ELSE

            DO k = 1, nzl
              DO j = 1, nyl
                DO i = 1, nxl

! ...             Exchange and Correlation Potentials
                  ROE   =   RHOE(i,j,k,1) + RHOE(i,j,k,2)
                  ETA   = ( RHOE(i,j,k,1) - RHOE(i,j,k,2) ) / ROE
                  !
                  !  rhoe(...1) = spin up   charge density
                  !  rhoe(...2) = spin down charge density
                  !
                  IF (roe .GT. small_rho) THEN
                    call xc_lsd(roe,eta,ex,ec,vxa(1),vca(1),vxa(2),vca(2))
                  ELSE
                    ex = 0.0d0; ec = 0.0d0; vxa = 0.0d0; vca = 0.0d0
                  END IF

                  vpot(i,j,k,1) = VXA(1) + VCA(1) 
                  vpot(i,j,k,2) = VXA(2) + VCA(2) 
                  SXC = SXC + (EX + EC) * ROE 

! ...             Gradient corrections to the Exchange and Correlation Potentials
                  IF ( tgc ) THEN
                    v2xc2(i,j,k) = (EX + EC) * ROE 
                    v2xc2(i,j,k) = v2xc2(i,j,k) - (VXA(1)+VCA(1))*RHOE(i,j,k,1)
                    v2xc2(i,j,k) = v2xc2(i,j,k) - (VXA(2)+VCA(2))*RHOE(i,j,k,2)

                    rhoa =  RHOE(i,j,k,1)
                    rhob =  RHOE(i,j,k,2)

                    grhoaa =          grho(i,j,k,1,1)**2
                    grhoaa = grhoaa + grho(i,j,k,2,1)**2
                    grhoaa = grhoaa + grho(i,j,k,3,1)**2

                    grhobb =          grho(i,j,k,1,2)**2
                    grhobb = grhobb + grho(i,j,k,2,2)**2
                    grhobb = grhobb + grho(i,j,k,3,2)**2

                    grhoab =          grho(i,j,k,1,1)* grho(i,j,k,1,2)
                    grhoab = grhoab + grho(i,j,k,2,1)* grho(i,j,k,2,2)
                    grhoab = grhoab + grho(i,j,k,3,1)* grho(i,j,k,3,2)

                    IF(roe .GT. small_rho) THEN
                      call GC_LSD(rhoa, rhob, grhoaa, grhobb, grhoab, sx, sc, &
                        v1x(1), v2x(1), v1x(2), v2x(2), v1c(1), v2c(1), v1c(2), &
                        v2c(2), v2x(3), v2c(3))
                    ELSE
                      sx  = 0.0d0; sc  = 0.0d0; v1x = 0.0d0
                      v2x = 0.0d0; v1c = 0.0d0; v2c = 0.0d0
                    END IF

                    v2xc(i,j,k,1) = v2x(1) + v2c(1) 
                    v2xc(i,j,k,2) = v2x(2) + v2c(2) 
                    v2xc(i,j,k,3) = v2x(3) + v2c(3) 

                    vpot(i,j,k,1) = vpot(i,j,k,1) + V1X(1) + V1C(1)  
                    vpot(i,j,k,2) = vpot(i,j,k,2) + V1X(2) + V1C(2)  

                    SXC = SXC + sx + sc 

                    v2xc2(i,j,k) = v2xc2(i,j,k) + sx + sc 
                    v2xc2(i,j,k) = v2xc2(i,j,k) - (V1X(1)+V1C(1)) * RHOA
                    v2xc2(i,j,k) = v2xc2(i,j,k) - (V1X(2)+V1C(2)) * RHOB

                  END IF

                END DO
              END DO
            END DO

          END IF

          RETURN
        END SUBROUTINE vofxc_lsd_gc


        SUBROUTINE vofxc_lsd_nogc(rhoe, vpot, sxc, vxc)
       
          REAL(dbl) :: rhoe(:,:,:,:)
          REAL(dbl) :: vpot(:,:,:,:)
          REAL(dbl) sxc,vxc

          integer    ipol,i,j,k,nxl,nyl,nzl,is, nspin
          REAL(dbl)     roe,rhoa,rhob,v1x(2),v1c(2),v2c(3),v2x(3)
          REAL(dbl)     ex,ec,sx,sc,eta,grhoaa,grhobb,grhoab
          REAL(dbl)     VXA(2),VCA(2)

          sxc = 0.0d0
          vxc = 0.0d0
          nxl = SIZE(rhoe,1)
          nyl = SIZE(rhoe,2)
          nzl = SIZE(rhoe,3)
          nspin = SIZE(rhoe,4)

          IF(nspin .EQ. 1) THEN
            do k=1,nzl
              do j=1,nyl
                do i=1,nxl
                  ROE   = RHOE(i,j,k,1)
                  eta   = 0.0d0
                  if(roe.gt.small_rho) then
                    call xc_lsd(roe,eta,ex,ec,vxa(1),vca(1),vxa(2),vca(2))
                  else
                    ex = 0.0d0; ec = 0.0d0; vxa = 0.0d0; vca = 0.0d0;
                  endif
                  VPOT(i,j,k,1)  =  VXA(1) + VCA(1)
                  SXC = SXC + (EX+EC)*ROE
                end do
              end do
            end do
          ELSE
            do k=1,nzl
              do j=1,nyl
                do i=1,nxl
                  ROE   =       RHOE(i,j,k,1)
                  ROE   = ROE + RHOE(i,j,k,2)
                  ETA   = ( RHOE(i,j,k,1) - RHOE(i,j,k,2) ) / ROE
                  IF (roe.GT.small_rho) THEN
                    call xc_lsd(roe,eta,ex,ec,vxa(1),vca(1),vxa(2),vca(2))
                  ELSE
                    ex = 0.0_dbl; ec = 0.0_dbl; vxa = 0.0_dbl; vca = 0.0_dbl;
                  ENDIF
                  vpot(i,j,k,1) = VXA(1) + VCA(1) 
                  vpot(i,j,k,2) = VXA(2) + VCA(2) 
                  SXC = SXC + (EX + EC) * ROE 
                end do
              end do
            end do
          END IF
          RETURN
        END SUBROUTINE vofxc_lsd_nogc


!     ==================================================================
      SUBROUTINE XC_LSD(RHO,ETA,EX,EC,VXA,VCA,VXB,VCB)
!     ==--------------------------------------------------------------==
!     ==  LSD EXCHANGE AND CORRELATION FUNCTIONALS                    ==
!     ==                                                              ==
!     ==  EXCHANGE  :  SLATER alpha                                   ==
!     ==  CORRELATION : PERDEW & ZUNGER                               ==
!     ==                VOSKO, WILK & NUSSAIR                         ==
!     ==                LEE, YANG & PARR                              ==
!     ==                PERDEW & WANG                                 ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)

! ... Gradient Correction & exchange and correlation
      REAL(dbl) :: SMALL, PI34, THIRD
      PARAMETER (SMALL=1.D-20)
      PARAMETER (PI34= 0.75D0 / 3.141592653589793D+00,THIRD=1.D0/3.D0)

!     ==--------------------------------------------------------------==
! ... Exchange Options
!     ==--------------------------------------------------------------==
      IF(MFXCX.EQ.1 .AND. RHO.GT.SMALL) THEN
        CALL LSD_SX(RHO,ETA,EX,VXA,VXB,SALPHA)
      ELSE
        EX=0.0D0
        VXA=0.0D0
        VXB=0.0D0
      ENDIF
!     ==--------------------------------------------------------------==
! ... Correlation Options
!     ==--------------------------------------------------------------==
      IF(RHO.LE.SMALL) THEN
        EC  = 0.0D0
        VCA = 0.0D0
        VCB = 0.0D0
        EX  = 0.0D0
        VXA = 0.0D0
        VXB = 0.0D0
      ELSE IF(MFXCC.EQ.1) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL LSD_PZ(RS,ETA,EC,VCA,VCB,IFLG)
      ELSEIF(MFXCC.EQ.2) THEN
        RS = (PI34/RHO)**THIRD
        CALL LSD_VWN(RS,EC,VC)
      ELSEIF(MFXCC.EQ.3) THEN
        CALL LSD_LYP(RHO,ETA,EC,VCA,VCB)
      ELSEIF(MFXCC.EQ.4) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL LSD_PW(RS,EC,VC,IFLG)
        VCA = VC
      ELSEIF(MFXCC.EQ.9) THEN
        CALL LSD_PADE(RHO,ETA,EC,VCA,VCB)
      ELSE
        EC=0.0D0
        VCA=0.0D0
        VCB=0.0D0
      ENDIF
      RETURN
      END SUBROUTINE XC_LSD
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE GC_LSD(RHOA,RHOB,GRHOAA,GRHOBB,GRHOAB,SX,SC,V1XA,V2XA, &
                        V1XB,V2XB,V1CA,V2CA,V1CB,V2CB,V2XAB,V2CAB)
!     ==--------------------------------------------------------------==
!     ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
!     ==                                                              ==
!     ==  EXCHANGE  :  BECKE88                                        ==
!     ==               GGAX                                           ==
!     ==  CORRELATION : PERDEW86                                      ==
!     ==                LEE, YANG & PARR                              ==
!     ==                GGAC                                          ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)

! ... Gradient Correction & exchange and correlation

      REAL(dbl) :: SMALL
      PARAMETER(SMALL=1.D-20)
!     ==--------------------------------------------------------------==
! ... Exchange Options
!     ==--------------------------------------------------------------==
      IF(MGCX.EQ.1) THEN
        CALL LSD_B88(BBETA,RHOA,RHOB,GRHOAA,GRHOBB, SX,V1XA,V2XA,V1XB,V2XB)
        V2XAB=0.0D0
      ELSEIF(MGCX.EQ.2) THEN
        CALL LSD_GGAX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.3) THEN
        CALL LSD_PBEX(RHOA,RHOB,GRHOAA,GRHOBB,SX,V1XA,V2XA,V1XB,V2XB)
        V2XAB=0.0D0
      ELSE
        SX=0.0D0
        V1XA=0.0D0
        V2XA=0.0D0
        V1XB=0.0D0
        V2XB=0.0D0
        V2XAB=0.0D0
      ENDIF
!     ==--------------------------------------------------------------==
! ... Correlation Options
!     ==--------------------------------------------------------------==
      SC=0.0D0
      V1CA=0.0D0
      V2CA=0.0D0
      V1CB=0.0D0
      V2CB=0.0D0
      V2CAB=0.0D0
      IF(MGCC.EQ.1) THEN
        RHO=RHOA+RHOB
        GRHO=GRHOAA+2.*GRHOAB+GRHOBB
        IF(RHO.GT.SMALL) THEN
          CALL LSD_P86(RHOA,RHOB,GRHO,SC,V1CA,V2CA,V1CB,V2CB,V2CAB)
        ENDIF
      ELSEIF(MGCC.EQ.2) THEN
        RHO=RHOA+RHOB
        IF(RHO.GT.SMALL)  CALL LSD_GLYP(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC,&
                  V1CA,V2CA,V1CB,V2CB,V2CAB)
      ELSEIF(MGCC.EQ.3) THEN
        CALL LSD_GGAC(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.4) THEN
        RHO=RHOA+RHOB
        IF(RHO.GT.SMALL)  CALL LSD_PBEC(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC, &
                        V1CA,V2CA,V1CB,V2CB,V2CAB)
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE GC_LSD




!     ==================================================================
      SUBROUTINE LSD_SX(RHO,ETA,EX,VXA,VXB,ALPHA)
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: F, THIRD, F43
      PARAMETER (F = -1.39578858466194911D0)
      PARAMETER (THIRD=1.D0/3.D0,F43=4.D0/3.D0)
!     ==--------------------------------------------------------------==
      RHOA=0.5D0*RHO*(1.D0+ETA)
      RHOB=0.5D0*RHO*(1.D0-ETA)
      RSA = RHOA**THIRD
      RSB = RHOB**THIRD
      EX = F*ALPHA*(RHOA*RSA+RHOB*RSB)/RHO
      VXA = F43*F*ALPHA*RSA
      VXB = F43*F*ALPHA*RSB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_SX



!     ==================================================================
      SUBROUTINE LSD_PZ(RS,ETA,EPZ,VPZA,VPZB,IFLG)
!     ==--------------------------------------------------------------==
!     ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: au, bu, cu, du, gcu, b1u, b2u
      REAL(dbl) :: ap, bp, cp, dp, gcp, b1p, b2p
      REAL(dbl) :: f13, f43
      PARAMETER (AU=0.0311D0,BU=-0.048D0,CU=0.0020D0,DU=-0.0116D0, &
                 GCU=-0.1423D0,B1U=1.0529D0,B2U=0.3334D0)
      PARAMETER (AP=0.01555D0,BP=-0.0269D0,CP=0.0007D0,DP=-0.0048D0, &
                 GCP=-0.0843D0,B1P=1.3981D0,B2P=0.2611D0)
      PARAMETER (F13=1.D0/3.D0,F43=4.D0/3.D0)
!     ==--------------------------------------------------------------==
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EU=AU*XLN+BU+CU*RS*XLN+DU*RS
        EP=AP*XLN+BP+CP*RS*XLN+DP*RS
        VU=AU*XLN+(BU-AU/3.D0)+2.D0/3.D0*CU*RS*XLN+ (2.D0*DU-CU)/3.D0*RS
        VP=AP*XLN+(BP-AP/3.D0)+2.D0/3.D0*CP*RS*XLN+ (2.D0*DP-CP)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OXU=1.D0+B1U*RS1+B2U*RS2
        DOXU=1.D0+7./6.*B1U*RS1+4./3.*B2U*RS2
        EU=GCU/OXU
        VU=EU*DOXU/OXU
        OXP=1.D0+B1P*RS1+B2P*RS2
        DOXP=1.D0+7./6.*B1P*RS1+4./3.*B2P*RS2
        EP=GCP/OXP
        VP=EP*DOXP/OXP
      ENDIF
      FETA=((1.+ETA)**F43+(1.-ETA)**F43-2.)/(2.**F43-2.)
      DFETA=F43*((1.+ETA)**F13-(1.-ETA)**F13)/(2.**F43-2.)
      EPZ=EU+FETA*(EP-EU)
      VPZA=VU+FETA*(VP-VU)+(EP-EU)*(1.-ETA)*DFETA
      VPZB=VU+FETA*(VP-VU)+(EP-EU)*(-1.-ETA)*DFETA
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_PZ



!     ==================================================================
      SUBROUTINE LSD_VWN(RS,EVWN,VVWN)
!     ==--------------------------------------------------------------==
!     ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
!     ==                 CAN. J. PHYS. 58 1200  (1980)                ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: a, b, c, x0, two
      PARAMETER (A=0.0310907,B=3.72744,C=12.9352,X0=-0.10498)
      PARAMETER (TWO=2.0D0)
!     ==--------------------------------------------------------------==
!      CALL STOPGM('LSD_VWN','NOT PROGRAMMED')
       CALL errore('LSD_VWN','NOT PROGRAMMED',-1)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_VWN



!     ==================================================================
      SUBROUTINE LSD_LYP(RHO,ETA,ELYP,VALYP,VBLYP)
!     ==--------------------------------------------------------------==
!     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
!     ==  THIS IS ONLY THE LDA PART                                   ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: a, b, c, d, cf, small
      PARAMETER (SMALL=1.D-24)
      PARAMETER (A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
      PARAMETER (CF=2.87123400018819108D0)
!     ==--------------------------------------------------------------==
      RA=RHO*0.5D0*(1.D0+ETA)
      RA=MAX(RA,SMALL)
      RB=RHO*0.5D0*(1.D0-ETA)
      RB=MAX(RB,SMALL)
      RM3=RHO**(-1.D0/3.D0)
      DR=(1.D0+D*RM3)
      E1=4.D0*A*RA*RB/RHO/DR
      OR=EXP(-C*RM3)/DR*RM3**11.D0
      DOR=-1.D0/3.D0*RM3**4*OR*(11.D0/RM3-C-D/DR)
      E2=2.D0**(11.D0/3.D0)*CF*A*B*OR*RA*RB*(RA**(8.d0/3.d0)+ RB**(8.d0/3.d0))
      ELYP=(-E1-E2)/RHO
      DE1A=-E1*(1.D0/3.D0*D*RM3**4/DR+1./RA-1./RHO)
      DE1B=-E1*(1.D0/3.D0*D*RM3**4/DR+1./RB-1./RHO)
      DE2A=-2.D0**(11.D0/3.D0)*CF*A*B*(DOR*RA*RB*(RA**(8.d0/3.d0)+ &
            RB**(8.d0/3.d0))+OR*RB*(11./3.*RA**(8.d0/3.d0)+ &
            RB**(8.d0/3.d0)))
      DE2B=-2.D0**(11.D0/3.D0)*CF*A*B*(DOR*RA*RB*(RA**(8.d0/3.d0)+ &
      RB**(8.d0/3.d0))+OR*RA*(11./3.*RB**(8.d0/3.d0)+ &
      RA**(8.d0/3.d0)))
      VALYP=DE1A+DE2A
      VBLYP=DE1B+DE2B
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_LYP



!     ==================================================================
      SUBROUTINE LSD_PADE(RHO,ETA,EC,VCA,VCB)
!     ==--------------------------------------------------------------==
!     ==  PADE APPROXIMATION                                          ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: a0, a1, a2, a3, b1, b2, b3, b4
      REAL(dbl) :: da0, da1, da2, da3, db1, db2, db3, db4
      REAL(dbl) :: rsfac, fsfac
      PARAMETER (A0=.4581652932831429d0,A1=2.217058676663745d0, &
           A2=0.7405551735357053d0,A3=0.01968227878617998d0)
      PARAMETER (B1=1.0D0,B2=4.504130959426697d0, &
                 B3=1.110667363742916d0,B4=0.02359291751427506d0)
      PARAMETER (DA0=.119086804055547D0,DA1=.6157402568883345d0, &
                 DA2=.1574201515892867d0,DA3=.003532336663397157d0)
      PARAMETER (DB1=0.0d0,DB2=.2673612973836267d0,  &
                 DB3=.2052004607777787d0,DB4=.004200005045691381d0)
      PARAMETER (RSFAC=.6203504908994000d0)
      PARAMETER (FSFAC=1.92366105093153617d0)
!     ==--------------------------------------------------------------==
      RS=RSFAC*RHO**(-1.d0/3.d0)
      FS=FSFAC*((1.d0+ETA)**(4.d0/3.d0)+(1.d0-ETA)**(4.d0/3.d0)-2.d0)
      DFS=FSFAC*4.d0/3.d0* ((1.d0+ETA)**(1.d0/3.d0)-(1.d0-ETA)**(1.d0/3.d0))
      DFSA=DFS*(1.d0-ETA)
      DFSB=DFS*(-1.d0-ETA)
      A0P=A0+FS*DA0
      A1P=A1+FS*DA1
      A2P=A2+FS*DA2
      A3P=A3+FS*DA3
      B1P=B1+FS*DB1
      B2P=B2+FS*DB2
      B3P=B3+FS*DB3
      B4P=B4+FS*DB4
      TOP=A0P+RS*(A1P+RS*(A2P+RS*A3P))
      DTOP=A1P+RS*(2.d0*A2P+RS*3.d0*A3P)
      TOPX=DA0+RS*(DA1+RS*(DA2+RS*DA3))
      BOT=RS*(B1P+RS*(B2P+RS*(B3P+RS*B4P)))
      DBOT=B1P+RS*(2.d0*B2P+RS*(3.d0*B3P+RS*4.d0*B4P))
      BOTX=RS*(DB1+RS*(DB2+RS*(DB3+RS*DB4)))
      EC=-TOP/BOT
      VC=EC+RS*(DTOP/BOT-TOP*DBOT/(BOT*BOT))/3.d0
      DX=-(TOPX/BOT-TOP*BOTX/(BOT*BOT))
      VCA=VC+DX*DFSA
      VCB=VC+DX*DFSB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_PADE


!     ==================================================================
      SUBROUTINE LSD_PW(RS,EPWC,VPWC,IFLG)
!     ==--------------------------------------------------------------==
!     ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z)
      INTEGER :: IFLG
      REAL(dbl) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
      PARAMETER (A=0.031091,A1=0.21370,B1=7.5957,B2=3.5876,B3=1.6382, &
        B4=0.49294,C0=A,C1=0.046644,C2=0.00664,C3=0.01043, D0=0.4335,D1=1.4408)
!     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-(2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1+A1*RS)*OLOG
        VPWC=-2.*A*(1.D0+2./3.*A1*RS)*OLOG-2./3.*A*(1.+A1*RS)*DOM/(OM*(OM+1.))
      ELSEIF(IFLG.EQ.3) THEN
!..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_PW


!     ==================================================================
      SUBROUTINE LSD_B88(B1,RHOA,RHOB,GRHOA,GRHOB, SX,V1XA,V2XA,V1XB,V2XB)
!     ==--------------------------------------------------------------==
! BECKE EXCHANGE: PRA 38, 3098 (1988)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: OB3, SMALL
      PARAMETER(OB3=1.D0/3.D0,SMALL=1.D-20)
!     ==--------------------------------------------------------------==
      SX=0.0D0
      V1XA=0.0D0
      V2XA=0.0D0
      V1XB=0.0D0
      V2XB=0.0D0
      IF(ABS(RHOA).GT.SMALL) THEN
        AA    = GRHOA
        A     = SQRT(AA)
        BR1   = RHOA**OB3
        BR2   = BR1*BR1
        BR4   = BR2*BR2
        XS    = A/BR4
        XS2   = XS*XS
        SA2B8 = SQRT(1.0D0+XS2)
        SHM1  = LOG(XS+SA2B8)
        DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
        DD2   = DD*DD
        DDD   = 6.0D0*B1*(SHM1+XS/SA2B8)
        GF    = -B1*XS2/DD
        DGF   = (-2.0D0*B1*XS*DD + B1*XS2*DDD)/DD2
        SX    = GF*BR4
        V1XA  = 4./3.*BR1*(GF-XS*DGF)
        V2XA  = DGF/A
      ENDIF
      IF(ABS(RHOB).GT.SMALL) THEN
        AA    = GRHOB
        A     = SQRT(AA)
        BR1   = RHOB**OB3
        BR2   = BR1*BR1
        BR4   = BR2*BR2
        XS    = A/BR4
        XS2   = XS*XS
        SA2B8 = SQRT(1.0D0+XS2)
        SHM1  = LOG(XS+SA2B8)
        DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
        DD2   = DD*DD
        DDD   = 6.0D0*B1*(SHM1+XS/SA2B8)
        GF    = -B1*XS2/DD
        DGF   = (-2.0D0*B1*XS*DD + B1*XS2*DDD)/DD2
        SX    = SX+GF*BR4
        V1XB  = 4./3.*BR1*(GF-XS*DGF)
        V2XB  = DGF/A
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_B88


!     ==================================================================
      SUBROUTINE LSD_GGAX(RHO,GRHO,SX,V1X,V2X)
!     ==--------------------------------------------------------------==
! J.P.PERDEW ET AL. PRB 46 6671 (1992)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: f1, f2, f3, f4, f5, pi
      PARAMETER(F1=0.19645,F2=7.7956,F3=0.2743,F4=0.1508,F5=0.004)
      PARAMETER(PI=3.141592653589793D+00)
!     ==--------------------------------------------------------------==
!     CALL STOPGM('LSD_GGAX','NOT PROGRAMMED')
      CALL errore('LSD_GGAX','NOT PROGRAMMED',-1)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_GGAX


!     ==================================================================
      SUBROUTINE LSD_PBEX(RHOA,RHOB,GRHOA,GRHOB,SX,V1XA,V2XA,V1XB,V2XB)
! J.P.PERDEW ET AL. PRL XX XXXX (1996)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: small, us, ax, um, uk, ul
      PARAMETER(SMALL=1.D-20)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0, &
                UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
!     ==--------------------------------------------------------------==
      SXA=0.0D0
      SXB=0.0D0
      V1XA=0.0D0
      V2XA=0.0D0
      V1XB=0.0D0
      V2XB=0.0D0
      IF(ABS(RHOA).GT.SMALL) THEN
        RHO   = 2.D0*RHOA
        AA    = 4.D0*GRHOA
        RR    = RHO**(-4./3.)
        EX    = AX/RR
        S2    = AA*RR*RR*US*US
        PO    = 1.D0/(1.D0 + UL*S2)
        FX    = UK-UK*PO
        SXA   = EX*FX
        DFX   = 2.D0*UK*UL*PO*PO
        V1XA  = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
        V2XA  = EX*DFX*(US*RR)**2
      ENDIF
      IF(ABS(RHOB).GT.SMALL) THEN
        RHO   = 2.D0*RHOB
        AA    = 4.D0*GRHOB
        RR    = RHO**(-4./3.)
        EX    = AX/RR
        S2    = AA*RR*RR*US*US
        PO    = 1.D0/(1.D0 + UL*S2)
        FX    = UK-UK*PO
        SXB   = EX*FX
        DFX   = 2.D0*UK*UL*PO*PO
        V1XB  = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
        V2XB  = EX*DFX*(US*RR)**2
      ENDIF
      SX    = 0.5*(SXA+SXB)
      V2XA  = 2.D0*V2XA
      V2XB  = 2.D0*V2XB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_PBEX


!     ==================================================================
      SUBROUTINE LSD_P86(RHOA,RHOB,GRHO,SC,V1CA,V2CA,V1CB,V2CB,V2CAB)
!     ==--------------------------------------------------------------==
! PERDEW CORRELATION: PRB 33, 8822 (1986)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: p1, p2, p3, p4, pc1, pc2, pci, ob3, fpi
      PARAMETER(P1=0.023266D0,P2=7.389D-6,P3=8.723D0,P4=0.472D0)
      PARAMETER(PC1=0.001667D0,PC2=0.002568,PCI=PC1+PC2)
      PARAMETER(OB3=1.D0/3.D0, FPI=4.0D0*3.141592653589793D0)
!     ==--------------------------------------------------------------==
      RHO=RHOA+RHOB
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      RS    = (3./(FPI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      CNA   = PC2+P1*RS+P2*RS2
      CNB   = 1.+P3*RS+P4*RS2+1.D4*P2*RS3
      CN    = PC1 + CNA/CNB
      DRS   = -OB3*(3./FPI)**OB3 / BR4
      DCNA  = (P1+2.*P2*RS)*DRS
      DCNB  = (P3+2.*P4*RS+3.D4*P2*RS2)*DRS
      DCN   = DCNA/CNB - CNA/(CNB*CNB)*DCNB
      S1    = SQRT(RHOA**(5./3.)+RHOB**(5./3.))
      S2    = 2.**OB3/RHO**(5./6.)
      D     = S1*S2
      DDA   = S2*5./6.*(RHOA**(2./3.)/S1-S1/RHO)
      DDB   = S2*5./6.*(RHOB**(2./3.)/S1-S1/RHO)
      PHI   = 0.192*PCI/CN*A*RHO**(-7./6.)
      EPHI  = EXP(-PHI)
      SC    = AA/BR4*CN*EPHI/D
      V1C   = SC*((1.+PHI)*DCN/CN -((4./3.)-(7./6.)*PHI)/RHO)
      V1CA  = -SC/D*DDA+V1C
      V1CB  = -SC/D*DDB+V1C
      V2C   = CN*EPHI/BR4*(2.-PHI)/D
      V2CA  = V2C
      V2CB  = V2C
      V2CAB = V2C
!     ==--------------------------------------------------------------==
      RETURN
      END  SUBROUTINE LSD_P86



!     ==================================================================
      SUBROUTINE LSD_GLYP(RA,RB,GRHOAA,GRHOAB,GRHOBB,SC,  &
                                            V1CA,V2CA,V1CB,V2CB,V2CAB)
!     ==--------------------------------------------------------------==
! LEE, YANG PARR: GRADIENT CORRECTION PART
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: a, b, c, d
      PARAMETER(A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
!     ==--------------------------------------------------------------==
      RHO=RA+RB
      RM3=RHO**(-1.D0/3.D0)
      DR=(1.D0+D*RM3)
      OR=EXP(-C*RM3)/DR*RM3**11.D0
      DOR=-1.D0/3.D0*RM3**4*OR*(11.D0/RM3-C-D/DR)
      DER=C*RM3+D*RM3/DR
      DDER=1./3.*(D*D*RM3**5/DR/DR-DER/RHO)
      DLAA=-A*B*OR*(RA*RB/9.*(1.-3*DER-(DER-11.)*RA/RHO)-RB*RB)
      DLAB=-A*B*OR*(RA*RB/9.*(47.-7.*DER)-4./3.*RHO*RHO)
      DLBB=-A*B*OR*(RA*RB/9.*(1.-3*DER-(DER-11.)*RB/RHO)-RA*RA)
      DLAAA=DOR/OR*DLAA-A*B*OR*(RB/9.*(1.-3*DER-(DER-11.)*RA/RHO)- &
            RA*RB/9.*((3.+RA/RHO)*DDER+(DER-11.)*RB/RHO/RHO))
      DLAAB=DOR/OR*DLAA-A*B*OR*(RA/9.*(1.-3.*DER-(DER-11.)*RA/RHO)- &
            RA*RB/9.*((3.+RA/RHO)*DDER-(DER-11.)*RA/RHO/RHO)-2.*RB)
      DLABA=DOR/OR*DLAB-A*B*OR*(RB/9.*(47.-7.*DER)-7./9.*RA*RB*DDER- &
            8./3.*RHO)
      DLABB=DOR/OR*DLAB-A*B*OR*(RA/9.*(47.-7.*DER)-7./9.*RA*RB*DDER- &
            8./3.*RHO)
      DLBBA=DOR/OR*DLBB-A*B*OR*(RB/9.*(1.-3.*DER-(DER-11.)*RB/RHO)- &
            RA*RB/9.*((3.+RB/RHO)*DDER-(DER-11.)*RB/RHO/RHO)-2.*RA)
      DLBBB=DOR/OR*DLBB-A*B*OR*(RA/9.*(1.-3*DER-(DER-11.)*RB/RHO)- &
            RA*RB/9.*((3.+RB/RHO)*DDER+(DER-11.)*RA/RHO/RHO))
      SC=DLAA*GRHOAA+DLAB*GRHOAB+DLBB*GRHOBB
      V1CA=DLAAA*GRHOAA+DLABA*GRHOAB+DLBBA*GRHOBB
      V1CB=DLAAB*GRHOAA+DLABB*GRHOAB+DLBBB*GRHOBB
      V2CA=2.*DLAA
      V2CB=2.*DLBB
      V2CAB=DLAB
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_GLYP


!     ==================================================================
      SUBROUTINE LSD_GGAC(RHO,GRHO,SC,V1C,V2C)
!     ==--------------------------------------------------------------==
! PERDEW & WANG GGA CORRELATION PART
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: al, pa, pb, pc, pd, cx, cxc0, cc0, ob3, pi
      PARAMETER(al=0.09,pa=0.023266D0,pb=7.389D-6,pc=8.723D0,pd=0.472D0)
      PARAMETER(CX=-0.001667D0,CXC0=0.002568,CC0=-CX+CXC0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
!     ==--------------------------------------------------------------==
!      CALL STOPGM('LSD_GGAC','NOT PROGRAMMED')
      CALL errore('LSD_GGAC','NOT PROGRAMMED',-1)
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE LSD_GGAC


!     ==================================================================
      SUBROUTINE LSD_PBEC(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC, &
                    V1CA,V2CA,V1CB,V2CB,V2CAB)
! PBE Correlation functional
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      REAL(dbl) :: be, ga, ob3, pi
      PARAMETER(BE=0.06672455060314922D0,GA=0.031090690869654895D0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
!     ==--------------------------------------------------------------==
      RHO=RHOA+RHOB
      ETA=(RHOA-RHOB)/RHO
      PHI=0.5D0*((1.d0+ETA)**(2./3.)+(1.d0-ETA)**(2./3.))
      PHI3=PHI*PHI*PHI
      GRHO=GRHOAA+2.D0*GRHOAB+GRHOBB
      CALL XC_LSD(RHO,ETA,EX,EC,VXA,VCA,VXB,VCB)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3./(4.*PI*RHO))**OB3
      XKF   = (9.*PI/4.)**OB3/RS
      XKS   = SQRT(4.*XKF/PI)
      T     = A/(2.*XKS*RHO*PHI)
      EXPE  = EXP(-EC/(PHI3*GA))
      AF    = BE/GA * (1./(EXPE-1.))
      Y     = AF*T*T
      XY    = (1.+Y)/(1.+Y+Y*Y)
      S1    = 1.+BE/GA*T*T*XY
      H0    = GA*PHI3 * LOG(S1)
      SC    = RHO*H0
!..
      ! DTDPHI= -T/PHI
      ! DPHIDE= 1./(3.*(1.d0+ETA)**OB3) - 1./(3.*(1.d0-ETA)**OB3)
      ! DEDRA = 2.*RHOB/(RHO*RHO)
      ! DEDRB = -2.*RHOA/(RHO*RHO)
      ! DPHIDA= DPHIDE*DEDRA
      ! DPHIDB= DPHIDE*DEDRB
      ! DTDRA = -T*(DPHIDA/PHI+7./(6.*RHO))
      ! DTDRB = -T*(DPHIDB/PHI+7./(6.*RHO))
      ! DADRA = AF*AF*EXPE/(-BE*PHI3)*(3.*EC/PHI*DPHIDA-(VCA-EC)/RHO)
      ! DADRB = AF*AF*EXPE/(-BE*PHI3)*(3.*EC/PHI*DPHIDB-(VCB-EC)/RHO)
      ! DSDA  = -BE/GA * AF * T**6 * (2.+Y) / (1.+Y+Y*Y)**2
      ! DSDT  = 2.*BE/GA * T * (1.+2.*Y) / (1.+Y+Y*Y)**2
      ! DSDRA = DSDA*DADRA + DSDT*DTDRA
      ! DSDRB = DSDA*DADRB + DSDT*DTDRB
      ! DHDT  = GA*PHI3/S1*DSDT
      ! DHDRA = 3.*H0/PHI*DPHIDA + GA*PHI3/S1*DSDRA
      ! DHDRB = 3.*H0/PHI*DPHIDB + GA*PHI3/S1*DSDRB

       if (ETA.eq.1.D0) THEN !if SIC then dPHI/dRHO=0
         DTDPHI= -T/PHI
         DEDRA = 2.*RHOB/(RHO*RHO)
         DEDRB = -2.*RHOA/(RHO*RHO)
         DTDRA = -T*7./(6.*RHO)
         DTDRB = -T*7./(6.*RHO)
         DADRA = AF*AF*EXPE/(-BE*PHI3)*((VCA-EC)/RHO)
         DADRB = AF*AF*EXPE/(-BE*PHI3)*((VCB-EC)/RHO)
         DSDA  = -BE/GA * AF * T**6 * (2.+Y) / (1.+Y+Y*Y)**2
         DSDT  = 2.*BE/GA * T * (1.+2.*Y) / (1.+Y+Y*Y)**2
         DSDRA = DSDA*DADRA + DSDT*DTDRA
         DSDRB = DSDA*DADRB + DSDT*DTDRB
         DHDT  = GA*PHI3/S1*DSDT
         DHDRA = GA*PHI3/S1*DSDRA
         DHDRB = GA*PHI3/S1*DSDRB
       else
         DTDPHI= -T/PHI
         DPHIDE= 1./(3.*(1.d0+ETA)**OB3) - 1./(3.*(1.d0-ETA)**OB3)
         DEDRA = 2.*RHOB/(RHO*RHO)
         DEDRB = -2.*RHOA/(RHO*RHO)
         DPHIDA= DPHIDE*DEDRA
         DPHIDB= DPHIDE*DEDRB
         DTDRA = -T*(DPHIDA/PHI+7./(6.*RHO))
         DTDRB = -T*(DPHIDB/PHI+7./(6.*RHO))
         DADRA = AF*AF*EXPE/(-BE*PHI3)*(3.*EC/PHI*DPHIDA-(VCA-EC)/RHO)
         DADRB = AF*AF*EXPE/(-BE*PHI3)*(3.*EC/PHI*DPHIDB-(VCB-EC)/RHO)
         DSDA  = -BE/GA * AF * T**6 * (2.+Y) / (1.+Y+Y*Y)**2
         DSDT  = 2.*BE/GA * T * (1.+2.*Y) / (1.+Y+Y*Y)**2
         DSDRA = DSDA*DADRA + DSDT*DTDRA
         DSDRB = DSDA*DADRB + DSDT*DTDRB
         DHDT  = GA*PHI3/S1*DSDT
         DHDRA = 3.*H0/PHI*DPHIDA + GA*PHI3/S1*DSDRA
         DHDRB = 3.*H0/PHI*DPHIDB + GA*PHI3/S1*DSDRB
       endif

!..
      V1CA  = H0 + RHO*DHDRA
      V2CA  = RHO*DHDT*T/AA
      V1CB  = H0 + RHO*DHDRB
      V2CB  = RHO*DHDT*T/AA
      V2CAB = RHO*DHDT*T/AA
!     ==--------------------------------------------------------------==
      RETURN
      END  SUBROUTINE LSD_PBEC
!     ==================================================================

      end module local_spin_density
