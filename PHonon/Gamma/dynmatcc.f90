!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE dynmatcc(dyncc)
  !--------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp, tau
  USE atom,       ONLY : rgrid
  USE constants,  ONLY : tpi
  USE cell_base,  ONLY : omega, tpiba2
  USE ener,       ONLY : etxc, vtxc
  USE uspp_param, ONLY : upf
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect,      ONLY : nl, ngm, igtongl, ngl, g, gg, gl
  USE scf,        ONLY : rho, rho_core, rhog_core
  USE wavefunctions_module,  ONLY: psic
  USE cgcom
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  IMPLICIT NONE
  real(DP):: dyncc(3*nat,nmodes)
  !
  INTEGER:: i,j,na,nb,nta,ntb,ir,ig,nt, nu_i,nu_j,mu_i,mu_j
  COMPLEX(DP), POINTER:: vxc(:), work1(:), gc(:,:)
  COMPLEX(DP) :: exc
  real(DP), ALLOCATABLE:: rhocg(:), dyncc1(:,:,:,:)
  real(DP) :: exg
  LOGICAL :: nlcc(ntyp)
  !
  !
  dyncc(:,:) = 0.d0
  !
  IF ( any( upf(1:ntyp)%nlcc ) ) GOTO 10
  RETURN
10 CONTINUE
  !
  work1 => psic
  vxc   => aux2
  ALLOCATE  ( dyncc1( 3,nat,3,nat))
  ALLOCATE  ( gc    ( dfftp%nnr, 3))
  ALLOCATE  ( rhocg( ngl))
  !
  CALL v_xc  (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  CALL fwfft ( 'Dense', vxc, dfftp )
  !
  dyncc1(:,:,:,:) = 0.d0
  ! temporary
  nlcc(1:ntyp) =  upf(1:ntyp)%nlcc
  DO na=1,nat
     nta=ityp(na)
     IF ( upf(nta)%nlcc ) THEN
        CALL drhoc (ngl, gl, omega, tpiba2, rgrid(nta)%mesh, rgrid(nta)%r, &
                    rgrid(nta)%rab, upf(nta)%rho_atc, rhocg)
        DO ig=1,ngm
           exg = tpi* ( g(1,ig)*tau(1,na) + &
                        g(2,ig)*tau(2,na) + &
                        g(3,ig)*tau(3,na) )
           exc = cmplx(cos(exg),-sin(exg),kind=DP)*tpiba2
           work1(ig)= rhocg(igtongl(ig))* exc * conjg(vxc(nl(ig)))
           gc(ig,1) = g(1,ig) * exc * (0.0d0,-1.0d0)
           gc(ig,2) = g(2,ig) * exc * (0.0d0,-1.0d0)
           gc(ig,3) = g(3,ig) * exc * (0.0d0,-1.0d0)
        ENDDO
        DO i=1,3
           DO j=1,3
              DO ig=1,ngm
                 dyncc1(i,na,j,na) = dyncc1(i,na,j,na) -  &
                      dble(work1(ig)) * g(i,ig) * g(j,ig)
              ENDDO
           ENDDO
        ENDDO
        DO i=1,3
           CALL dvb_cc  (nlcc, nt, ngm, dfftp%nnr, &
                nl,igtongl,rhocg,dmuxc,gc(1,i),aux3,gc(1,i))
        ENDDO
        DO nb=1,nat
           ntb=ityp(nb)
           IF ( upf(ntb)%nlcc ) THEN
              CALL drhoc (ngl, gl, omega, tpiba2, rgrid(ntb)%mesh, &
                          rgrid(ntb)%r, rgrid(ntb)%rab, upf(ntb)%rho_atc,&
                          rhocg)
              DO ig=1,ngm
                 exg = tpi* ( g(1,ig)*tau(1,nb) + &
                              g(2,ig)*tau(2,nb) + &
                              g(3,ig)*tau(3,nb) )
                 exc = -cmplx(sin(exg),cos(exg),kind=DP)
                 work1(ig) = exc * rhocg(igtongl(ig))
              ENDDO
              DO i=1,3
                 DO j=1,3
                    DO ig=1,ngm
                       dyncc1(i,na,j,nb) = dyncc1(i,na,j,nb) +      &
                            dble( work1(ig)*conjg(gc(ig,i)))*g(j,ig)
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  DEALLOCATE(rhocg)
  DEALLOCATE(gc)
#if defined(__MPI)
  CALL mp_sum( dyncc1, intra_pool_comm )
#endif
  CALL dscal(3*nat*3*nat,-omega,dyncc1,1)
  !
  ! dyncc1 contains the entire dynamical matrix (core-correction part)
  ! in cartesian coordinates: transform to generic modes
  !
  DO nu_i=1,nmodes
     IF ( has_equivalent((nu_i-1)/3+1)==0 ) THEN
        DO nu_j=1,nmodes
           DO mu_i=1,3*nat
              na=(mu_i-1)/3+1
              i = mu_i-3*(na-1)
              DO mu_j=1,3*nat
                 nb=(mu_j-1)/3+1
                 j = mu_j-3*(nb-1)
                 dyncc(nu_i,nu_j) = dyncc(nu_i,nu_j) +              &
                      dyncc1(i,na,j,nb)*u(mu_i,nu_i)*u(mu_j,nu_j)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  DEALLOCATE(dyncc1)
  !
  RETURN
END SUBROUTINE dynmatcc
