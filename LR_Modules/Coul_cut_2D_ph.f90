!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! 
! This module contains the variables and subroutines necessary to cutoff
! the quantities used in the PHonon code. Details of the implementation can be found in:
!
! Sohier, T., Calandra, M., & Mauri, F. (2017), 
! "Density functional perturbation theory for gated two-dimensional heterostructures: 
! Theoretical developments and application to flexural phonons in graphene." 
! Physical Review B, 96(7), 75448. https://doi.org/10.1103/PhysRevB.96.075448
!
!----------------------------------------------------------------------------
MODULE Coul_cut_2D_ph
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the variables and subroutines needed for the 
  ! ... 2D Coulomb cutoff 
  !
  USE kinds, ONLY :  DP
  USE constants, ONLY : tpi, pi
  SAVE
  !
  REAL(DP), ALLOCATABLE :: cutoff_2D_qg(:)
  ! cutoff factor at q+G
  COMPLEX(DP), ALLOCATABLE :: lr_Vlocq(:,:)
  ! cutoff version of the long-range pat of the local ionic potential
  ! 
  !
  PUBLIC :: cutoff_fact_qg, cutoff_lr_Vlocq, cutoff_localq
  
CONTAINS
!
!----------------------------------------------------------------------
subroutine cutoff_fact_qg ()
  !----------------------------------------------------------------------
  !
  ! ...This routine calculates the cutoff factor at q+G and stores it in 
  !  a vector called cutoff_2D(:), to be re-used in various routines.
  !  see Eq. (24) of PRB 96, 075448
  !
  USE kinds
  USE gvect,     ONLY : g, ngm
  USE cell_base, ONLY : alat, celldm, at
  USE qpoint,    ONLY : xq 
  USE io_global, ONLY : stdout
  implicit none
  !
  !
  integer :: ng, i
  real(DP) :: qGzlz, qGplz, lz
  !
  IF (.not.ALLOCATED(cutoff_2D_qg)) ALLOCATE( cutoff_2D_qg (ngm) )  
  !  Check that material is in the x-y plane
  do i=1,2
     IF (abs(at(3,i))>1d-8) WRITE(stdout, *) "2D CODE WILL NOT WORK, 2D MATERIAL NOT IN X-Y PLANE!!"
  enddo
  ! redefine cutoff distance
  lz=0.5d0*at(3,3)*alat
  do ng = 1, ngm
     qGplz=SQRT( ( xq(1)+ g (1, ng) )**2 + &
         ( xq(2)+ g (2, ng) )**2 )*tpi*lz/alat
     qGzlz= ( xq(3)+ g (3, ng))*tpi*lz/alat
     ! sufficient if xq3=0
     cutoff_2D_qg(ng)= 1.0d0- exp(-qGplz)*cos(qGzlz)
     ! otherwise one needs to use. But I don see why would want 
     ! out-of-plane q for 2D material.
     !cutoff_2D_qg(ng)= 1.0d0+ exp(-qGplz)* &
     !        (qGzlz/qGplz*sin(qGzlz) -cos(qGzlz))
  enddo
  return
end subroutine cutoff_fact_qg
!
!----------------------------------------------------------------------
subroutine cutoff_lr_Vlocq ( )
  !----------------------------------------------------------------------
  !
  !  This routine calculates the long-range part of vloc(q+G) for 2D calculations.
  !  see Eq. (69) of PRB 96, 075448
  !
  USE kinds
  USE constants,   ONLY : fpi, e2, eps8
  USE gvect,      ONLY :  g, ngm
  !USE gvecs,      ONLY : ngms
  USE ions_base,  ONLY : zv, nsp
  USE uspp_param, ONLY : upf
  USE cell_base,  ONLY : omega, tpiba2
  USE qpoint,    ONLY : xq
  implicit none
  ! local variable
  integer :: ig, nt 
  REAL(DP)::fac, g2a
  ! counter over G vectors
  IF (.not.ALLOCATED(lr_Vlocq)) ALLOCATE( lr_Vlocq (ngm,nsp) )
  lr_Vlocq (:,:)=0.0d0
  DO nt = 1, nsp
     fac= upf(nt)%zp * e2 / tpiba2
     DO ig = 1, ngm
        g2a = (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) **2 + &
                (xq (3) + g (3, ig) ) **2  
        if (g2a<eps8) then 
           lr_Vlocq (ig, nt)=0.0d0
        else
           lr_Vlocq (ig, nt)= - fpi / omega* fac * cutoff_2D_qg(ig)* &
                    exp ( - g2a * tpiba2 * 0.25d0) / g2a
        endif 
     END DO
  END DO
  return
end subroutine cutoff_lr_Vlocq
!
!----------------------------------------------------------------------
subroutine cutoff_localq (dvlocin, fact,  u1, u2, u3, gu0, nt, na)
  !--------------------------------------------------------------------
  !
  !  This subroutine is called to re-add the long-range part of the 
  !  derivative of teh local part of the ionic potential, using lr_Vlocq 
  !  computed in above routine.
  !  see Eq. (68) of PRB 96, 075448
  !
  USE kinds
  USE fft_base,  ONLY : dffts
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g
  USE gvecs,   ONLY : ngms, nls
  implicit none
  !
  complex(DP), INTENT(INOUT) :: dvlocin (dffts%nnr)
  COMPLEX(DP), INTENT(IN) :: fact,  u1, u2, u3, gu0
  INTEGER, INTENT(IN)     :: nt, na
  ! input/output: derivative of vloc
  ! input : some local vatriables from the routines calling this one
  !
  ! local variables
  integer :: ig
  complex(DP) :: gtau, gu
  !
  do ig = 1, ngms
     gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
              eigts3 (mill(3,ig), na)
     gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
     dvlocin (nls (ig) ) = dvlocin (nls (ig) ) + lr_Vlocq (ig, nt) &
              * gu * fact * gtau
  enddo  
  return
end subroutine cutoff_localq
!
!----------------------------------------------------------------------
subroutine cutoff_dv_of_drho (dvaux, is, dvscf)
  !----------------------------------------------------------------------
  !
  ! This subroutine cuts off the Hartree potential generated by the 
  ! induced change in electronic density 
  ! see Eq. (70) of PRB 96, 075448
  !
  USE kinds
  USE constants,        ONLY : fpi, e2
  USE cell_base,        ONLY : tpiba2
  USE fft_base,         ONLY: dfftp
  USE noncollin_module, ONLY : nspin_mag
  USE gvect,            ONLY : g, ngm, nl
  USE qpoint,           ONLY : xq 
  implicit none
  !
  COMPLEX(DP), INTENT(INOUT) :: dvaux( dfftp%nnr,  nspin_mag)
  COMPLEX(DP), INTENT(IN) :: dvscf (dfftp%nnr, nspin_mag)
  INTEGER, INTENT(IN)     :: is
  ! input / output : dv
  ! input : drho
  ! input : spin counter
  !
  !local variables
  INTEGER :: ig
  REAL(DP) :: qg2
  !
  do ig = 1, ngm
     qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
     if (qg2 > 1.d-8) then
        dvaux(nl(ig),is) = dvaux(nl(ig),is) + cutoff_2D_qg(ig)*&
                           e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
     endif
  enddo
  return
end subroutine cutoff_dv_of_drho
!
!----------------------------------------------------------------------
subroutine cutoff_dynmat0 (dynwrk, rhog)
  !----------------------------------------------------------------------
  !
  ! This subroutine cuts off d^2Vloc/dudu * rho^0 
  ! See Eq. (66) in PRB 96, 075448
  !
  USE kinds
  USE constants,   ONLY : tpi, eps8
  USE cell_base,   ONLY : omega, tpiba2
  USE fft_base,    ONLY: dfftp
  USE gvect,       ONLY : g, ngm, nl, gg
  USE Coul_cut_2D, ONLY : lr_Vloc
  USE ions_base,   ONLY : nat, ityp, ntyp => nsp, tau
  implicit none
  !
  COMPLEX(DP), INTENT(INOUT) :: dynwrk (3 * nat, 3 * nat)
  COMPLEX(DP), INTENT(IN) :: rhog  ( dfftp%nnr)
  !
  integer :: ng, na, icart, jcart, na_icart, na_jcart
  REAL(DP) :: gtau, fac
  !
  DO na = 1, nat
     DO icart = 1, 3
        na_icart = 3 * (na - 1) + icart
        DO jcart = 1, 3
           na_jcart = 3 * (na - 1) + jcart
           DO ng = 1, ngm
              gtau = tpi * (g (1, ng) * tau (1, na) + &
                            g (2, ng) * tau (2, na) + &
                            g (3, ng) * tau (3, na) )
              fac = omega * lr_Vloc ( ng , ityp (na) ) * tpiba2 * &
                   ( DBLE (rhog (nl (ng) ) ) * COS (gtau) - &
                    AIMAG (rhog (nl (ng) ) ) * SIN (gtau) )
                  dynwrk (na_icart, na_jcart) = dynwrk (na_icart, na_jcart) - &
                   fac * g (icart, ng) * g (jcart, ng)
           ENDDO
        ENDDO
     ENDDO
  ENDDO  
  return
end subroutine cutoff_dynmat0
END MODULE Coul_cut_2D_ph
