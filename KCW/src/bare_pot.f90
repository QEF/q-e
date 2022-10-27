! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
  !---------------------------------------------------------------------
  !
  !! This routine compute the Hxc potential due to the (peridoc part of the) wannier
  !! charge density. V^{0n}_{Hxc}(r) = \int dr' f_Hxc(r,r') w^{0n}(r')
  !  
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvecs,                ONLY : ngms
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : tpiba2, omega
  USE control_kcw,          ONLY : spin_component, kcw_iverbosity, x_q 
  USE gvect,                ONLY : g
  USE qpoint,               ONLY : xq
  USE constants,            ONLY : e2, fpi
  USE eqv,                  ONLY : dmuxc
  USE control_lr,           ONLY : lrpa
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE io_global,            ONLY : stdout
  !USE exx_base,             ONLY : g2_convolution
  USE coulomb,             ONLY : g2_convolution
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT (OUT) :: rhog(ngms)
  ! ... periodic part of wannier density in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg(ngms,nspin)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg_(ngms,nspin)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (OUT) :: vh_rhog(ngms)
  ! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  ! 
  COMPLEX(DP), INTENT (OUT) :: delta_vr(dffts%nnr,nspin)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vr_(dffts%nnr,nspin)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (IN)  :: rhor(dffts%nnr)
  ! ... periodic part of wannier density in r-space
  !
  INTEGER, INTENT (IN)      :: iq
  ! ... q-point index
  !
  COMPLEX(DP)               :: aux(dffts%nnr), aux_(dffts%nnr)
  ! ... auxiliary vectors 
  !
  COMPLEX(DP), ALLOCATABLE  :: vaux(:)
  ! ... auxiliary vector
  !
  INTEGER                   :: ig, is, ir
  ! ... counters 
  !
  REAL(DP)                  :: qg2, eh_corr
  ! ... |q+G|^2, g=0 correction to the hartree energy 
  !
  LOGICAL                   :: lgamma
  !
  REAL(DP)                  :: xkq(3), xk(3)
  ! ... coordinate of k and k+q 
  !
  REAL(DP)                  :: fac(ngms)
  ! ... Coulomb kernel possibly with the special treatment of the q+g+0 component 
  !
  COMPLEX(DP)               :: vh_rhog_g0eq0(ngms)
  ! ... The hartree potential with th q+g=0 component set to zero
  !
  aux(:) = rhor(:)/omega
  CALL fwfft ('Rho', aux, dffts)  ! NsC: Dense or smooth grid?? I think smooth is the right one. 
  rhog(:) = aux(dffts%nl(:))
  !! The periodic part of the orbital density in g space  
  !
  delta_vr = ZERO
  delta_vr_ = ZERO
  aux      = ZERO
  !
  xk(:) = 0.D0
  xkq(:) = -x_q(:,iq)
  ! ... auxiliary coordinate to pass to g2_convoltion 
  ! 
  ! ... The Hartree contribution first 
  !
  CALL g2_convolution(ngms, g, xk, xkq, fac)
  ! ... the hartree kernel (eventually within the 
  ! ... Gygi-Baldereschi sheme, see setup_coulomb) 
  !
  DO ig = 1, ngms
    !
    qg2 = SUM ( (g(:,ig)+x_q(:,iq))**2 )
    !
    vh_rhog_g0eq0(ig) =  e2 * fpi * rhog(ig) / (tpiba2 * qg2)
    IF (qg2 .lt. 1e-8) vh_rhog_g0eq0(ig) = (0.D0, 0.D0)
    ! ... set to zero the q+g=0 component
    !
    vh_rhog(ig) =  rhog(ig) * cmplx(fac(ig), 0.d0)
    ! ... the hartree potential possibly with the special treatment of the q+g=0 component  
    !
  ENDDO
  !
  ! ... eventually add MT corrections
  !
  lgamma = (xq(1)==0.D0.AND.xq(2)==0.D0.AND.xq(3)==0.D0)
  IF (do_comp_mt .AND. lgamma) then
     ! ... this make sense only for a GAMMA only calculation 
     ! ... (do_compt_mt = .false if more than one k-point (see kcw_readin.f90) 
     !
     IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(5x, " ADDING Martyna-Tuckerman correction" ,/)')
     !
     ALLOCATE( vaux( ngms ) )
     CALL wg_corr_h (omega, ngms, rhog, vaux, eh_corr)
     vh_rhog(1:ngms) = vh_rhog(1:ngms) +  vaux(1:ngms)
     vh_rhog_g0eq0(1:ngms) = vh_rhog_g0eq0(1:ngms) +  vaux(1:ngms)
     DEALLOCATE( vaux )
  ENDIF
  !
  ! ... Go to r-space 
  !
  aux=(0.d0,0.d0)
  aux_=(0.d0,0.d0)
  aux(dffts%nl(:)) = vh_rhog(:)
  aux_(dffts%nl(:)) = vh_rhog_g0eq0(:)
  CALL invfft ('Rho', aux, dffts)
  CALL invfft ('Rho', aux_, dffts)
  !
  DO is = 1, nspin
     !
     delta_vr(:,is) = aux(:)
     delta_vr_(:,is) = aux_(:)
     !
  END DO
  !
  ! .. Add the xc contribution
  !
  IF (.NOT. lrpa) THEN
     !
     DO is = 1, nspin
        !
        DO ir = 1, dffts%nnr
           !
           delta_vr(ir,is) = delta_vr(ir,is) + dmuxc(ir,is,spin_component) * rhor(ir)/omega
           delta_vr_(ir,is) = delta_vr_(ir,is) + dmuxc(ir,is,spin_component) * rhor(ir)/omega
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  ! ... Back to g-space
  !
  DO is = 1, nspin
     !
     aux(:) = delta_vr(:,is)
     aux_(:) = delta_vr_(:,is) 
     !
     CALL fwfft ('Rho', aux, dffts)
     CALL fwfft ('Rho', aux_, dffts)
     !
     delta_vg(:,is) = aux(dffts%nl(:))
     delta_vg_(:,is) = aux_(dffts%nl(:))
     !
  ENDDO
  !
  !
END subroutine bare_pot
