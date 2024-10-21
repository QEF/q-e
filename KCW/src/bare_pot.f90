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
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvecs,                ONLY : ngms
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : tpiba2, omega
  USE control_kcw,          ONLY : spin_component, kcw_iverbosity, x_q, nrho
  USE gvect,                ONLY : g
  USE qpoint,               ONLY : xq
  USE constants,            ONLY : e2, fpi
  USE eqv,                  ONLY : dmuxc
  USE control_lr,           ONLY : lrpa
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE io_global,            ONLY : stdout
  !USE exx_base,             ONLY : g2_convolution
  USE coulomb,             ONLY : g2_convolution
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
  ! GC LR suff
  USE xc_lib,            ONLY : xclib_dft_is
  USE scf,               ONLY : rho, rho_core
  USE uspp,              ONLY : nlcc_any
  USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE dv_of_drho_lr,     ONLY : dv_of_drho_xc


  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT (OUT) :: rhog(ngms,nrho)
  ! ... periodic part of wannier density in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg_(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (OUT) :: vh_rhog(ngms)
  ! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  ! 
  COMPLEX(DP), INTENT (OUT) :: delta_vr(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vr_(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (IN)  :: rhor(dffts%nnr,nrho)
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
  INTEGER                   :: ig, is, ir, iss, ip
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
  COMPLEX(DP), ALLOCATABLE :: rhor_(:,:)
  !
  !! The periodic part of the orbital density in g space  
  DO ip=1,nrho !<---CONSIDER TO SUBSTITUTE WITH nspin_mag
      aux(:) = rhor(:,ip)/omega
      CALL fwfft ('Rho', aux, dffts)  
      rhog(:,ip) = aux(dffts%nl(:))
  END DO
  !
  delta_vr = ZERO
  delta_vr_ = ZERO
  aux      = ZERO
  !
  ! .. First the xc contribution
  !
  IF (.NOT. lrpa) THEN 
    ALLOCATE ( rhor_(dffts%nnr,nspin_mag) )
    rhor_ = CMPLX(0.D0,0.D0,kind=DP)
    IF (nspin_mag == 4) THEN
      rhor_=rhor/omega
    ELSE
      rhor_(:,spin_component) = rhor(:,1)/omega
    ENDIF
    CALL dv_of_drho_xc(delta_vr, rhor_)
    DEALLOCATE (rhor_)
  ENDIF 
  !
  delta_vr_ = delta_vr 
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
    vh_rhog_g0eq0(ig) =  e2 * fpi * rhog(ig,1) / (tpiba2 * qg2)
    IF (qg2 .lt. 1e-8) vh_rhog_g0eq0(ig) = (0.D0, 0.D0)
    ! ... set to zero the q+g=0 component
    !
    vh_rhog(ig) =  rhog(ig,1) * cmplx(fac(ig), 0.d0)
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
  aux(dffts%nl(:))  = vh_rhog(:)                    ! G-space components of the Hartree potential
  aux_(dffts%nl(:)) = vh_rhog_g0eq0(:)
  CALL invfft ('Rho', aux, dffts)
  CALL invfft ('Rho', aux_, dffts)
  !
  IF (nspin_mag==4 .and. domag) THEN    ! Perturbing potential due to Hartree
    delta_vr(:,1) = delta_vr(:,1)   + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
  ELSEIF (nspin_mag==4 .and. .NOT. domag) THEN
    delta_vr(:,1) = delta_vr(:,1) + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
  ELSE
    DO is = 1, nspin_mag
      delta_vr(:,is)  = delta_vr(:,is)   + aux(:)
      delta_vr_(:,is) = delta_vr_(:,is) + aux_(:)
    END DO
  END IF
  !
  DO is = 1, nspin_mag
    !
    aux(:) = delta_vr(:,is)
    aux_(:) = delta_vr_(:,is) 
    !
    CALL fwfft ('Rho', aux, dffts)
    CALL fwfft ('Rho', aux_, dffts)
    !
    delta_vg(:,is)  = aux(dffts%nl(:))
    delta_vg_(:,is) = aux_(dffts%nl(:))
    !
  ENDDO
  !
  !
END subroutine bare_pot
