!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... These are PHONON-specific modules (Conjugate Gradient version)
!
MODULE phunits
  !
  SAVE
  !
  CHARACTER(len=80) :: &
      fildyn
  CHARACTER(len=75) :: &
      title_ph
  INTEGER :: &
      iuwfc,     &!
      iubar,     &!
      iudwf,     &!
      iuscf,     &!
      iuvkb,     &!
      lrwfc,     &!
      lrbar,     &!
      lrdwf,     &!
      lrscf       !
  !
END MODULE phunits
!
!
MODULE flags
  !
  SAVE
  !
  LOGICAL :: &
      trans,     &!
      epsil,     &!
      raman,     &!
      equil,     &!
      nlcc_any,  &!
      asr         !
  !    
END MODULE flags
!
!
MODULE dielectric
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  REAL(kind=DP) :: &
      epsilon0(3,3)  
  REAL(kind=DP), ALLOCATABLE :: &
      zstar(:,:,:)
  !    
END MODULE dielectric
!
!
MODULE modes1
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  INTEGER :: &
      nmodes
  REAL(kind=DP), ALLOCATABLE :: & 
      dyn(:,:),   &!
      u(:,:)       !
  !
END MODULE modes1
!
!
MODULE cgconv
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  INTEGER :: &
      niter_ph
  REAL(kind=DP) :: &
      tr2_ph
  !    
END MODULE cgconv
!
!
MODULE AA
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       aux2(:),   &!
       aux3(:)     !
  REAL(KIND=DP), ALLOCATABLE, TARGET :: &
       auxr(:)     !
  !     
END MODULE AA
!
!
MODULE dmu
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  REAL(KIND=DP), ALLOCATABLE:: &
       dmuxc(:),        &!  d V_xc / d rho
       grho(:,:,:),     &!  gradient of the unperturbed density
       dvxc_rr(:,:,:),  &!
       dvxc_sr(:,:,:),  &!  derivatives of the E_xc functional w.r.t.
       dvxc_ss(:,:,:),  &!  r=rho and s=|grad(rho)|
       dvxc_s (:,:,:)
  !    
END MODULE dmu
!
!
MODULE phon
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE:: &
       dvpsi(:,:),      &!
       dpsi(:,:)         !
  !     
END MODULE phon
!
!
MODULE symmetry
  !
  SAVE
  !
  INTEGER :: &
       n_diff_sites,      &!
       nasr
  INTEGER, ALLOCATABLE :: &
       equiv_atoms(:,:),  &!
       n_equiv_atoms(:)    !
  INTEGER, ALLOCATABLE :: &
       has_equivalent(:)   !
  !
END MODULE symmetry
!
!
MODULE diffs
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  INTEGER :: &
       nderiv,     &!
       first,      &!
       last         !
  REAL(kind=DP) :: &
       deltatau     !
  !     
END MODULE diffs
!
!
MODULE cgcom
  USE cgconv
  USE phunits
  USE flags
  USE modes1
  USE AA
  USE phon
  USE diffs
  USE dmu
  USE symmetry
  USE dielectric
END MODULE cgcom
