!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddens_pulay (drhop, dbecsum, mode0, npe)
  !----------------------------------------------------------------------
  !! Non-self-consistent, Pulay-like ultrasoft-augmentation contribution
  !! to drho for atomic-displacement (phonon) perturbations. Computes
  !! Eq. B31 of Ref [1]: the term that arises because the projectors
  !! |beta> move with the ions, built from the ground-state becsum and
  !! alphasum together with the mode displacement u(:,mode).
  !!
  !! Distinct from LR_Modules/lr_addusddens.f90, which adds the
  !! self-consistent augmentation contribution computed from dbecsum
  !! (the dbecsum built from the Sternheimer dpsi). For phonons both
  !! contributions are needed; for electric-field / EELS / Hubbard
  !! perturbations only lr_addusddens applies, because the projectors
  !! do not move.
  !!
  !! dbecsum and drhop entering here already contain the
  !! orthogonalisation contribution to the change of the wavefunctions;
  !! the terms with alphasum and becsum are added on top.
  !! The contribution of the change of the Fermi energy is not
  !! calculated here but added later by ef_shift.
  !! [1] PRB 64, 235118 (2001).
  !
  !
  USE kinds, only : DP
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  USE gvect,  ONLY : ngm, g, eigts1, eigts2, eigts3, mill
  USE uspp,     ONLY : okvan, becsum
  USE cell_base, ONLY : tpiba
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE wavefunctions,  ONLY: psic
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE paw_variables, ONLY : okpaw
  USE modes,     ONLY : u
  USE phus,    ONLY : becsumort, alphasum
  USE noncollin_module, ONLY : nspin_mag
  USE qpoint,     ONLY : xq, eigqts

  implicit none
  !
  !   the dummy variables
  !

  integer :: npe
  !! input: the number of perturbations
  complex(DP) :: drhop (dfftp%nnr, nspin_mag, npe)
  !! inp/out: change of the charge density
  complex(DP) :: dbecsum (nhm*(nhm+1)/2, nat, nspin_mag, npe)
  !! input: sum over kv of bec
  integer ::  mode0
  !! input:the mode of the representation
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, mu, mode, ipert, is, ijh
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! pointer on modes
  ! pointer on the mode
  ! counter on perturbations
  ! counter on spin
  ! counter on combined beta functions

  real(DP), allocatable  :: qmod (:), qpg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP) :: fact, zsum, bb, alpha, u1, u2, u3, fact_alpha_bb
  ! auxiliary variables
  complex(DP), allocatable ::  sk (:), qgm(:), aux (:,:,:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  
  ! Arrays to store pre-computed mode-independent quantities  
  complex(DP), allocatable :: alpha_array(:)
  ! pre-computed alpha values for vectorization  
  complex(DP), allocatable :: qgm_all(:,:)
  integer, allocatable :: ijh_map(:,:)
  integer :: n_entries, ientry
  ! qgm values for all (nt,ih,jh) combinations  
  ! mapping: (nt, ijh) for each entry
  ! total number of (nt,ih,jh) entries
  ! counter for entries
  !
  IF (.not.okvan) return
  !
  CALL start_clock ('addusddens_pulay')
  !
  ! Pre-compute mode-independent quantities
  ! First count the number of (nt,ih,jh) entries
  n_entries = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp) then
        do ih = 1, nh(nt)
           do jh = ih, nh(nt)
              n_entries = n_entries + 1
           enddo
        enddo
     endif
  enddo
  !
  ! Allocate arrays for pre-computed quantities
  allocate (aux(  ngm , nspin_mag , npe))
  allocate (alpha_array(ngm))
  allocate (qgm_all(ngm, n_entries))
  allocate (ijh_map(2, n_entries))  ! stores (nt, ijh)
  allocate (sk (  ngm))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm(  ngm))
  allocate (qmod( ngm))
  allocate (qpg( 3  , ngm))
  !      WRITE( stdout,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !  These computations are independent of modes, so compute once for all
  !
  call setqmod (ngm, xq, g, qmod, qpg)
  call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (qmod (ig) ) * tpiba
  enddo
  !
  ! Pre-compute qgm for all (nt,ih,jh) combinations
  ientry = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp) then
        ijh = 0
        do ih = 1, nh(nt)
           do jh = ih, nh(nt)
              ijh = ijh + 1
              ientry = ientry + 1
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              !
              ! Store mapping information and qgm
              ijh_map(1, ientry) = nt
              ijh_map(2, ientry) = ijh
              do ig = 1, ngm
                 qgm_all(ig, ientry) = qgm(ig)
              enddo
           enddo
        enddo
     endif
  enddo
  !
  fact = cmplx (0.d0, - tpiba, kind=DP)
  aux(:,:,:) = (0.d0, 0.d0)
  !
  ! Use pre-computed qgm and compute structure factors on-the-fly
  !
  do ientry = 1, n_entries
     nt = ijh_map(1, ientry)
     ijh = ijh_map(2, ientry)
     !
     ! Loop over all atoms of this type
     do na = 1, nat
        if (ityp(na) .eq. nt) then
           mu = 3 * (na - 1)
           !
           ! Calculate structure factor for this (ientry, na) combination
           do ig = 1, ngm
              sk(ig) = eigts1 (mill(1,ig), na) * &
                       eigts2 (mill(2,ig), na) * &
                       eigts3 (mill(3,ig), na) * &
                       eigqts (na) * qgm_all(ig, ientry)
           enddo
           !
           ! Process perturbations - preliminary sums approach
           !
           do ipert = 1, npe
              mode = mode0 + ipert
              u1 = u (mu + 1, mode)
              u2 = u (mu + 2, mode)
              u3 = u (mu + 3, mode)
              
              ! Check if displacement is significant for u-dependent terms
              if (abs(u1) + abs(u2) + abs(u3) >= 1.0E-12_DP) then
                 ! Pre-compute alpha array only if displacement is significant
                 alpha_array(1:ngm) = qpg(1,1:ngm)*u1 + qpg(2,1:ngm)*u2 + qpg(3,1:ngm)*u3
              endif
              
              do is = 1, nspin_mag
                 zsum = dbecsum (ijh, na, is, ipert)
                 !
                 bb = becsum (ijh, na, is)
                 
                 ! Add u-dependent terms only if displacement is significant
                 if (abs(u1) + abs(u2) + abs(u3) >= 1.0E-12_DP) then
                    zsum = zsum + &
                         ( alphasum (ijh, 1, na, is) * u1 &
                         + alphasum (ijh, 2, na, is) * u2 &
                         + alphasum (ijh, 3, na, is) * u3)
                 endif
                 !
                 call start_clock('compute_auxiliary')
                 ! Vectorized update of aux array - only if displacement is significant
                 if (abs(u1) + abs(u2) + abs(u3) >= 1.0E-12_DP) then
                    fact_alpha_bb = fact * bb
                    do ig = 1, ngm
                       aux(ig,is,ipert) = aux(ig,is,ipert) + fact_alpha_bb * alpha_array(ig) * sk(ig)
                    enddo
                 endif
                 call zaxpy (ngm, zsum, sk, 1, aux(1,is,ipert), 1)
                 IF (okpaw) becsumort(ijh,na,is,mode) = zsum
              enddo
           enddo
        endif
     enddo
  enddo
  !
  !     convert aux to real space
  !
  do ipert = 1, npe
     do is = 1, nspin_mag
        psic(:) = (0.d0, 0.d0)
        do ig = 1, ngm
           psic (dfftp%nl (ig) ) = aux (ig, is, ipert)
        enddo
        CALL invfft('Rho', psic, dfftp)
        call zaxpy(dfftp%nnr, (1.d0, 0.d0), psic, 1, drhop(1,is,ipert), 1)
     enddo
  enddo
  deallocate (qpg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (sk)
  deallocate (aux)
  deallocate (alpha_array)
  deallocate (qgm_all)
  deallocate (ijh_map)
  !
  call stop_clock ('addusddens_pulay')
  !
end subroutine addusddens_pulay
