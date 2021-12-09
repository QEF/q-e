!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine symdvscf (nper, irr, dvtosym)
  !---------------------------------------------------------------------
  !! Symmetrize the self-consistent potential of the perturbations
  !! belonging to an irreducible representation.  
  !! The routine is generalized to include, in the noncollinear 
  !! magnetic case, also the symmetry operations that require the 
  !! time-reversal operator (meaning that TS is a symmetry of the 
  !! crystal).  
  !
  USE kinds, only : DP
  USE constants, ONLY: tpi
  USE fft_base,  ONLY: dfftp
  USE cell_base, ONLY : at
  USE symm_base, ONLY : s, ft, t_rev
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag
  USE modes,   ONLY : t, tmq
  USE lr_symm_base, ONLY : minus_q, irotmq, nsymq, gi, gimq 
  !
  implicit none
  !
  integer :: nper
  !! the number of perturbations
  integer :: irr
  !! the representation under conside
  complex(DP) :: dvtosym(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x,nspin_mag,nper)
  !! the potential to be symmetrized
  !
  ! ... local variables
  !
  integer :: ftau(3,nsymq), s_scaled(3,3,nsymq)
  integer :: is, ri, rj, rk, i, j, k, ipert, jpert, ipol, isym, &
       irot
  !  counters
  real(DP) :: gf(3), n(3)
  !  temp variables
  complex(DP), allocatable :: dvsym (:,:,:,:), add_dvsym(:)
  ! the symmetrized potential
  complex(DP) ::  aux2, term (3, 48), phase (48)
  ! auxiliary space
  ! the multiplication factor
  ! the phase factor

  if (nsymq == 1.and. (.not.minus_q) ) return
  call start_clock ('symdvscf')

  allocate (dvsym(  dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , nper))
  allocate (add_dvsym(nper))
  !
  ! if necessary we symmetrize with respect to  S(irotmq)*q = -q + Gi
  !
  n(1) = tpi / DBLE (dfftp%nr1)
  n(2) = tpi / DBLE (dfftp%nr2)
  n(3) = tpi / DBLE (dfftp%nr3)

  CALL scale_sym_ops( nsymq, s, ft, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
       s_scaled, ftau )

  if (minus_q) then
     gf(:) =  gimq (1) * at (1, :) * n(:) + &
              gimq (2) * at (2, :) * n(:) + &
              gimq (3) * at (3, :) * n(:)
     term (:, 1) = CMPLX(cos (gf (:) ), sin (gf (:) ) ,kind=DP)
     do is = 1, nspin_lsda
        phase (1) = (1.d0, 0.d0)
        do k = 1, dfftp%nr3
           do j = 1, dfftp%nr2
              do i = 1, dfftp%nr1
                 CALL rotate_grid_point(s_scaled(1,1,irotmq), ftau(1,irotmq), &
                      i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                 do ipert = 1, nper
                    aux2 = (0.d0, 0.d0)
                    do jpert = 1, nper
                       aux2 = aux2 + tmq (jpert, ipert, irr) * &
                            dvtosym (ri, rj, rk, is, jpert) * phase (1)
                    enddo
                    dvsym (i, j, k, ipert) = (dvtosym (i, j, k, is, ipert) +&
                         CONJG(aux2) ) * 0.5d0
                 enddo
                 phase (1) = phase (1) * term (1, 1)
              enddo
              phase (1) = phase (1) * term (2, 1)
           enddo
           phase (1) = phase (1) * term (3, 1)
        enddo
        do ipert = 1, nper
           dvtosym(:, :, :, is, ipert) = dvsym (:, :, :, ipert)
        enddo
     enddo
  endif
  !
  ! Here we symmetrize with respect to the small group of q
  !
  do isym = 1, nsymq
     gf(:) =  gi (1,isym) * at (1, :) * n(:) + &
              gi (2,isym) * at (2, :) * n(:) + &
              gi (3,isym) * at (3, :) * n(:)
     term (:, isym) = CMPLX(cos (gf (:) ), sin (gf (:) ) ,kind=DP)
  enddo
  
  do is = 1, nspin_lsda
     dvsym(:,:,:,:) = (0.d0, 0.d0)
     do isym = 1, nsymq
        phase (isym) = (1.d0, 0.d0)
     enddo
     do k = 1, dfftp%nr3
        do j = 1, dfftp%nr2
           do i = 1, dfftp%nr1
              do isym = 1, nsymq
                 irot = isym
                 CALL rotate_grid_point(s_scaled(1,1,irot), ftau(1,irot), &
                   i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                 add_dvsym(:) = (0.d0, 0.d0)
                 do ipert = 1, nper
                    do jpert = 1, nper
                       add_dvsym(ipert) = add_dvsym(ipert) + t (jpert, ipert, irot, irr) * &
                                 dvtosym (ri, rj, rk, is, jpert) * phase (isym)
                       !dvsym (i, j, k, ipert) = dvsym (i, j, k, ipert) + &
                       !     t (jpert, ipert, irot, irr) * &
                       !     dvtosym (ri, rj, rk, is, jpert) * phase (isym)
                    enddo
                 enddo
                 if (t_rev(isym)==0) then 
                    dvsym (i, j, k, :) = dvsym (i, j, k, :) + add_dvsym(:)
                 else
                    dvsym (i, j, k, :) = dvsym (i, j, k, :) + conjg(add_dvsym(:))
                 end if
              enddo
              do isym = 1, nsymq
                 phase (isym) = phase (isym) * term (1, isym)
              enddo
           enddo
           do isym = 1, nsymq
              phase (isym) = phase (isym) * term (2, isym)
           enddo
        enddo
        do isym = 1, nsymq
           phase (isym) = phase (isym) * term (3, isym)
        enddo
     enddo

     do ipert = 1, nper
        dvtosym(:,:,:,is,ipert) = dvsym(:,:,:,ipert) / DBLE (nsymq)
     enddo

  enddo
  deallocate (dvsym)
  deallocate (add_dvsym)

  call stop_clock ('symdvscf')
  return
end subroutine symdvscf
