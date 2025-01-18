!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine sym_dmag (dmagtosym)
  !---------------------------------------------------------------------
  !! Symmetrize the change of the magnetization density belonging to
  !! an irreducible representation.  
  !! The routine is generalized to include also the symmetry operations
  !! that require the time-reversal operator (meaning that TS is a 
  !! symmetry of the crystal).  
  !! For a more complete explanation, please see: 
  !! Phys. Rev. B 100, 045115 (2019).
  !
  USE kinds, only : DP
  USE constants, ONLY: tpi
  USE fft_base, ONLY: dfftp
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, ft, t_rev, sname, invs
  USE noncollin_module, ONLY: nspin_mag
  USE lr_symm_base, ONLY : minus_q, nsymq, gi, upert, lr_npert

  implicit none

  complex(DP) :: dmagtosym(dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag, lr_npert)
  !! the magnetization to symmetrize (only 2:4 components)
  !
  ! ... local variables
  !
  integer :: ftau(3,nsymq), s_scaled(3,3,nsymq)
  integer :: is, ri, rj, rk, i, j, k, ipert, jpert, ipol, isym, kpol
  !  counter on spin polarizations
  !
  !  the rotated points
  !
  !
  !  counter on mesh points
  !
  ! counter on perturbations
  ! counter on perturbations
  ! counter on polarizations
  ! counter on symmetries
  ! the rotation

  real(DP) :: g1 (nsymq), g2 (nsymq), g3 (nsymq), in1, in2, in3
  ! used to construct the phases
  ! auxiliary variables

  complex(DP), allocatable :: dmagsym (:,:,:,:,:), dmags(:,:)
  ! the symmetrized potential
  complex(DP) ::  term (3, nsymq), phase (nsymq), mag(3), magrot(3)
  ! auxiliary space
  ! the multiplication factor
  ! the phase factor

  ! For noncollinear magnetism, time_reversal is false, and thus minus_q is always false.
  ! This is because all symmetries are combined spatial-and-time-reversal symmetries, with
  ! the time-reversal component stored in t_rev(isym). So, symmetry q -> -q is not treated
  ! separately for noncollinear magnetism. (See phq_setup.f90 and set_small_group_of_q.f90.)
  !
  ! This routine is called only under noncollinear magnetism. Thus, minus_q must be false.
  ! If not, there is a problem in the symmetry setup.
  IF (minus_q) THEN
     CALL errore("sym_dmag", "minus_q must be false for noncollinear magnetism."&
                            &"sym_dmag must not be called.", 1)
  ENDIF
  !
  if (nsymq == 1.and. (.not.minus_q) ) return
  !
  call start_clock ('sym_dmag')
  !
  allocate (dmagsym(  dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , 3, lr_npert))
  allocate (dmags( 3, lr_npert))
  !
  in1 = tpi / DBLE (dfftp%nr1)
  in2 = tpi / DBLE (dfftp%nr2)
  in3 = tpi / DBLE (dfftp%nr3)
  !
  CALL scale_sym_ops( nsymq, s, ft, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
       s_scaled, ftau )
  !
  ! Here we symmetrize with respect to the small group of q
  !
  do isym = 1, nsymq
     g1 (isym) = 0.d0
     g2 (isym) = 0.d0
     g3 (isym) = 0.d0
     do ipol = 1, 3
        g1 (isym) = g1 (isym) + gi (ipol, isym) * at (ipol, 1)
        g2 (isym) = g2 (isym) + gi (ipol, isym) * at (ipol, 2)
        g3 (isym) = g3 (isym) + gi (ipol, isym) * at (ipol, 3)
     enddo
     g1 (isym) = NINT(g1(isym))*in1
     g2 (isym) = NINT(g2(isym))*in2
     g3 (isym) = NINT(g3(isym))*in3
     term (1, isym) = CMPLX(cos (g1 (isym) ), sin (g1 (isym) ) ,kind=DP)
     term (2, isym) = CMPLX(cos (g2 (isym) ), sin (g2 (isym) ) ,kind=DP)
     term (3, isym) = CMPLX(cos (g3 (isym) ), sin (g3 (isym) ) ,kind=DP)
  enddo

  dmagsym(:,:,:,:,:) = (0.d0, 0.d0)
  do isym = 1, nsymq
     phase (isym) = (1.d0, 0.d0)
  enddo
  do k = 1, dfftp%nr3
     do j = 1, dfftp%nr2
        do i = 1, dfftp%nr1
           do isym = 1, nsymq
              ! rotate_grid_point finds the rotated of i,j,k with S^-1
              CALL rotate_grid_point(s_scaled(1,1,isym), ftau(1,isym), &
                   i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
              dmags=(0.d0,0.d0)
              do ipert = 1, lr_npert
                 do jpert = 1, lr_npert
                    do is=2,4
                       dmags(is-1,ipert)=dmags(is-1,ipert) + &
                            upert(jpert, ipert, isym) * &
                            dmagtosym (ri, rj, rk, is, jpert) * phase (isym)
                    enddo
                 enddo
                 do kpol = 1, 3
                    mag(kpol)=bg(1,kpol)*dmags(1,ipert) + &
                              bg(2,kpol)*dmags(2,ipert) + &
                              bg(3,kpol)*dmags(3,ipert)
                 enddo
                 ! rotate the magnetic moment
                 do kpol = 1, 3
                    magrot(kpol) = s(1,kpol,invs(isym))*mag(1) + &
                                   s(2,kpol,invs(isym))*mag(2) + &
                                   s(3,kpol,invs(isym))*mag(3)
                 enddo
                 if (sname(isym)(1:3)=='inv') magrot=-magrot
                 if(t_rev(isym).eq.1) magrot=-magrot
                 ! go back to cartesian coordinates
                 do kpol = 1, 3
                    mag(kpol)=at(kpol,1)*magrot(1) + &
                              at(kpol,2)*magrot(2) + &
                              at(kpol,3)*magrot(3)
                 enddo
                 if (t_rev(isym) == 1) then 
                    mag(:) = conjg(mag(:))
                 end if
                 dmagsym(i,j,k,1,ipert)=dmagsym(i,j,k,1,ipert)+mag(1)
                 dmagsym(i,j,k,2,ipert)=dmagsym(i,j,k,2,ipert)+mag(2)
                 dmagsym(i,j,k,3,ipert)=dmagsym(i,j,k,3,ipert)+mag(3)
              enddo
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

  do is=2,4
     do ipert = 1, lr_npert
        dmagtosym(:,:,:,is,ipert) = dmagsym(:,:,:,is-1,ipert) / DBLE (nsymq)
     enddo
  enddo

  deallocate (dmags)
  deallocate (dmagsym)

  call stop_clock ('sym_dmag')
  return
end subroutine sym_dmag
