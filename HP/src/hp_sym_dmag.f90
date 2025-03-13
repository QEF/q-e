!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine hp_sym_dmag (dmagtosym)
    !---------------------------------------------------------------------
    ! symmetrize the change of the magnetization density.
    ! The routine is generalized to include also the 
    ! symmetry operations that require the time-reversal 
    ! operator (meaning that TS is a symmetry of the crystal).
    ! For a more complete explanation, please see: 
    ! Phys. Rev. B 100, 045115 (2019).
    !
    ! Inspired by PHonon/PH/sym_dmag.f90 
    !
    USE kinds, ONLY : DP
    USE constants, ONLY: tpi
    USE fft_base, ONLY: dfftp
    USE cell_base, ONLY : at, bg
    USE symm_base, ONLY : s, ft, t_rev, sname, invs
    USE noncollin_module, ONLY: nspin_mag
    USE ions_base,  ONLY : tau
    USE lr_symm_base, ONLY : minus_q, irotmq, nsymq, gi, gimq
    USE ldaU_hp, ONLY : nah_pert
    USE qpoint,           ONLY : xq
    USE symm_base,        ONLY : s
    USE cell_base,        ONLY : bg, at
    USE noncollin_module, ONLY : noncolin
  
    implicit none
  
    complex(DP) :: dmagtosym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag)
    ! the magnetization to symmetrize (only 2:4 components)
    real(DP) :: gi_t(3, 48), aq (3), raq (3), wrk (3)
    integer :: ftau(3,nsymq), s_scaled(3,3,nsymq)
    integer :: is, ri, rj, rk, i, j, k, ipol, isym, irot, kpol
    !  counter on spin polarizations
    !
    !  the rotated points
    !
    !
    !  counter on mesh points
    !
    ! counter on polarizations
    ! counter on symmetries
    ! the rotation
  
    real(DP) :: g1 (48), g2 (48), g3 (48), in1, in2, in3, ggf(48)
    ! used to construct the phases
    ! auxiliary variables
  
    complex(DP), allocatable :: dmagsym (:,:,:,:), dmags(:)
    ! the symmetrized potential
    complex(DP) ::  aux2(3), term (3, 48), phase (48), phase2 (48), mag(3), magrot(3)
    ! auxiliary space
    ! the multiplication factor
    ! the phase factor
  
    if (nsymq == 1.and. (.not.minus_q) ) return
    call start_clock ('hp_sym_dmag')
  
    allocate (dmagsym(  dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , 3))
    allocate (dmags( 3))
    !
    ! if necessary we symmetrize with respect to  S(irotmq)*q = -q + Gi
    !
    in1 = tpi / DBLE (dfftp%nr1)
    in2 = tpi / DBLE (dfftp%nr2)
    in3 = tpi / DBLE (dfftp%nr3)
  
    CALL scale_sym_ops( nsymq, s, ft, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
         s_scaled, ftau )
    
    if (minus_q) then
       g1 (1) = 0.d0
       g2 (1) = 0.d0
       g3 (1) = 0.d0
       !
       do ipol = 1, 3
          g1 (1) = g1 (1) + gimq (ipol) * in1 * at (ipol, 1)
          g2 (1) = g2 (1) + gimq (ipol) * in2 * at (ipol, 2)
          g3 (1) = g3 (1) + gimq (ipol) * in3 * at (ipol, 3)
       enddo
       !
       term (1, 1) = CMPLX(cos (g1 (1) ), sin (g1 (1) ) ,kind=DP)
       term (2, 1) = CMPLX(cos (g2 (1) ), sin (g2 (1) ) ,kind=DP)
       term (3, 1) = CMPLX(cos (g3 (1) ), sin (g3 (1) ) ,kind=DP)
       !
       phase (1) = (1.d0, 0.d0)
       !
       do k = 1, dfftp%nr3
         !
          do j = 1, dfftp%nr2
            !
             do i = 1, dfftp%nr1
                !
                CALL rotate_grid_point(s_scaled(1,1,irotmq), ftau(1,irotmq), &
                     i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                !
                aux2 = (0.d0, 0.d0)
                do is=2,4
                   aux2(is-1) = aux2(is-1) + dmagtosym (ri, rj, rk, is) * phase (1)
                enddo
                !
                do kpol = 1, 3
                   mag(kpol)=bg(1,kpol)*aux2(1) + bg(2,kpol)*aux2(2) + &
                             bg(3,kpol)*aux2(3)
                enddo
                !
                ! rotate the magnetic moment
                do kpol = 1, 3
                   magrot(kpol) = s(1,kpol,invs(irotmq))*mag(1) + &
                                  s(2,kpol,invs(irotmq))*mag(2) + &
                                  s(3,kpol,invs(irotmq))*mag(3)
                enddo
                !
                if (sname(irotmq)(1:3)=='inv') magrot=-magrot
                if(t_rev(irotmq).eq.1) magrot=-magrot
                !
                ! go back to cartesian coordinates
                do kpol = 1, 3
                   mag(kpol)=at(kpol,1)*magrot(1) + &
                             at(kpol,2)*magrot(2) + &
                             at(kpol,3)*magrot(3)
                   dmagsym(i,j,k,kpol)=(dmagtosym(i,j,k,kpol+1)+&
                                     CONJG(mag(kpol)) ) * 0.5d0
                enddo
                !
                phase (1) = phase (1) * term (1, 1)
                !
             enddo
             !
             phase (1) = phase (1) * term (2, 1)
             !
          enddo
          !
          phase (1) = phase (1) * term (3, 1)
          !
       enddo
       !
       do is = 2, 4
         dmagtosym(:, :, :, is) = dmagsym (:, :, :, is-1)
       end do
       !
    endif
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
    !
    ! LB: could it be written as in sym_dmag.f90 ? 
    IF (noncolin) THEN
      gi_t = 0.d0 
      aq   = xq
      call cryst_to_cart (1, aq, at, - 1)
      do isym = 1, nsymq
         raq = 0.d0
         do i = 1, 3
            do j = 1, 3
               raq (i) = raq (i) + DBLE (s (i, j, isym) ) * &
                  aq (j)
         enddo
         enddo
         do i = 1, 3
            IF (t_rev(isym)==1) wrk (i) = aq (i) - raq (i)
         enddo
         call cryst_to_cart (1, wrk, bg, 1)
         gi_t (:, isym) = wrk (:)
      enddo
    ENDIF
   !
   DO isym = 1, nsymq
      ggf (isym) = 0.d0
      IF (t_rev(isym) == 1) then
         do ipol = 1, 3
            ggf (isym) = ggf (isym) + (gi_t (ipol, isym) * tau(ipol,nah_pert)) * tpi
         enddo
      ELSE
         do ipol = 1, 3
            ggf (isym) = ggf (isym) + (gi (ipol, isym) * tau(ipol,nah_pert)) * tpi
         enddo
      ENDIF
      phase2 (isym) = CMPLX( cos( ggf(isym) ), -sin( ggf(isym) ), kind=DP)
    ENDDO
    !
    !
    dmagsym(:,:,:,:) = (0.d0, 0.d0)
    do isym = 1, nsymq
       phase (isym) = (1.d0, 0.d0)
    enddo
    !
    do k = 1, dfftp%nr3
      !
       do j = 1, dfftp%nr2
         !
          do i = 1, dfftp%nr1
            !
             do isym = 1, nsymq
                irot = isym
                CALL rotate_grid_point(s_scaled(1,1,irot), ftau(1,irot), &
                  i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                dmags=(0.d0,0.d0)
                !
                do is=2,4
                   dmags(is-1)= dmags(is-1) + &
                        dmagtosym (ri, rj, rk, is) * phase (isym) 
                enddo
                !
                do kpol = 1, 3
                   mag(kpol)=bg(1,kpol)*dmags(1) + &
                             bg(2,kpol)*dmags(2) + &
                             bg(3,kpol)*dmags(3)
                enddo
                !
                ! rotate the magnetic moment
                do kpol = 1, 3
                   magrot(kpol) = s(1,kpol,invs(irot))*mag(1) + &
                                  s(2,kpol,invs(irot))*mag(2) + &
                                  s(3,kpol,invs(irot))*mag(3)
                enddo
                !
                if (sname(irot)(1:3)=='inv') magrot=-magrot
                if(t_rev(irot).eq.1) magrot=-magrot
                !
                ! go back to cartesian coordinates
                do kpol = 1, 3
                   mag(kpol)=at(kpol,1)*magrot(1) + &
                             at(kpol,2)*magrot(2) + &
                             at(kpol,3)*magrot(3)
                enddo
                !
                if (t_rev(isym) == 1) then 
                   mag(:) = conjg(mag(:)) * phase2 (isym)
                end if
                !
                dmagsym(i,j,k,1)=dmagsym(i,j,k,1)+mag(1)
                dmagsym(i,j,k,2)=dmagsym(i,j,k,2)+mag(2)
                dmagsym(i,j,k,3)=dmagsym(i,j,k,3)+mag(3)
                !
             enddo
             !
             do isym = 1, nsymq
                phase (isym) = phase (isym) * term (1, isym)
             enddo
             !
          enddo
          !
          do isym = 1, nsymq
             phase (isym) = phase (isym) * term (2, isym)
          enddo
          !
       enddo
       !
       do isym = 1, nsymq
          phase (isym) = phase (isym) * term (3, isym)
       enddo
       !
    enddo
  
    do is=2,4
      dmagtosym(:,:,:,is) = dmagsym(:,:,:,is-1) / DBLE (nsymq)
    enddo
  
    deallocate (dmags)
    deallocate (dmagsym)
  
    call stop_clock ('hp_sym_dmag')
    return
  end subroutine hp_sym_dmag
  
