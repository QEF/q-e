!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

#include "f_defs.h"

module exx

  USE kinds,    ONLY : DP
  implicit none

  real (DP):: exxalfa=0.d0 ! 1 if exx, 0 elsewhere
  real (DP):: yukawa = 0.d0
  integer:: iunexx
  integer :: exx_nwordwfc
  !
  ! variables defining the auxiliary k-point grid used in X BZ integration
  !
  integer :: nq1=1, nq2=1, nq3=1 ! integers defining the X integration mesh
  integer :: nqs=1               ! number of points in the q-gridd
  integer :: nkqs                ! total number of different k+q
  real (DP), allocatable :: &
             xkq(:,:)            ! xkq(3,nkqs) the auxiliary k+q set
  real (DP), allocatable :: &
             x_occupation(:,:)   ! x_occupation(nbnd,nks) the weight of 
                                 ! auxiliary functions in the density matrix
  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  ! 
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  integer, allocatable :: index_xkq(:,:) ! index_xkq(nks,nqs) 
  integer, allocatable :: index_xk(:)    ! index_xk(nkqs)  
  integer, allocatable :: index_sym(:)   ! index_sym(nkqs)

  real (DP) :: exxdiv ! 1 if exx, 0 elsewhere
  logical   :: x_gamma_extrapolation =.true.
  logical   :: on_double_grid =.false.
  real (DP) :: grid_factor = 8.d0/7.d0 ! 
  real (DP) :: eps =1.d-6


contains
  !------------------------------------------------------------------------
  subroutine init_h_wfc()
  !------------------------------------------------------------------------

     USE kinds,      ONLY : DP
     USE constants,  ONLY : pi
     USE io_files,   ONLY : iunigk, nwordwfc, iunwfc
     USE klist,      ONLY : xk, nks, wk, ngk
     USE wvfct,      ONLY : nbnd, npw, npwx, igk, g2kin, wg
     USE gvect,      ONLY : g
     USE cell_base,  ONLY : tpiba2, omega
     USE wavefunctions_module, ONLY : evc
     USE buffers,    ONLY : save_buffer

     implicit none
     integer :: ik, ig
     REAL(DP) :: alpha, norm

     wg = 0.d0
     wg(1,1:nks/2) = 1.d0 * wk(1:nks/2)

     IF ( nks > 1 ) REWIND( iunigk )
     do ik=1,nks
        IF ( nks > 1 ) READ( iunigk ) igk
        npw=ngk(ik)
        DO ig = 1, npw
           !
           g2kin(ig) = ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                       ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                       ( xk(3,ik) + g(3,igk(ig)) )**2
        END DO
        g2kin(:) = g2kin(:) * tpiba2
        alpha = 1.d0

        norm = 0.d0
        DO ig = 1, npw
           evc(ig,1) =  sqrt(alpha**3/pi) * 8.d0 * pi * alpha /  &
                        (alpha**2 + g2kin(ig))**2 / sqrt(omega)
!           norm = norm + abs(evc(ig,1))**2
        end do
!        write (*,*) "NORM  ", ik, norm
!        evc(:,1) = evc(:,1)/sqrt(norm)

        CALL save_buffer ( evc, nwordwfc, iunwfc, ik)

     END DO

     return
  end subroutine init_h_wfc
  
  !------------------------------------------------------------------------
  subroutine exx_grid_init()
  !------------------------------------------------------------------------
  USE symme,     ONLY : nsym, s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE klist,     ONLY : xk
  USE wvfct,     ONLY : nbnd
  USE io_global, ONLY : stdout
!
! parallel stuff
!
  USE mp_global,  ONLY : nproc, npool
  USE klist,      ONLY : nkstot 

  implicit none
  integer :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
  integer, allocatable :: temp_index_xk(:), temp_index_sym(:)
  integer, allocatable :: temp_index_ikq(:), new_ikq(:)
  real (DP), allocatable :: temp_xkq(:,:)
  logical:: xk_not_found
  real (DP) :: sxk(3), dxk(3), xk_cryst(3)
  real (DP) :: dq1, dq2, dq3
  logical:: no_pool_para

  call start_clock ('exx_grid')

!  write (stdout,'(a)') " EXX: unshifted q-grid "
  write (stdout,'(a,3i4)') " EXX : q-grid dimensions are ",nq1,nq2,nq3
  if (x_gamma_extrapolation) then
     write (stdout,'(a)') " EXX : q->0 dealt with 8/7 -1/7 trick"
     grid_factor = 8.d0/7.d0
  else 
     write (stdout,'(a)') " EXX : q->0 term not estimated "
     grid_factor = 1.d0
  end if

  nqs = nq1 * nq2 * nq3
  !
  ! all processors need to have access to all k+q points
  !
  no_pool_para =  (nq1*nq2*nq3 /= 1) .and. (npool>nspin) .and. ( nsym > 1 )
  if (no_pool_para ) then
     write(stdout,'(5x,a)') &
         '(nq1*nq2*nq3 /= 1) .and. (npool>nspin) .and. ( nsym > 1 )'
     call errore('exx_grid',&
              'pool parallelization not possible in this case', npool)
  end if
  !
  ! set a safe limit as the maximum number of auxiliary points we may need
  ! and allocate auxiliary arrays
  !
  max_nk = nkstot * min(48, 2 * nsym)
  allocate ( temp_index_xk(max_nk), temp_index_sym(max_nk) )
  allocate ( temp_index_ikq(max_nk), new_ikq(max_nk) )
  allocate ( temp_xkq(3,max_nk) )
  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
  temp_nkqs = 0
  do isym=1,nsym
     do ik =1, nkstot
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
        sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                 s(:,2,isym)*xk_cryst(2) + &
                 s(:,3,isym)*xk_cryst(3)
        ! add sxk to the auxiliary list if it is not already present
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) = isym 
        end if

        sxk(:) = - sxk(:)
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) =-isym 
        end if

     end do
  end do

  !
  ! define the q-mesh step-sizes
  !
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)
  !
  ! allocate and fill the array index_xkq(nkstot,nqs)
  !
  if ( nspin == 2 ) then
     allocate ( index_xkq(2*nkstot,nqs) )
     allocate ( x_occupation(nbnd,2*nkstot) )
  else
     allocate ( index_xkq(nkstot,nqs) )
     allocate ( x_occupation(nbnd,nkstot) )
  end if
  nkqs = 0
  new_ikq(:) = 0
  do ik=1,nkstot 
     xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1

            xk_not_found = .true.
            do ikq=1, temp_nkqs
               if ( xk_not_found ) then
                  dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                  if ( abs(dxk(1)).le.eps .and. &
                       abs(dxk(2)).le.eps .and. &
                       abs(dxk(3)).le.eps ) then

                       xk_not_found = .false.

                       if( new_ikq(ikq) == 0) then
                           nkqs = nkqs + 1
                           temp_index_ikq(nkqs) = ikq
                           new_ikq(ikq) = nkqs
                       end if
                       index_xkq(ik,iq) = new_ikq(ikq)

                  end if
               end if
            end do
            if (xk_not_found) then
               write (*,*) ik, iq, temp_nkqs
               write (*,*) sxk(:)
               call errore('exx_grid_init', &
                           ' k + q is not an S*k ', (ik-1) * nqs + iq )
            end if

         end do
       end do
     end do

  end do
  !
  ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
  !
  if ( nspin == 2 ) then
     allocate ( xkq(3,2*nkqs), index_xk(2*nkqs), index_sym(2*nkqs) )
  else
     allocate ( xkq(3,nkqs), index_xk(nkqs), index_sym(nkqs) )
  end if

  do ik =1, nkqs
     ikq = temp_index_ikq(ik)
     xkq(:,ik) = bg(:,1)*temp_xkq(1,ikq) + &
                 bg(:,2)*temp_xkq(2,ikq) + &
                 bg(:,3)*temp_xkq(3,ikq)
     index_xk(ik)  = temp_index_xk(ikq)
     index_sym(ik) = temp_index_sym(ikq)
  end do
  !
  ! clean up
  !
  deallocate (temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
  !
  ! check that everything is what it should be
  !
  call exx_grid_check () 
  !
  call stop_clock ('exx_grid')
  !
  return
  end subroutine exx_grid_init

  !------------------------------------------------------------------------
  subroutine exx_grid_check ( )
  !------------------------------------------------------------------------
  USE symme,     ONLY : nsym, s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : stdout
  USE klist
  implicit none
  real (DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  integer :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
  real (DP) :: eps, dq1, dq2, dq3
  eps = 1.d-6
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)

  do ik =1, nkstot
     xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1
            
            ikq  = index_xkq(ik,iq) 
            ikk  = index_xk(ikq)
            isym = index_sym(ikq)

            xkk_cryst(:)=at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk)+at(3,:)*xk(3,ikk)
            if (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
            isym = abs (isym)
            dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                     s(:,2,isym)*xkk_cryst(2) + &
                     s(:,3,isym)*xkk_cryst(3) - sxk(:)
            dxk(:) = dxk(:) - nint(dxk(:))
            if ( .not. ( abs(dxk(1)).le.eps .and. &
                         abs(dxk(2)).le.eps .and. &
                         abs(dxk(3)).le.eps )   ) then
                 write(*,*) ik,iq
                 write(*,*) ikq,ikk,isym
                 write(*,*) dxk(:)
                 call errore('exx_grid_check', &
                             'something wrong', 1 )
            end if

         end do
       end do
     end do
  end do

  write (stdout,*) ' EXX GRID CHECK SUCCESSFUL '

  return

  end subroutine exx_grid_check

  !------------------------------------------------------------------------
  subroutine exxinit()
  !------------------------------------------------------------------------

    !This subroutine is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : find_free_unit, nwordwfc, &
                                     prefix, tmp_dir, iunwfc, iunigk
    USE io_global,            ONLY : stdout
    USE buffers,              ONLY : get_buffer
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
    USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                     nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, gamma_only
    USE klist,                ONLY : wk, ngk, nks
    USE symme,                ONLY : nsym, s, ftau

    use mp_global,            ONLY : nproc_pool, me_pool
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active

    implicit none
    integer :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    integer :: h_ibnd, half_nbnd
    COMPLEX(DP),allocatable :: temppsic(:), psic(:), tempevc(:,:)
#ifdef __PARA
    integer nxxs
    COMPLEX(DP),allocatable :: temppsic_all(:), psic_all(:)
#endif
    logical, allocatable :: present(:)
    logical :: exst
    integer, allocatable :: rir(:,:)

    call start_clock ('exxinit')

#ifdef __PARA
    nxxs = nrx1s * nrx2s * nrx3s
    allocate(psic_all(nxxs), temppsic_all(nxxs) )
#endif
    allocate(present(nsym),rir(nrx1s*nrx2s*nrx3s,nsym))
    allocate(temppsic(nrxxs), psic(nrxxs), tempevc( npwx, nbnd ))

    exx_nwordwfc=2*nrxxs
    if (.not.exx_is_active()) then 
       iunexx = find_free_unit()
       call diropn(iunexx,'exx', exx_nwordwfc, exst) 
       exxdiv = exx_divergence() 
       exxalfa = get_exx_fraction()
       write (stdout,*) " ! EXXALFA SET TO ", exxalfa
       call start_exx
    endif

    IF ( nks > 1 ) REWIND( iunigk )

    present(1:nsym) = .false.
    do ikq =1,nkqs
       isym = abs(index_sym(ikq))
       if (.not. present(isym) ) then
          present(isym) = .true.
          if ( mod (s (2, 1, isym) * nr1s, nr2s) .ne.0 .or. &
               mod (s (3, 1, isym) * nr1s, nr3s) .ne.0 .or. &
               mod (s (1, 2, isym) * nr2s, nr1s) .ne.0 .or. &
               mod (s (3, 2, isym) * nr2s, nr3s) .ne.0 .or. &
               mod (s (1, 3, isym) * nr3s, nr1s) .ne.0 .or. &
               mod (s (2, 3, isym) * nr3s, nr2s) .ne.0 ) then
             call errore ('exxinit',' EXX + smooth grid is not working',isym)
          end if
          do ir=1, nrx1s * nrx2s * nrx3s
             rir(ir,isym) = ir
          end do
          do k = 1, nr3s
             do j = 1, nr2s
                do i = 1, nr1s
                   call ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                                  nr1s, nr2s, nr3s, ri, rj , rk )
                   ir =   i + ( j-1)*nrx1s + ( k-1)*nrx1s*nrx2s
                   rir(ir,isym) = ri + (rj-1)*nrx1s + (rk-1)*nrx1s*nrx2s
                end do
             end do
          end do

       end if
    end do

    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer  (tempevc, nwordwfc, iunwfc, ik)
       ELSE
          tempevc(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
       ENDIF

       if (gamma_only) then
          half_nbnd = ( nbnd + 1 )/2
          h_ibnd = 0
          do ibnd =1, nbnd, 2     
             h_ibnd = h_ibnd + 1
             !
             temppsic(:) = ( 0.D0, 0.D0 )
             !
             if (ibnd < nbnd) then
                temppsic(nls (igk(1:npw))) = tempevc(1:npw,ibnd)  &
                          + ( 0.D0, 1.D0 ) * tempevc(1:npw,ibnd+1)
                temppsic(nlsm(igk(1:npw))) = CONJG( tempevc(1:npw,ibnd) ) &
                          + ( 0.D0, 1.D0 ) * CONJG( tempevc(1:npw,ibnd+1) )
             else
                temppsic(nls (igk(1:npw))) = tempevc(1:npw,ibnd) 
                temppsic(nlsm(igk(1:npw))) = CONJG( tempevc(1:npw,ibnd) ) 
             end if
             CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )

             do ikq=1,nkqs
                if (index_xk(ikq) .ne. ik) cycle

                isym = abs(index_sym(ikq) )
#ifdef __PARA
                call cgather_smooth(temppsic,temppsic_all)
                IF ( me_pool == 0 ) &
                     psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                call cscatter_smooth(psic_all,psic)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                if (index_sym(ikq) < 0 ) &
                   call errore('exxinit','index_sym < 0 with gamma_only (!?)',1)

                CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*half_nbnd+h_ibnd,1)
             end do
          end do
       else
          do ibnd =1, nbnd     
             temppsic(:) = ( 0.D0, 0.D0 )
             temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
             CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )

             do ikq=1,nkqs
                if (index_xk(ikq) .ne. ik) cycle

                isym = abs(index_sym(ikq) )
#ifdef __PARA
                call cgather_smooth(temppsic,temppsic_all)
                IF ( me_pool == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                call cscatter_smooth(psic_all,psic)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                if (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))

                CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*nbnd+ibnd,1)
             end do
          end do
       end if
    end do
    deallocate(temppsic, psic, tempevc)
    deallocate(present,rir)
#ifdef __PARA
    deallocate(temppsic_all, psic_all)
#endif 

    ! set appropriately the x_occupation
    do ik =1,nks
       x_occupation(1:nbnd,ik) = wg (1:nbnd, ik) / wk(ik) 
    end do

    call stop_clock ('exxinit')  

  end subroutine exxinit

  !-----------------------------------------------------------------------
  subroutine vexx(lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
    !This routine calculates V_xx \Psi

    ! ... This routine computes the product of the Hamiltonian
    ! ... matrix with m wavefunctions contained in psi
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi, spsi, hpsi
    ! ...    n     true dimension of psi, spsi, hpsi
    ! ...    m     number of states psi
    ! ...    psi
    !
    ! ... output:
    ! ...    hpsi  Vexx*psi
    !
    USE constants, ONLY : fpi, e2
    USE cell_base, ONLY : alat, omega, bg, at
    USE symme,     ONLY : nsym, s
    USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm
    USE gsmooth,   ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                           nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, current_k, gamma_only
    USE klist,     ONLY : xk
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl

    implicit none
    INTEGER          :: lda, n, m
    COMPLEX(DP) :: psi(lda,m) 
    COMPLEX(DP) :: hpsi(lda,m)

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP), allocatable :: rhoc(:), vc(:)
    real (DP),   allocatable :: fac(:)
    integer          :: ibnd, ik, im , ig, ikq, iq, isym
    integer          :: h_ibnd, half_nbnd
    real(DP) :: x1, x2
    real(DP) :: tpiba2, qq, xk_cryst(3), sxk(3), xkq(3), x, q(3)

    call start_clock ('vexx')

    allocate (tempphic(nrxxs), temppsic(nrxxs), result(nrxxs), &
              rhoc(nrxxs), vc(nrxxs), fac(ngm) )

    tpiba2 = (fpi / 2.d0 / alat) **2

    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n

    do im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       if (gamma_only) temppsic(nlsm(igk(1:npw))) = CONJG(psi(1:npw,im))
       CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
       
       result(:) = (0.d0,0.d0)

       do iq = 1, nqs
          ikq  = index_xkq(current_k,iq)
          ik   = index_xk(ikq)
          isym = abs(index_sym(ikq))

          xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
          if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3) 
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

          do ig=1,ngm
             q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
             q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
             q(3)= xk(3,current_k) - xkq(3) + g(3,ig)
             qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) ) 
             if (x_gamma_extrapolation) then
                on_double_grid = .true.
                x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
             end if             

             if (qq.gt.1.d-8) then
                fac(ig)=e2*fpi/(tpiba2*qq + yukawa ) * grid_factor
                if (on_double_grid) fac(ig) = 0.d0
             else
                fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
 !                         - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
                if (yukawa .gt. 1.d-8 .and. .not.x_gamma_extrapolation) then
                   fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
                end if
             end if
          end do
          if (gamma_only) then
             half_nbnd = ( nbnd + 1 ) / 2
             h_ibnd = 0
             do ibnd=1,nbnd, 2 !for each band of psi
                h_ibnd = h_ibnd + 1
                x1 = x_occupation(ibnd,  ik)
                if (ibnd < nbnd) then
                   x2 = x_occupation(ibnd + 1,ik)
                else
                   x2 = 0.d0
                end if
                if ( abs(x1) < 1.d-6 .and.  abs(x2) < 1.d-6 ) cycle
                !
                !loads the phi from file
                !
                CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
                                    (ikq-1)*half_nbnd+h_ibnd, -1 )
                !calculate rho in real space
                rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                !brings it to G-space
                CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
                vc(:) = ( 0.D0, 0.D0 )
                vc(nls(1:ngm))  = fac(1:ngm) * rhoc(nls(1:ngm))
                vc(nlsm(1:ngm)) = fac(1:ngm) * rhoc(nlsm(1:ngm))
                !brings back v in real space
                CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

                vc = CMPLX( x1 * DBLE (vc), x2 * AIMAG(vc) )/ nqs

                !accumulates over bands and k points
                result(1:nrxxs) = result(1:nrxxs) + &
                                  DBLE( vc(1:nrxxs) * tempphic(1:nrxxs) )
             end do
          else
             do ibnd=1,nbnd !for each band of psi
                if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                !
                !loads the phi from file
                !
                CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
                                        (ikq-1)*nbnd+ibnd, -1 )
                !calculate rho in real space
                rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                !brings it to G-space
                CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
                vc(:) = ( 0.D0, 0.D0 )
                vc(nls(1:ngm)) = fac(1:ngm) * rhoc(nls(1:ngm))
                if (gamma_only) vc(nlsm(1:ngm)) = fac(1:ngm) * rhoc(nlsm(1:ngm))
                vc = vc * x_occupation(ibnd,ik) / nqs

                !brings back v in real space
                CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

                !accumulates over bands and k points
                result(1:nrxxs)=result(1:nrxxs)+vc(1:nrxxs)*tempphic(1:nrxxs)
             end do
          end if
       end do

       !brings back result in G-space
       CALL cft3s( result, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
       !adds it to hpsi
       hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
    end do

    deallocate (tempphic, temppsic, result, rhoc, vc, fac )

    call stop_clock ('vexx')

     end subroutine vexx

  function exxenergy ()
    ! This function is called to correct the deband value and have 
    ! the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,    ONLY : get_buffer
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, current_k, gamma_only
    USE gvect,      ONLY : gstart
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk
    USE klist,      ONLY : ngk, nks

    implicit none
    REAL (DP)   :: exxenergy,  energy
    INTEGER          :: ibnd, ik
    COMPLEX(DP) :: vxpsi ( npwx, nbnd ), psi(npwx,nbnd)
    COMPLEX(DP) :: ZDOTC

    call start_clock ('exxenergy')

    energy=0.d0
    IF ( nks > 1 ) REWIND( iunigk )
    do ik=1,nks
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer  (psi, nwordwfc, iunwfc, ik)
       ELSE
          psi(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
       END IF
       vxpsi(:,:) = (0.d0, 0.d0)
       call vexx(npwx,npw,nbnd,psi,vxpsi)
       do ibnd=1,nbnd
          energy = energy + &
                   wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
       if (gamma_only .and. gstart == 2) then
           do ibnd=1,nbnd
              energy = energy - &
                       0.5d0 * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd)
           end do
       end if
    end do

    if (gamma_only) energy = 2.d0 * energy


    call reduce ( 1, energy)
    call poolreduce(1,energy)

    exxenergy = energy

    call stop_clock ('exxenergy')
  end function exxenergy

  !-----------------------------------------------------------------------
  function exxenergy2()
  !-----------------------------------------------------------------------
    !
    USE constants, ONLY : fpi, e2
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at
    USE symme,     ONLY : nsym, s
    USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm
    USE gsmooth,   ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                          nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k, gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl

    implicit none
    REAL (DP)   :: exxenergy2,  energy

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:)
    COMPLEX(DP), allocatable :: rhoc(:)
    real (DP),   allocatable :: fac(:)
    integer          :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    integer          :: half_nbnd, h_ibnd
    real(DP)    :: x1, x2
    real(DP) :: tpiba2, qq, xk_cryst(3), sxk(3), xkq(3), vc, x, q(3)

    call start_clock ('exxen2')


    energy=0.d0

    tpiba2 = (fpi / 2.d0 / alat) **2

    allocate (tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )

    IF ( nks > 1 ) REWIND( iunigk )
    do ikk=1,nks
       current_k = ikk
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF

       do jbnd=1, nbnd !for each band of psi (the k cycle is outside band)
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
          if(gamma_only) temppsic(nlsm(igk(1:npw))) = CONJG(evc(1:npw,jbnd))

          CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
       
          do iq = 1, nqs
             ikq  = index_xkq(current_k,iq)
             ik   = index_xk(ikq)
             isym = abs(index_sym(ikq))

             xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
             if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
             sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                      s(:,2,isym)*xk_cryst(2) + &
                      s(:,3,isym)*xk_cryst(3) 
             xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

             do ig=1,ngm
                q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                q(3)= xk(3,current_k) - xkq(3) + g(3,ig)
                qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )
                if (x_gamma_extrapolation) then
                   on_double_grid = .true.
                   x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                   x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                   x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                end if             
                if (qq.gt.1.d-8) then
                   fac(ig)=e2*fpi/(tpiba2*qq + yukawa ) * grid_factor
                   if (gamma_only) fac(ig) = 2.d0 * fac(ig)
                   if (on_double_grid) fac(ig) = 0.d0
                else
                   fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
 !                            - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
                   if (yukawa.gt.1.d-8 .and. .not. x_gamma_extrapolation) then
                      fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
                   end if
                end if
             end do

             if (gamma_only) then
                half_nbnd = ( nbnd + 1) / 2
                h_ibnd = 0
                do ibnd=1,nbnd, 2 !for each band of psi
                   h_ibnd = h_ibnd + 1
                   x1 = x_occupation(ibnd,ik)
                   if ( ibnd < nbnd ) then
                      x2 = x_occupation(ibnd+1,ik)
                   else
                      x2 = 0.d0
                   end if
                   if ( abs(x1) < 1.d-6 .and. abs(x2) < 1.d-6 ) cycle
                   !
                   !loads the phi from file
                   !
                   CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                                          (ikq-1)*half_nbnd+h_ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
                   vc = 0.D0
                   do ig=1,ngm
                      vc = vc + fac(ig) * x1 * &
                                abs( rhoc(nls(ig)) + CONJG(rhoc(nlsm(ig))) )**2
                      vc = vc + fac(ig) * x2 * &
                                abs( rhoc(nls(ig)) - CONJG(rhoc(nlsm(ig))) )**2
                   end do
                   vc = vc * omega * 0.25d0 / nqs

                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do
             else
                do ibnd=1,nbnd !for each band of psi
                   if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                   !
                   !loads the phi from file
                   !
                   CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                                          (ikq-1)*nbnd+ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
                   vc = 0.D0
                   do ig=1,ngm
                      vc = vc + fac(ig) * rhoc(nls(ig)) * CONJG(rhoc(nls(ig)))
                   end do
                   vc = vc * omega * x_occupation(ibnd,ik) / nqs

                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do
             end if
          end do
       end do
    end do

    deallocate (tempphic, temppsic, rhoc, fac )

    call reduce ( 1, energy)
    call poolreduce(1,energy)

    exxenergy2 = energy

    call stop_clock ('exxen2')

  end function  exxenergy2

  function exx_divergence ()

     USE constants, ONLY : fpi, e2
     USE cell_base, ONLY : bg, at, alat, omega
     USE gvect,     ONLY : ngm, g, ecutwfc
     USE wvfct,     ONLY : gamma_only
     USE io_global, ONLY : stdout

     implicit none
     real(DP) :: exx_divergence

     ! local variables
     integer :: iq1,iq2,iq3, ig
     real(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)

     integer :: nqq, iq
     real(DP) :: aa, dq

     call start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10.d0 * tpiba2 / ecutwfc

     dq1= 1.d0/DBLE(nq1)
     dq2= 1.d0/DBLE(nq2) 
     dq3= 1.d0/DBLE(nq3) 

     div = 0.d0
     do iq1=1,nq1
        do iq2=1,nq2
           do iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              do ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) ) 
                 if (x_gamma_extrapolation) then
                    on_double_grid = .true.
                    x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 end if
                 if (.not.on_double_grid) then
                    if ( qq.gt.1.d-8 .or. yukawa .gt. 1.d-8) then
                       div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                    else
                       div = div - alpha ! or maybe something else
                    end if
                 end if
              end do
           end do
        end do
     end do
     call reduce (1, div)
     if (gamma_only) then
        div = 2.d0 * div
        if ( .not. x_gamma_extrapolation ) then
           if (yukawa .le. 1.d-8) then
              div = div + alpha
           else
              div = div - tpiba2/yukawa
           end if
        end if
     end if

     div = div * e2 * fpi / tpiba2 / nqs

     alpha = alpha / tpiba2

     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0.d0
     do iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
     end do
     aa = aa * 8.d0 /fpi
     aa = aa + 1.d0/sqrt(alpha*0.25d0*fpi) 
     write (stdout,*) aa, 1.d0/sqrt(alpha*0.25d0*fpi)
    
     div = div - e2*omega * aa

!     div = div - e2*omega/sqrt(alpha*0.25d0*fpi)

     exx_divergence = div * nqs

     write (stdout,'(a,i4,a,3f12.4)') 'EXX divergence (',nq1,')= ', &
                                  div, alpha

     call stop_clock ('exx_div')

     call print_clock ('exx_div')

     return
  end function exx_divergence 

end module exx
