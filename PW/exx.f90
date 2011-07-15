!
! Copyright (C) 2005-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx
  !--------------------------------------
  !
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE io_global,            ONLY : ionode, stdout
  !
  IMPLICIT NONE
  SAVE

  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0.d0                ! 1 if exx, 0 elsewhere
  INTEGER :: exx_nwordwfc, ji

  !
  ! variables defining the auxiliary k-point grid 
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-gridd
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq(:,:)   ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)           
                                         ! x_occupation(nbnd,nks) the weight of 
                                         ! auxiliary functions in the density matrix
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporay buffer to store wfc 

  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  ! 
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  INTEGER, ALLOCATABLE :: index_xkq(:,:) ! index_xkq(nks,nqs) 
  INTEGER, ALLOCATABLE :: index_xk(:)    ! index_xk(nkqs)  
  INTEGER, ALLOCATABLE :: index_sym(:)   ! index_sym(nkqs)

  !
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL (DP)         :: eps =1.d-6
  REAL (DP)         :: exxdiv = 0.0
  CHARACTER(80)     :: exxdiv_treatment 
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL (DP)         :: grid_factor = 8.d0/7.d0 
  !
  ! Gygi-Baldereschi 
  LOGICAL           :: use_regularization = .TRUE.
  !
  ! yukawa method
  LOGICAL           :: use_yukawa = .FALSE.
  REAL (DP)         :: yukawa = 0.d0
  !
  ! erfc screening
  LOGICAL           :: use_erfc_simple_div = .FALSE.
  REAL (DP)         :: erfc_scrlen = 0.d0
  !
  ! cutoff techniques
  LOGICAL           :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL           :: use_coulomb_vcut_spheric = .FALSE.
  REAL (DP)         :: ecutvcut
  TYPE(vcut_type)   :: vcut

  !
  ! energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP, & !   sum <phi|Vx(phi)|phi>
              fock1 = 0.0_DP, & !   sum <psi|vx(phi)|psi>
              fock2 = 0.0_DP, & !   sum <psi|vx(psi)|psi>
              dexx  = 0.0_DP    !   fock1  - 0.5*(fock2+fock0)

CONTAINS

  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
  !------------------------------------------------------------------------
  IF ( ALLOCATED (index_xkq) ) DEALLOCATE (index_xkq)
  IF ( ALLOCATED (index_xk ) ) DEALLOCATE (index_xk )
  IF ( ALLOCATED (index_sym) ) DEALLOCATE (index_sym)
  IF ( ALLOCATED (x_occupation) ) DEALLOCATE (x_occupation)
  IF ( ALLOCATED (xkq) ) DEALLOCATE (xkq)
  IF ( ALLOCATED (exxbuff) ) DEALLOCATE (exxbuff)
  !
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  subroutine exx_grid_init()
  !------------------------------------------------------------------------
  !
  USE symm_base,  ONLY : nsym, s
  USE cell_base,  ONLY : bg, at, alat
  USE lsda_mod,   ONLY : nspin
  USE klist,      ONLY : xk
  USE wvfct,      ONLY : nbnd
  USE io_global,  ONLY : stdout
  !
  USE mp_global,  ONLY : nproc, npool, nimage
  USE klist,      ONLY : nkstot 
  !
  IMPLICIT NONE
  !
  CHARACTER(13) :: sub_name='exx_grid_init'
  integer       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
  integer,   allocatable :: temp_index_xk(:), temp_index_sym(:)
  integer,   allocatable :: temp_index_ikq(:), new_ikq(:)
  real (DP), allocatable :: temp_xkq(:,:)
  logical       :: xk_not_found
  real (DP)     :: sxk(3), dxk(3), xk_cryst(3)
  real (DP)     :: dq1, dq2, dq3
  logical       :: no_pool_para

  CALL start_clock ('exx_grid')

  !
  ! definitions and checks
  !
#ifdef EXXDEBUG
  IF (ionode) WRITE(stdout,'(/,2x,a,3i4)') "EXX : q-grid dimensions are ", nq1,nq2,nq3
#endif
  !
  grid_factor = 1.d0
  !
  IF (x_gamma_extrapolation) THEN
      !
#ifdef EXXDEBUG
      IF (ionode) WRITE (stdout,'(2x,a)') "EXX : q->0 dealt with 8/7 -1/7 trick"
#endif
      grid_factor = 8.d0/7.d0
      !
  ENDIF
 
  !
  ! parallelism over q
  !
  nqs = nq1 * nq2 * nq3
!
! old check for parallelization on images  
!  IF( MOD(nqs,nimage) /= 0 ) THEN
!      CALl errore(sub_name, 'The total number of q points must be multiple of nimage', mod(nqs,nimage))
!  ENDIF
!
  !
  ! all processors need to have access to all k+q points
  !
  no_pool_para =  (nq1*nq2*nq3 /= 1) .and. (npool>nspin) .and. ( nsym > 1 )
  if (no_pool_para ) then
     write(stdout,'(5x,a)') &
         '(nq1*nq2*nq3 /= 1) .and. (npool>nspin) .and. ( nsym > 1 )'
     call errore(sub_name,&
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
     if(.not.ALLOCATED(index_xkq)) allocate ( index_xkq(2*nkstot,nqs) )
     if(.not.ALLOCATED(x_occupation)) allocate ( x_occupation(nbnd,2*nkstot) )
  else
     if(.not.ALLOCATED(index_xkq)) allocate ( index_xkq(nkstot,nqs) )
     if(.not.ALLOCATED(x_occupation)) allocate ( x_occupation(nbnd,nkstot) )
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

  CALL stop_clock ('exx_grid')
  !
  RETURN
  END SUBROUTINE exx_grid_init

  !------------------------------------------------------------------------
  subroutine exx_div_check()
  !------------------------------------------------------------------------
  !
  USE cell_base,  ONLY : bg, at, alat
  USE io_global,  ONLY : stdout
  USE funct,      ONLY : get_screening_parameter
  !
  IMPLICIT NONE
  !
  REAL (DP)     :: atws(3,3)
  CHARACTER(13) :: sub_name='exx_div_check'

  !
  ! EXX singularity treatment
  !
  !
  SELECT CASE ( TRIM(exxdiv_treatment) ) 
  CASE ( "gygi-baldereschi", "gygi-bald", "g-b" )
     !
     use_regularization = .TRUE.
     !
  CASE ( "yukawa" ) 
     !
     use_yukawa = .TRUE.
     IF ( yukawa <= 0.0 ) CALL errore(sub_name,'invalid yukawa parameter', 1)
     !
  CASE ( "erfc_simple" )
     !
     use_erfc_simple_div = .TRUE.
     erfc_scrlen = get_screening_parameter()
     write(0,*) "erfc_scrlen", erfc_scrlen
     IF ( erfc_scrlen <= 0.0 ) CALL errore(sub_name,'invalid screening length', 1)
     !
  CASE ( "vcut_ws" )
     !
     use_coulomb_vcut_ws = .TRUE.
     IF ( x_gamma_extrapolation ) &
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_ws', 1)
     !
  CASE ( "vcut_spheric" ) 
     !
     use_coulomb_vcut_spheric = .TRUE.
     IF ( x_gamma_extrapolation ) &
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_spheric', 1)
     !
  CASE ( "none" )
     use_regularization = .FALSE.
     !
  CASE DEFAULT
     CALL errore(sub_name,'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1)
  END SELECT
  !
#ifdef EXXDEBUG
  IF ( ionode ) WRITE (stdout,'(2x,"EXX : q->0 dealt with ",a, " trick" )') &
                TRIM(exxdiv_treatment)
#endif
  !
  ! <AF>
  ! Set variables for Coulomb vcut
  ! NOTE: some memory is allocated inside this routine (in the var vcut)
  !       and should be deallocated somewehre, at the end of the run
  !
  IF ( use_coulomb_vcut_ws .OR. use_coulomb_vcut_spheric ) THEN
      !
      ! build the superperiodicity direct lattice
      !
      atws = alat * at
      !
      atws(:,1) = atws(:,1) * nq1
      atws(:,2) = atws(:,2) * nq2
      atws(:,3) = atws(:,3) * nq3
      !
      !CALL start_clock ('exx_vcut_init')
      CALL vcut_init( vcut, atws, ecutvcut )
      !CALL stop_clock ('exx_vcut_init')
      !
      IF ( ionode ) CALL vcut_info( stdout, vcut )
      !          
  ENDIF
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : exx div treatment check successful')")
#endif 
  RETURN
  END SUBROUTINE exx_div_check 


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( )
  !------------------------------------------------------------------------
  USE symm_base, ONLY : s
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
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : grid check successful')")
#endif
  return

  end subroutine exx_grid_check

  !------------------------------------------------------------------------
  subroutine exx_restart(l_exx_was_active)
  !------------------------------------------------------------------------
    !This subroutine is called when restarting an exx calculation
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter
    USE fft_base,             ONLY : dffts
    USE io_global,            ONLY : stdout

    implicit none
    logical, intent(in) :: l_exx_was_active
    logical :: exst

    if (.not. l_exx_was_active ) return ! nothing had happpened yet
    !!
    exx_nwordwfc=2*dffts%nnr
    !iunexx = find_free_unit()
    !call diropn(iunexx,'exx', exx_nwordwfc, exst) 
    erfc_scrlen = get_screening_parameter()
    exxdiv = exx_divergence() 
    exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
    write (stdout,*) " ! EXXALFA SET TO ", exxalfa
#endif
    call start_exx
    call weights()
    call exxinit()
    fock0 = exxenergy2()
 
    return
  end subroutine exx_restart
  !------------------------------------------------------------------------
  subroutine exxinit()
  !------------------------------------------------------------------------

    !This subroutine is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    !
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, &
                                     tmp_dir, prefix
    USE io_global,            ONLY : stdout
    USE buffers,              ONLY : get_buffer
    USE gvecs,              ONLY : nls, nlsm, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : wk, ngk, nks
    USE symm_base,            ONLY : nsym, s, ftau

    use mp_global,            ONLY : nproc_pool, me_pool
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter 
    use fft_base,             ONLY : cgather_smooth, cscatter_smooth, dffts
    use fft_interfaces,       ONLY : invfft

    implicit none
    integer :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    integer :: h_ibnd, half_nbnd
    COMPLEX(DP),allocatable :: temppsic(:), psic(:), tempevc(:,:)
    integer :: nxxs, nrxxs
#ifdef __PARA
    COMPLEX(DP),allocatable :: temppsic_all(:), psic_all(:)
#endif
    logical, allocatable :: present(:)
    logical :: exst
    integer, allocatable :: rir(:,:)

    call start_clock ('exxinit')

    ! Beware: not the same as nrxxs in parallel case
    nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
    nrxxs= dffts%nnr
#ifdef __PARA
    allocate(psic_all(nxxs), temppsic_all(nxxs) )
#endif
    allocate(present(nsym),rir(nxxs,nsym))
    allocate(temppsic(nrxxs), psic(nrxxs),tempevc( npwx, nbnd ))

    if( .not. allocated( exxbuff ) ) allocate( exxbuff( nrxxs, nkqs, nbnd ) )

    !write(*,*) 'debug: exxbuff size=',size(exxbuff)
    exx_nwordwfc=2*nrxxs
    if (.not.exx_is_active()) then 
       !iunexx = find_free_unit()
       !call diropn(iunexx,'exx', exx_nwordwfc, exst) 
       erfc_scrlen = get_screening_parameter()
       exxdiv = exx_divergence() 
       exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
       write (stdout,*) " ! EXXALFA SET TO ", exxalfa
#endif
       call start_exx
    endif

    IF ( nks > 1 ) REWIND( iunigk )

    present(1:nsym) = .false.
    do ikq =1,nkqs
       isym = abs(index_sym(ikq))
       if (.not. present(isym) ) then
          present(isym) = .true.
          if ( mod (s (2, 1, isym) * dffts%nr1, dffts%nr2) .ne.0 .or. &
               mod (s (3, 1, isym) * dffts%nr1, dffts%nr3) .ne.0 .or. &
               mod (s (1, 2, isym) * dffts%nr2, dffts%nr1) .ne.0 .or. &
               mod (s (3, 2, isym) * dffts%nr2, dffts%nr3) .ne.0 .or. &
               mod (s (1, 3, isym) * dffts%nr3, dffts%nr1) .ne.0 .or. &
               mod (s (2, 3, isym) * dffts%nr3, dffts%nr2) .ne.0 ) then
             call errore ('exxinit',' EXX + smooth grid is not working',isym)
          end if
          do ir=1, nxxs
             rir(ir,isym) = ir
          end do
          do k = 1, dffts%nr3
             do j = 1, dffts%nr2
                do i = 1, dffts%nr1
                   call ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                                  dffts%nr1,dffts%nr2,dffts%nr3, ri, rj , rk )
                   ir =   i + ( j-1)*dffts%nr1x + ( k-1)*dffts%nr1x*dffts%nr2x
                   rir(ir,isym) = ri + (rj-1)*dffts%nr1x + (rk-1)*dffts%nr1x*dffts%nr2x
                end do
             end do
          end do

       end if
    end do

    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (tempevc, nwordwfc, iunwfc, ik)
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
             CALL invfft ('Wave', temppsic, dffts)

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

                do ji=1, nrxxs
                   exxbuff(ji,ikq,h_ibnd)=psic(ji)
                enddo

                !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*half_nbnd+h_ibnd,1)
             end do
          end do
       else
          do ibnd =1, nbnd     
             temppsic(:) = ( 0.D0, 0.D0 )
             temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
             CALL invfft ('Wave', temppsic, dffts)

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
                do ji=1, nrxxs
                   exxbuff(ji,ikq,ibnd)=psic(ji)
                   
                enddo
                !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*nbnd+ibnd,1)
             end do
          end do
       end if
    end do

    deallocate(temppsic, psic,tempevc)
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
  SUBROUTINE vexx(lda, n, m, psi, hpsi)
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
    USE constants, ONLY : fpi, e2, pi
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base, ONLY : nsym, s
    USE gvect,     ONLY : ngm
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, current_k
    USE control_flags, ONLY : gamma_only
    USE klist,     ONLY : xk
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft

    USE parallel_include  
    USE mp_global, ONLY : nproc, inter_image_comm, me_image, my_image_id,&
         & nimage, nproc_image
    USE mp,        ONLY : mp_sum


    IMPLICIT NONE
    INTEGER          :: lda, n, m, nqi, myrank, mysize
    COMPLEX(DP) :: psi(lda,m) 
    COMPLEX(DP) :: hpsi(lda,m)

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP), allocatable :: rhoc(:), vc(:)
    REAL (DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ig, ikq, iq, isym, iqi
    INTEGER          :: h_ibnd, half_nbnd, ierr, nrxxs
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), x, q(3)
    ! <LMS> temp array for vcut_spheric
    REAL(DP) :: atws(3,3)

    CALL start_clock ('vexx')
    nrxxs = dffts%nnr
    ALLOCATE (tempphic(nrxxs), temppsic(nrxxs), result(nrxxs), &
              rhoc(nrxxs), vc(nrxxs), fac(ngm) )

    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n
!
! Was used for parallelization on images
!    nqi=nqs/nimage
    nqi=nqs
!

    DO im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       IF (gamma_only) temppsic(nlsm(igk(1:npw))) = CONJG(psi(1:npw,im))
       CALL invfft ('Wave', temppsic, dffts)
       
       result(:) = (0.d0,0.d0)

       DO iqi=1,nqi
!
! Was used for parallelization on images
!          iq=iqi+nqi*my_image_id
          iq=iqi
!
          ikq  = index_xkq(current_k,iq)
          ik   = index_xk(ikq)
          isym = ABS(index_sym(ikq))

          xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
          IF (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
               s(:,2,isym)*xk_cryst(2) + &
               s(:,3,isym)*xk_cryst(3) 
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)
          !
          ! calculate the 1/|r-r'| factor and place it in fac
          !
          CALL g2_convolution(ngm, g, xk(:,current_k), xkq, fac)
          !
          IF (gamma_only) THEN
             half_nbnd = ( nbnd + 1 ) / 2
             h_ibnd = 0
             DO ibnd=1,nbnd, 2 !for each band of psi
                h_ibnd = h_ibnd + 1
                x1 = x_occupation(ibnd,  ik)
                IF (ibnd < nbnd) THEN
                   x2 = x_occupation(ibnd + 1,ik)
                ELSE
                   x2 = 0.d0
                END IF
                IF ( ABS(x1) < 1.d-6 .AND.  ABS(x2) < 1.d-6 ) CYCLE
                !
                !loads the phi from file
                !
#ifdef EXXDEBUG
WRITE(stdout,*) "from vexx, nrxxs: ", nrxxs, ikq, h_ibnd 
call flush_unit(stdout)
#endif
              DO ji=1, nrxxs
                 tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                   
              ENDDO
              !CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
              !                    (ikq-1)*half_nbnd+h_ibnd, -1 )
              !calculate rho in real space
              rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
              !brings it to G-space
              CALL fwfft ('Smooth', rhoc, dffts)
   
              vc(:) = ( 0.D0, 0.D0 )
              vc(nls(1:ngm))  = fac(1:ngm) * rhoc(nls(1:ngm))
              vc(nlsm(1:ngm)) = fac(1:ngm) * rhoc(nlsm(1:ngm))
              !brings back v in real space
              CALL invfft ('Smooth', vc, dffts) 

              vc = CMPLX( x1 * DBLE (vc), x2 * AIMAG(vc) ,kind=DP)/ nqs

              !accumulates over bands and k points
              result(1:nrxxs) = result(1:nrxxs) + &
                   DBLE( vc(1:nrxxs) * tempphic(1:nrxxs) )

           END DO

        ELSE
           DO ibnd=1,nbnd !for each band of psi
              IF ( ABS(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
              !
              !loads the phi from file
              !
#ifdef EXXDEBUG
WRITE(stdout,*) "from vexx, nrxxs: ", nrxxs, ikq, ibnd
if(allocated(tempphic)) write(stdout,*) "tempphic allcated"
if(allocated(exxbuff)) write(stdout,*) "exxbuff allocated"
call flush_unit(stdout)
#endif
              DO ji=1, nrxxs
                 tempphic(ji)=exxbuff(ji,ikq,ibnd)                
              ENDDO
                
              !CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
              !                        (ikq-1)*nbnd+ibnd, -1 )
              !calculate rho in real space
              rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
              !brings it to G-space
              CALL fwfft ('Smooth', rhoc, dffts)
   
              vc(:) = ( 0.D0, 0.D0 )
              vc(nls(1:ngm)) = fac(1:ngm) * rhoc(nls(1:ngm))
              IF (gamma_only) vc(nlsm(1:ngm)) = fac(1:ngm) * rhoc(nlsm(1:ngm))
              vc = vc * x_occupation(ibnd,ik) / nqs

              !brings back v in real space
              CALL invfft ('Smooth', vc, dffts) 
                
              !accumulates over bands and k points
              result(1:nrxxs)=result(1:nrxxs)+vc(1:nrxxs)*tempphic(1:nrxxs)
           END DO
        END IF
     END DO
     !
     ! Was used for parallelization on images
     !       CALL mp_sum( result(1:nrxxs), inter_image_comm )
     !brings back result in G-space
     CALL fwfft ('Wave', result, dffts)
     !adds it to hpsi
     hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
  END DO
    
  DEALLOCATE (tempphic,temppsic, result, rhoc, vc, fac )
   
  CALL stop_clock ('vexx')

END SUBROUTINE vexx
!-----------------------------------------------------------------------
SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
!-----------------------------------------------------------------------
  ! This routine calculates the 1/|r-r'| part of the exact exchange 
  ! expression in reciprocal space (the G^-2 factor).
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba, at
  USE constants, ONLY : fpi, e2, pi
  
  IMPLICIT NONE
  
  INTEGER,  INTENT(IN)    :: ngm   ! Number of G vectors
  REAL(DP), INTENT(IN)    :: g(3,ngm) ! Cartesian components of G vectors
  REAL(DP), INTENT(IN)    :: xk(3) ! current k vector
  REAL(DP), INTENT(IN)    :: xkq(3) ! current q vector
  
  REAL(DP), INTENT(INOUT) :: fac(ngm) ! Calculated convolution
  
  
  !Local variables
  INTEGER :: ig !Counters 
  REAL(DP) :: q(3), qq, x
  
  !CALL start_clock ('vexx_ngmloop')
  DO ig=1,ngm
     !
     q(:)= xk(:) - xkq(:) + g(:,ig)
     !
     q = q * tpiba
     !
     qq = SUM(q(:)**2.0_DP) 
     !
     IF (x_gamma_extrapolation) THEN
        on_double_grid = .TRUE.
        x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
        x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
        x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
     ENDIF
     
     IF ( use_coulomb_vcut_ws ) THEN
        !
        fac(ig) = vcut_get(vcut,q)
        !
     ELSE IF ( use_coulomb_vcut_spheric ) THEN
        !
        fac(ig) = vcut_spheric_get(vcut,q)
        !
     ELSE IF (qq.GT.1.d-8) THEN
        !
        IF ( erfc_scrlen > 0  ) THEN
           fac(ig)=e2*fpi/qq*(1-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor
        ELSE
           fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
        END IF
        IF (on_double_grid) fac(ig) = 0.d0
        !
     ELSE
        !
        fac(ig)= - exxdiv ! or rather something else (see F.Gygi)
        !
        IF ( use_yukawa .AND. .NOT. x_gamma_extrapolation ) &
             fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
        !
     ENDIF
     !
  ENDDO
  !CALL stop_clock ('vexx_ngmloop')
END SUBROUTINE g2_convolution
!-----------------------------------------------------------------------

  function exxenergy ()
    ! This function is called to correct the deband value and have 
    ! the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,    ONLY : get_buffer
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE gvect,      ONLY : gstart
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk
    USE klist,      ONLY : ngk, nks
    USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
    USE mp,         ONLY : mp_sum

    implicit none
    REAL (DP)   :: exxenergy,  energy
    INTEGER          :: ibnd, ik
    COMPLEX(DP) :: vxpsi ( npwx, nbnd ), psi(npwx,nbnd)
    COMPLEX(DP) :: zdotc

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
                   wg(ibnd,ik) * zdotc(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
       if (gamma_only .and. gstart == 2) then
           do ibnd=1,nbnd
              energy = energy - &
                       0.5d0 * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd)
           end do
       end if
    end do

    if (gamma_only) energy = 2.d0 * energy


    call mp_sum( energy, intra_pool_comm )
    call mp_sum( energy, inter_pool_comm )

    exxenergy = energy

    call stop_clock ('exxenergy')
  end function exxenergy

  !-----------------------------------------------------------------------
  function exxenergy2()
  !-----------------------------------------------------------------------
    !
    USE constants, ONLY : fpi, e2, pi
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base,ONLY : nsym, s
    USE gvect,     ONLY : ngm
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    USE mp_global, ONLY : inter_pool_comm, intra_pool_comm, inter_image_comm
    USE mp_global, ONLY : my_image_id, nimage
    USE mp,        ONLY : mp_sum
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft

    IMPLICIT NONE
    REAL (DP)   :: exxenergy2,  energy

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:)
    COMPLEX(DP), allocatable :: rhoc(:)
    real (DP),   allocatable :: fac(:)
    integer          :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    integer          :: half_nbnd, h_ibnd, nqi, iqi, nrxxs
    real(DP)    :: x1, x2
    real(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc, x, q(3)
    ! temp array for vcut_spheric
    real(DP) :: atws(3,3) 

    call start_clock ('exxen2')

    nrxxs = dffts%nnr
    energy=0.d0

    ALLOCATE (tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
!
! Was used for parallelization on images    
!    nqi=nqs/nimage
    nqi=nqs
!
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

          CALL invfft ('Wave', temppsic, dffts)
       
          do iqi=1,nqi
!
! Was used for parallelization on images
!             iq=iqi+nqi*my_image_id
             iq=iqi
!
             ikq  = index_xkq(current_k,iq)
             ik   = index_xk(ikq)
             isym = abs(index_sym(ikq))

             xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
             if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
             sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                      s(:,2,isym)*xk_cryst(2) + &
                      s(:,3,isym)*xk_cryst(3) 
             xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

             !CALL start_clock ('exxen2_ngmloop')
             DO ig=1,ngm
                !
                q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                q(3)= xk(3,current_k) - xkq(3) + g(3,ig)
                !
                q = q * tpiba
                ! 
                qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )
                !
                IF (x_gamma_extrapolation) THEN
                   on_double_grid = .true.
                   x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                   x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                   x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                   on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                ENDIF

                IF ( use_coulomb_vcut_ws ) THEN
                   !
                   fac(ig) = vcut_get(vcut, q)
                   IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)
                   !
                ELSE IF ( use_coulomb_vcut_spheric ) THEN
                   !
                   fac(ig) = vcut_spheric_get(vcut, q)
                   IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 
                   ! 
                ELSE IF (qq > 1.d-8) THEN
                   !
                   IF ( erfc_scrlen > 0 ) THEN
                      fac(ig)=e2*fpi/qq*(1-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                   ELSE
                      fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                   END IF
                   IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                   IF (on_double_grid) fac(ig) = 0.d0
                   !
                ELSE
                   !
                   fac(ig)= - exxdiv ! or rather something else (see F.Gygi)
                   !
                   IF ( use_yukawa .AND. .NOT. x_gamma_extrapolation) &
                      fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                ENDIF
                !
             ENDDO
             !CALL stop_clock ('exxen2_ngmloop')

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
                   do ji=1, nrxxs
                      tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                   enddo
                
                   !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                   !                       (ikq-1)*half_nbnd+h_ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL fwfft ('Smooth', rhoc, dffts)
   
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
                   do ji=1, nrxxs
                      tempphic(ji)=exxbuff(ji,ikq,ibnd)

                   enddo

                   !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                   !                       (ikq-1)*nbnd+ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL fwfft ('Smooth', rhoc, dffts)
   
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
!
! Was used for image parallelization
!    call mp_sum( energy, inter_image_comm )
!
    call mp_sum( energy, intra_pool_comm )
    call mp_sum( energy, inter_pool_comm )

    exxenergy2 = energy

    call stop_clock ('exxen2')

  end function  exxenergy2

  function exx_divergence ()

     USE constants, ONLY : fpi, e2, pi
     USE cell_base, ONLY : bg, at, alat, omega
     USE gvect,     ONLY : ngm, g
     USE wvfct,     ONLY : ecutwfc
     USE io_global, ONLY : stdout
     USE control_flags, ONLY : gamma_only
     USE mp_global, ONLY : intra_pool_comm
     USE mp,        ONLY : mp_sum

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

     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0.d0
        return
     END IF

     IF ( use_erfc_simple_div .AND. erfc_scrlen > 0 .AND. &
          .NOT. x_gamma_extrapolation ) THEN
        exx_divergence = -e2*pi/(erfc_scrlen**2)
#ifdef EXXDEBUG
        write(stdout,*) 'exx_divergence', exx_divergence
#endif
        return
     END IF

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
                    if ( qq > 1.d-8 .OR. use_yukawa ) then
                       if ( erfc_scrlen > 0 ) then
                          div = div + exp( -alpha * qq) / qq * &
                                (1-exp(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
                       else

                          div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                       endif
                    else
                       div = div - alpha ! or maybe something else
                    end if
                 end if
              end do
           end do
        end do
     end do
     call mp_sum(  div, intra_pool_comm )
     if (gamma_only) then
        div = 2.d0 * div
        if ( .not. x_gamma_extrapolation ) then
           if ( use_yukawa ) then
              div = div - tpiba2/yukawa
           else
              div = div + alpha
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
        if ( erfc_scrlen > 0 ) then
           aa = aa  -exp( -alpha * qq) * exp(-qq/4.d0/erfc_scrlen**2) * dq
        else
           aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
        end if
     end do
     aa = aa * 8.d0 /fpi
     aa = aa + 1.d0/sqrt(alpha*0.25d0*fpi) 
#ifdef EXXDEBUG
     write (stdout,*) aa, 1.d0/sqrt(alpha*0.25d0*fpi)
#endif   
     div = div - e2*omega * aa

!     div = div - e2*omega/sqrt(alpha*0.25d0*fpi)

     exx_divergence = div * nqs
#ifdef EXXDEBUG
     write (stdout,'(a,i4,a,3f12.4)') 'EXX divergence (',nq1,')= ', &
                                  div, alpha
#endif
     call stop_clock ('exx_div')

     return
  end function exx_divergence 
  



  !-----------------------------------------------------------------------
  FUNCTION exx_stress()
  !-----------------------------------------------------------------------
  !
  ! This is Eq.(10) of PRB 73, 125120 (2006).
  !
  USE constants, ONLY : fpi, e2, pi, tpi
  USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
  USE buffers,   ONLY : get_buffer
  USE cell_base, ONLY : alat, omega, bg, at, tpiba
  USE symm_base,ONLY : nsym, s
  USE gvect,     ONLY : ngm
  USE gvecs,   ONLY : nls, nlsm, doublegrid
  USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k
  USE control_flags, ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE klist,     ONLY : xk, ngk, nks
  USE lsda_mod,  ONLY : lsda, current_spin, isk
  USE gvect,     ONLY : g, nl
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp_global, ONLY : my_image_id, nimage
  USE mp,        ONLY : mp_sum 
  use fft_base,  ONLY : dffts
  use fft_interfaces, ONLY : fwfft, invfft
  ! ---- local variables -------------------------------------------------
  IMPLICIT NONE
  real (dp)   :: exx_stress(3,3), exx_stress_(3,3)
  complex(dp), allocatable :: tempphic(:), temppsic(:)
  complex(dp), allocatable :: rhoc(:)
  real(dp), allocatable :: fac(:), fac_tens(:,:,:), fac_stress(:)
  integer :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
  integer :: half_nbnd, h_ibnd, nqi, iqi, beta, nrxxs
  real(dp) :: x1, x2
  real(dp) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
  ! temp array for vcut_spheric
  real(dp) :: atws(3,3) 
  real(dp) :: delta(3,3)
 
  call start_clock ('exx_stress')

  nrxxs = dffts%nnr
  delta = reshape( (/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/), (/3,3/))
  exx_stress_ = 0.d0
  allocate( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
  allocate( fac_tens(3,3,ngm), fac_stress(ngm) )

  if ( nks > 1 ) rewind( iunigk )
!
! Was used for image parallelization
!
  nqi=nqs
!  nqi = nqs/nimage
!
  ! loop over k-points
  do ikk = 1, nks
      current_k = ikk
      if (lsda) current_spin = isk(ikk)
      npw = ngk(ikk)

      if (nks > 1) then
          read(iunigk) igk
          call get_buffer(evc, nwordwfc, iunwfc, ikk)
      end if

      ! loop over bands
      do jbnd = 1, nbnd
          temppsic(:) = ( 0.d0, 0.d0 )
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
          if(gamma_only) temppsic(nlsm(igk(1:npw))) = conjg(evc(1:npw,jbnd))
          CALL invfft ('Wave', temppsic, dffts)       

          do iqi = 1, nqi
!
! Was used for image parallelization
!
!              iq = iqi + nqi*my_image_id
              iq=iqi
!
              ikq  = index_xkq(current_k,iq)
              ik   = index_xk(ikq)
              isym = abs(index_sym(ikq))

              xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
              if (index_sym(ikq) < 0) xk_cryst = -xk_cryst
              sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                       s(:,2,isym)*xk_cryst(2) + &
                       s(:,3,isym)*xk_cryst(3) 
              xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

              !CALL start_clock ('exxen2_ngmloop')
              do ig = 1, ngm
                 q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                 q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                 q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                 q = q * tpiba
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                 do beta = 1, 3
                     fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                 enddo

                 if (x_gamma_extrapolation) then
                    on_double_grid = .true.
                    x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 endif

                 if (use_coulomb_vcut_ws) then
                    fac(ig) = vcut_get(vcut, q)
                    fac_stress(ig) = 0.d0   ! not implemented
                    if (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                 else if ( use_coulomb_vcut_spheric ) then
                    fac(ig) = vcut_spheric_get(vcut, q)
                    fac_stress(ig) = 0.d0   ! not implemented
                    if (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 

                 else if (qq > 1.d-8) then
                    if ( erfc_scrlen > 0 ) then
                       fac(ig)=e2*fpi/qq*(1-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                       fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                           (1.d0+qq/4.d0/erfc_scrlen**2)*exp(-qq/4.d0/erfc_scrlen**2) - 1.d0) * &
                           grid_factor
                    else
                       fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                       fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                    end if

                    if (gamma_only) fac(ig) = 2.d0 * fac(ig)
                    if (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                    if (on_double_grid) fac(ig) = 0.d0
                    if (on_double_grid) fac_stress(ig) = 0.d0

                 else
                    fac(ig)= -exxdiv ! or rather something else (see f.gygi)
                    fac_stress(ig) = 0.d0  ! or -exxdiv_stress (not yet implemented)
                    if ( use_yukawa .and. .not. x_gamma_extrapolation) then
                       fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                       fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                    endif
                    if (erfc_scrlen > 0.d0 .and. .not. x_gamma_extrapolation) then
                       fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)
                       fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                    endif
                 endif
              enddo
              !CALL stop_clock ('exxen2_ngmloop')

              if (gamma_only) then
                  half_nbnd = (nbnd + 1) / 2
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

                      ! loads the phi from file
                      do ji=1, nrxxs
                        tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                      enddo
                
                      ! calculate rho in real space
                      rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                      ! brings it to G-space
                      CALL fwfft ('Smooth', rhoc, dffts)
   
                      vc = 0.d0
                      do ig=1,ngm
                        vc(:,:) = vc(:,:) + fac(ig) * x1 * &
                                  abs( rhoc(nls(ig))+CONJG(rhoc(nlsm(ig))))**2 * &
                                  (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        vc(:,:) = vc(:,:) + fac(ig) * x2 * &
                                  abs( rhoc(nls(ig))-CONJG(rhoc(nlsm(ig))))**2 * &
                                  (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      enddo
                      vc = vc / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                  enddo
              else
                  do ibnd=1,nbnd !for each band of psi
                     if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle

                     ! loads the phi from file
                     do ji=1, nrxxs
                       tempphic(ji)=exxbuff(ji,ikq,ibnd)
                     enddo

                     ! calculate rho in real space
                     rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega

                     ! brings it to G-space
                     CALL fwfft ('Smooth', rhoc, dffts)

                     vc = 0.d0
                     do ig=1,ngm
                       vc(:,:) = vc(:,:) + rhoc(nls(ig))*CONJG(rhoc(nls(ig))) * &
                                 (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                     end do
                     vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                     exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                  end do
              endif ! gamma or k-points

          enddo ! iqi
      enddo ! jbnd
  enddo ! ikk

  deallocate (tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
!
! Was used for image parallelization
!  call mp_sum( exx_stress_, inter_image_comm )
!
  call mp_sum( exx_stress_, intra_pool_comm )
  call mp_sum( exx_stress_, inter_pool_comm )

  exx_stress = exx_stress_

  call stop_clock ('exx_stress')

  END FUNCTION exx_stress



end module exx
