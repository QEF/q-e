  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------
  SUBROUTINE wann_run
  !---------------------------------------------------------------------
  !
  !  This is the subroutine which controls the w90 run.  Primarily,        
  !  we get the phases to remove degeneracies in the wfs, and 
  !  call pw2wan90epw 
  !  
  !---------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout, ionode_id
  USE wvfct,          ONLY : nbnd
  USE ions_base,      ONLY : nat
  USE start_k,        ONLY : nk1, nk2, nk3
  USE pwcom,          ONLY : nkstot !, npk
  USE epwcom,         ONLY : xk_cryst           
  USE wannier,        ONLY : mp_grid, n_wannier, kpt_latt
  USE mp,             ONLY : mp_bcast
  USE mp_world,       ONLY : world_comm
  !
  implicit none
  !
  ! intend in variables
  real(kind=DP) :: zero_vect(3)
  real(kind=DP), ALLOCATABLE :: tau_wan(:,:)
  integer :: num_kpts
  logical :: gamma_only_wan
  !
  !
  CALL start_clock( 'WANNIER' )
  ! 
  ! For small systems we do not do gamma ONLY.  
  gamma_only_wan = .false.
  !
  zero_vect = 0.0
  !
  mp_grid(1) = nk1
  mp_grid(2) = nk2
  mp_grid(3) = nk3
  num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
  !
  IF (num_kpts .ne. nkstot ) call errore('wannierize','inconsistent nscf and elph k-grids',1) 
  IF (nbnd .lt. n_wannier )  call errore('wannierize','Must have as many or more bands than Wannier functions',1) 
  !ALLOCATE ( kpt_latt(3,npk), tau_wan(3,nat) )
  ALLOCATE ( kpt_latt(3,num_kpts), tau_wan(3,nat) )
  !
  WRITE(stdout, '(5x,a)') repeat("-",67)
  WRITE(stdout, '(a, i2,a,i2,a,i2,a)') "     Wannierization on ", nk1, " x ", nk2, " x ", nk3 , " electronic grid"
  WRITE(stdout, '(5x,a)') repeat("-",67)
  !
  kpt_latt = xk_cryst(:,1:num_kpts)
  CALL mp_bcast(kpt_latt, ionode_id, world_comm)
  !
  ! write the short input file for the wannier90 code
  !
  CALL write_winfil
  !
  ! run the wannier90 code to create MLWFs
  !
  CALL pw2wan90epw
  !
  ! project the Wannier functions onto energy space
  !
!  CALL proj_w90
  !
  WRITE(stdout, '(5x,a)') repeat("-",67)
  CALL print_clock( 'WANNIER' )
  WRITE(stdout, '(5x,a)') repeat("-",67)
  !
  end SUBROUTINE wann_run
  !------------------------------------------------------------
  !
  !------------------------------------------------------------
  SUBROUTINE write_winfil
  !------------------------------------------------------------
  !!
  !!
  !!  This SUBROUTINE write the prefix.win file which wannier90.x
  !!  needs to run.  Primarily it contains information about the 
  !!  windows USEd for the disentanglement, and the initial projections.
  !!  JN - 10/2008  projections now in elph.in file  
  !------------------------------------------------------------
  !
  USE io_files,    ONLY : prefix
  USE io_epw,      ONLY : iuwinfil
  USE io_global,   ONLY : meta_ionode
  USE epwcom,      ONLY : nbndsub, nwanxx, proj, iprint, dis_win_min, &
                          dis_win_max, dis_froz_min, dis_froz_max, num_iter, &
                          wdata 
  !
  implicit none
  !
  integer :: i
  !
  logical :: random
  !
  IF (meta_ionode) THEN
    !
    IF (nbndsub .gt. nwanxx) call errore('write_winfil',"Too many wannier bands",nbndsub)
    !
    OPEN (unit = iuwinfil, file = trim(prefix)//".win", form = 'formatted')
    !    
    !  more input and options for interfacing with w90 can/will be added later
    WRITE (iuwinfil,'(a)') "begin projections"
    !
    random = .true.
    DO i = 1, nbndsub+1
       IF (proj(i) .ne. ' ') THEN
          WRITE (iuwinfil,*) trim(proj(i))
          random = .false.
       ENDIF
    ENDDO
    !
    IF (random) WRITE(iuwinfil,*) 'random' 
    !
    WRITE (iuwinfil,'(a)') "end projections"
    !
    WRITE (iuwinfil,'("num_wann ",i3)') nbndsub
    WRITE (iuwinfil,'("iprint ",i3)') iprint
    !
    WRITE (iuwinfil, '("dis_win_min ", f9.3)')  dis_win_min
    WRITE (iuwinfil, '("dis_win_max ", f9.3)')  dis_win_max
    WRITE (iuwinfil, '("dis_froz_min ", f9.3)') dis_froz_min
    WRITE (iuwinfil, '("dis_froz_max ", f9.3)') dis_froz_max
    WRITE (iuwinfil, '("num_iter ", i7)')       num_iter
    !
    ! write any extra parameters to the prefix.win file
    DO i = 1, nwanxx
       IF (wdata(i) .ne. ' ') WRITE(iuwinfil,*) wdata(i)
    ENDDO
    !
    CLOSE (iuwinfil)
    !
  ENDIF
  !
  END SUBROUTINE write_winfil



!------------------------------------------------------------
  SUBROUTINE proj_w90
!------------------------------------------------------------
  !
  ! This SUBROUTINE computes the energy projections of
  ! the computed Wannier functions
  ! 07/2010  Needs work.  Right now this sub is nearly worthless  
  !------------------------------------------------------------
  !
  USE io_files,    ONLY : prefix 
  USE io_global,   ONLY : stdout
  USE io_epw,      ONLY : iuprojfil
  USE mp_global,   ONLY : inter_pool_comm
  USE io_global,   ONLY : stdout, meta_ionode
  USE mp,          ONLY : mp_sum
  USE epwcom,      ONLY : DP, dis_win_max, dis_win_min
  USE constants_epw, ONLY : ryd2ev 
  USE wannier,     ONLY : n_wannier
  USE wvfct,       ONLY : nbnd, et
  USE klist,       ONLY : nks, nkstot
  !
  implicit none
  !
  integer :: ik, ibnd, ne, ie, iwann
  real(kind=DP)    :: dE, sigma, argv, en, xxq(3)
  real(kind=DP), ALLOCATABLE    ::  proj_wf(:,:)
  complex(kind=DP), ALLOCATABLE ::  cu(:,:,:), cuq(:,:,:)
  !
  logical :: lwin( nbnd, nks ), lwinq( nbnd, nks )
  ! FG: introduced after extensive compiler tests
  !
  WRITE(stdout,'(5x,"Computing energy projections")')
  ! dummy var
  xxq = 0.d0
  !
  ! tmp value
  dE = 0.05d0
  sigma = 2.d0 * dE
  !
  ! maxvalue = dis_win_max + 1
  ! minvalue = dis_win_min - 1
  ne = int( (dis_win_max - dis_win_min + 1) / dE )
  IF (ne .lt. 1)  CALL errore('proj_wan','Problem with disentanglement window',1)
  !
  ALLOCATE (proj_wf(n_wannier, ne+1))
  proj_wf = 0.d0
  !
  ALLOCATE (cu (nbnd, n_wannier, nks) )
  ALLOCATE (cuq(nbnd, n_wannier, nks) )
  !
  CALL loadumat(nbnd, n_wannier, nks, nkstot, xxq, cu, cuq, lwin, lwinq) 
  ! FG: introduced after ifort checks

  !
  !
  DO iwann = 1, n_wannier
     !
     DO ie = 1, ne
        en = dble(ie)/dble(ne) * (dis_win_max - dis_win_min + 1.0) + dis_win_min
        !
        DO ik = 1, nks
           DO ibnd = 1, nbnd
              !
              argv = ( et(ibnd,ik)*ryd2ev - en ) **2/ (2 * sigma **2)
              proj_wf(iwann, ie ) = proj_wf(iwann, ie ) + exp(-argv) * real (  cu(ibnd, iwann,ik ) * conjg( cu(ibnd, iwann,ik) ))
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! sum the contributions from all k-points
  CALL mp_sum(proj_wf, inter_pool_comm)
  !
  IF (meta_ionode) THEN
    !
    OPEN (unit = iuprojfil, file = trim(prefix)//".projw90", form = 'formatted')
    !
    WRITE(iuprojfil, '(5x,"Wannier energy projections")')
    !
    DO ie = 1, ne
       en =  dble(ie)/dble(ne) * (dis_win_max - dis_win_min + 1) + dis_win_min
       WRITE(iuprojfil, '(f9.3, 25f8.4)' )  en , proj_wf(:, ie)
    ENDDO
    !
    CLOSE (iuprojfil)
  ENDIF
  !
  IF ( ALLOCATED(proj_wf)) DEALLOCATE(proj_wf)
  IF ( ALLOCATED(cu))      DEALLOCATE(cu)
  IF ( ALLOCATED(cuq))     DEALLOCATE(cuq)
  !
!------------------------------------------------------------
  END  SUBROUTINE proj_w90
!------------------------------------------------------------
