!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!----------------------------------------------------------------------
MODULE exphonon
!----------------------------------------------------------------------
  !!
  !! This module contains various printing routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE prepare_ex_ph_g()
  !-----------------------------------------------------------------------
    USE input,         ONLY : nqc1, nqc2, nqc3, nkc1, nkc2, nkc3, nbndv_explrn, nbndc_explrn, negnv_explrn
    USE global_var,    ONLY : nktotf, nqtotf, nkf, eigval_ex, wf_temp, wf, A_cv_all, epf17      
    USE modes,         ONLY : nmodes
    USE pwcom,         ONLY : nks
    USE io_global,     ONLY : stdout, ionode_id, meta_ionode_id
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : npool, my_pool_id, world_comm
    USE ep_constants,  ONLY : czero
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    nktotf=nkc1*nkc2*nkc3       
    nqtotf=nqc1*nqc2*nqc3
    nkf=nks
    !
    ! WRITE(stdout,'(/5x,a,I10)')'nktotf: ', nktotf
    ! WRITE(stdout,'(/5x,a,I10)')'nqtotf: ', nqtotf
    ! WRITE(stdout,'(/5x,a,I10)')'nkf : ', nkf 
    !
    ALLOCATE(wf(nmodes,nqc1*nqc2*nqc3), STAT=ierr)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error allocating wf', 1) 
    !
    ALLOCATE(eigval_ex(negnv_explrn, nqtotf), STAT=ierr)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error allocating eigval_ex', 1) 
    eigval_ex = 0.d0
    !
    ALLOCATE(epf17(nbndc_explrn+nbndv_explrn,nbndc_explrn+nbndv_explrn,nks,nmodes),STAT=ierr)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error allocating epf17(nbndfst,nbndfst,nktotf,nmodes)', 1)
    epf17 = czero
    !
    IF(my_pool_id == meta_ionode_id) THEN
      ALLOCATE(A_cv_all(nbndv_explrn,nbndc_explrn,nktotf,negnv_explrn,nktotf),STAT=ierr)
      A_cv_all = 0.d0
      CALL read_Acv_all()
    ENDIF
    !
    IF( my_pool_id .EQ. meta_ionode_id ) THEN
      CALL system('mkdir -p G_full_epmatq') 
    ENDIF
    !
    CALL mp_barrier(world_comm)
  !-----------------------------------------------------------------------
  END SUBROUTINE prepare_ex_ph_g
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE release_ex_ph_g()
  !-----------------------------------------------------------------------
    USE input,         ONLY : nqc1, nqc2, nqc3, nkc1, nkc2, nkc3, nbndv_explrn, nbndc_explrn, negnv_explrn
    USE global_var,    ONLY : nktotf, nqtotf, nkf, eigval_ex, wf_temp, wf, A_cv_all, epf17      
    USE modes,         ONLY : nmodes
    USE pwcom,         ONLY : nks
    USE io_global,     ONLY : stdout, ionode_id
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : npool, my_pool_id, world_comm
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! IF(ALLOCATED(wf)) DEALLOCATE(wf, STAT=ierr)
    DEALLOCATE(wf, STAT=ierr)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error deallocating wf', 1) 
    !
    ! IF(ALLOCATED(eigval_ex)) DEALLOCATE(eigval_ex,STAT=ierr)
    DEALLOCATE(eigval_ex,STAT=ierr)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error deallocating eigval_ex', 1) 
    !
    ! IF(ALLOCATED(A_cv_all)) DEALLOCATE(A_cv_all)
    IF(my_pool_id == ionode_id) THEN
      DEALLOCATE(A_cv_all)
      IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error deallocating A_cv_all', 1) 
    ENDIF
    !
    ! IF(ALLOCATED(epf17)) DEALLOCATE(epf17)
    DEALLOCATE(epf17)
    IF(ierr /=0) call errore('elphon_shuffle_wrap', 'Error deallocating epf17', 1) 
  !-----------------------------------------------------------------------
  END SUBROUTINE release_ex_ph_g
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE read_Acv_all()
  !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE input,         ONLY : nbndc_explrn, nbndv_explrn, negnv_explrn, nqc1, nqc2, nqc3                 
    USE global_var,    ONLY : A_cv_all                                               
    USE ep_constants,  ONLY : ci                                        
    !
    IMPLICIT NONE
    !
    INTEGER                    :: iq
    !! Index for exciton momentum
    INTEGER                    :: ik
    !! Index for crystal momentm
    INTEGER                    :: ibnd
    !! Index for electron band index
    INTEGER                    :: nktotf_f 
    !! Number of exciton momentum in the grid.
    INTEGER                    :: ikminusiexq
    !! Index for ik - iq
    INTEGER                    :: dims(7)
    !! Auxiliary array for reading BGW output
    REAL(KIND=8), ALLOCATABLE  :: Aq_temp(:,:,:,:,:,:,:)
    !! Auxiliary array for reading BGW output
    CHARACTER(LEN = 256)       :: filename
    !! Filenames
    CHARACTER(LEN = 256)       :: grpname
    !! Group names. For HDF5
    CHARACTER(LEN = 256)       :: dsetname
    !! Data set names. For HDF5
    CHARACTER(LEN = 256)       :: fileexq
    !! Temporary name to convert integer to character
    !
    nktotf_f = nqc1 * nqc2 * nqc3
    !
    DO iq = 1, nktotf_f
      WRITE(fileexq,"(I10)") iq - 1 !ZD: 0-index. Consistent with ST convention.
      !
      filename = './eigv/q_' // trim(adjustl(fileexq)) // '/eigenvectors.h5'
      grpname="exciton_data"
      dsetname="eigenvectors"
      !    
      dims=(/2,1,nbndv_explrn,nbndc_explrn,nktotf_f,negnv_explrn,1/)
      !
      ALLOCATE(Aq_temp(2,1,nbndv_explrn,nbndc_explrn,nktotf_f,negnv_explrn,1))
      Aq_temp=0.d0
      !
#if defined(__HDF5)
      CALL READ_HDF5_SLAB(filename,grpname,dsetname,Aq_temp(:,:,:,:,:,:,:),7,negnv_explrn,dims)
#endif
      ! Broadcast the A_cvq over pools.
      !!ZD: Reorder both k and vband and this step
      !
      DO ik=1,nktotf_f
        CALL indexing_q2(iq, ik, nqc1, nqc2, nqc3, -1, ikminusiexq)
        DO ibnd=1,nbndv_explrn
          A_cv_all(ibnd,:,ikminusiexq,:,iq)= &
          Aq_temp(1,1, nbndv_explrn-ibnd + 1,:,ik,:,1) + ci*Aq_temp(2,1,nbndv_explrn-ibnd + 1,:,ik,:,1)
        ENDDO
      ENDDO
      !
      DEALLOCATE(Aq_temp)
      !
    ENDDO    
  !-----------------------------------------------------------------------
  END SUBROUTINE read_Acv_all
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE calc_G_epmat_subbnd(iq, ph_eigvec, iq_last)
  !-----------------------------------------------------------------------
    !! ZD: Adapted from STs qdabs module
    !! The routine that reads the BSE eigenvector a_sQ+q and a_s'Q and prepare the calculation of ex-ph g.
    !! Remember the ordering of vband and the convention for Q.
    !! (In BGW, specifying Q in kernel.inp and absorption.inp actually gives a exciton with -Q.)
    !! We consier every Q and q in G_ss'\nu(Q,q). So, full_G_expmatq(iQ, iex, iex', inu) for each iq.
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id, meta_ionode, meta_ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fsthick, degaussw,                           &
                              start_mesh,mode_res,nbndv_explrn,nbndc_explrn, negnv_explrn,              &
                              nqc1, nqc2, nqc3, nbndsub
    USE global_var,    ONLY : etf, ibndmin, nkf, epf17, wkf, wf, eigval_ex,        &
                              gtemp,nkqtotf,nkqf, totf, nktotf, xkf, xqf, A_cvq,   &
                              G_epmatq, A_cvqplusq, G_full_epmatq, A_cv_all, nkf                       
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci,  &
                              eps6, czero, eps8
    USE mp,            ONLY : mp_barrier, mp_sum, mp_gather, mp_bcast
    USE mp_global,     ONLY : npool,my_pool_id, world_comm
    USE cell_base,     ONLY : omega,bg,at
    USE mp_world,      ONLY : mpime
    USE bzgrid,        ONLY : kpmq_map
    USE low_lvl,       ONLY : init_random_seed  
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)             :: iq ! 
    !! iIdices for phonon q
    INTEGER, INTENT(in)             :: iq_last 
    !! The index of the last iq in symmetry shuffled kgrid
    COMPLEX(KIND = DP), INTENT(in)  :: ph_eigvec(nmodes, nmodes)
    !! The rotated eigenvectors, for the current q in the star
    INTEGER                         :: iexq, iexqplusq, ikminusiexq 
    !! Indices for exciton Q, Q+q, k-Q, respectively
    INTEGER                         :: dims(7)
    !! Auxiliary array for reading BGW output
    CHARACTER(LEN = 256)            :: filename
    !! Filenames
    CHARACTER(LEN = 256)            :: grpname
    !! Group names. For HDF5
    CHARACTER(LEN = 256)            :: dsetname
    !! Data set names. For HDF5
    CHARACTER(LEN = 20)             :: tp
    !! Temporary string for integer conversion
    CHARACTER(LEN = 20)             :: fileexq
    !! BGW eigenvector file with exciton momentum Q
    INTEGER                         :: imode, imode1
    !! Counter on phonon modes
    INTEGER                         :: ierr
    !! Error status
    INTEGER                         :: ibnd, ik, sbnd, s1bnd
    !! Indices for electronic bands, crystal momentum, exciton bands
    REAL(KIND=8),ALLOCATABLE :: Aq(:,:,:,:,:,:,:), temp_eigval(:)
    !! Auxiliary array for reading BGW output
    REAL :: start, finish, middle
    !! Variables for timing
    !
    CALL CPU_TIME(start)
    WRITE(stdout,'(/5x,a)')'Computing exciton-phonon coupling.'     
    !! ZD: Computing the full exciton-phonon matrix elements. 
    !! Aiming the parallelization over phonon q, not exciton momentum Q, due to the limitation of EPW
    ALLOCATE(wkf(5),STAT=ierr)
    IF(ierr/=0) CALL errore('qdabs', 'Error allocating wkf_write(tot)', 1)
    wkf(:)=REAL(1.0/(nktotf))
    !
    ! IF(ALLOCATED(G_full_epmatq)) DEALLOCATE(G_full_epmatq)
    ALLOCATE(G_full_epmatq(negnv_explrn,negnv_explrn,nktotf, nmodes),STAT=ierr)
    IF(ierr /=0) CALL errore('qdabs_BGW','error allocating G_full_epmatq',1)
    G_full_epmatq=0.d0
    !
    ALLOCATE(temp_eigval(nbndv_explrn*nbndc_explrn*nktotf), STAT=ierr)
    !
    DO iexq = 1,nktotf
      WRITE(fileexq,"(I10)") iexq - 1 !ZD: 0-index. Consistent with ST convention.
      ! CALL indexing_q(iq, iexq, nqx, 1, iexqplusq)
      CALL indexing_q2(iq, iexq, nqc1, nqc2, nqc3, 1, iexqplusq)
      !
      !! Read and distribute BSE eigvectors for Q
      !
      filename = './eigv/q_' // trim(adjustl(fileexq)) // '/eigenvectors.h5'
      grpname="exciton_data"
      dsetname="eigenvalues"
      !
      IF(my_pool_id==ionode_id .AND. iq .EQ. 1) THEN
        temp_eigval = 0.d0             
#if defined(__HDF5)
        CALL READ_HDF5_ENE(filename,grpname,dsetname,temp_eigval,1)
#endif
        eigval_ex(:, iexq) = temp_eigval(1:negnv_explrn)
      ENDIF
      !
      ALLOCATE(A_cvq(nbndv_explrn,nbndc_explrn,nktotf,negnv_explrn),STAT=ierr)
      A_cvq = 0.d0
      !
      ALLOCATE(A_cvqplusq(nbndv_explrn,nbndc_explrn,nktotf,negnv_explrn),STAT=ierr)
      A_cvqplusq = 0.d0
      !
      IF(my_pool_id==ionode_id) THEN
        A_cvq(:,:,:,:) = A_cv_all(:,:,:,:,iexq)
        A_cvqplusq(:,:,:,:) = A_cv_all(:,:,:,:,iexqplusq)
      ENDIF
      !
      CALL mp_barrier(world_comm)
      CALL mp_bcast(A_cvq, meta_ionode_id, world_comm)
      CALL mp_bcast(A_cvqplusq, meta_ionode_id, world_comm)
      CALL mp_barrier(world_comm)       
      !
      ! ZD: Start the construction of the ex-ph matrix: G_ss'nu(Q,q) for fixed q
      ! WRITE(stdout,'(/5x,a, I10)')'Start computing exciton-phonon matrix.', iexq
      !
      CALL calculate_epmat_oneexq(iexq, iexqplusq)
      !
      CALL mp_barrier(world_comm)
      !
      DEALLOCATE(A_cvq)
      DEALLOCATE(A_cvqplusq)
      CALL mp_barrier(world_comm)
    ENDDO
    !
    IF (iq == 1) CALL mp_bcast(eigval_ex, meta_ionode_id, world_comm)
    !
    ! Formatted writting. The quantities can be used to some analyses 
    ! Write phonon frequencies
    IF( iq_last==nktotf .AND. my_pool_id .EQ. ionode_id ) THEN 
      filename = './G_full_epmatq/phband.fmt'
      OPEN(UNIT=iuexphg, FILE = filename, ACTION = 'write', IOSTAT = ierr)
      DO iexq=1,nktotf
        DO imode=1,nmodes-1
          WRITE(iuexphg, '(e20.10)',advance='no') wf(imode, iexq)
        ENDDO
        WRITE(iuexphg, '(e20.10)') wf(nmodes, iexq)
      ENDDO
      CLOSE(iuexphg)
    ENDIF
    !
    ! Write exciton energies
    IF( iq==1 .AND. my_pool_id .EQ. ionode_id ) THEN 
      filename = './G_full_epmatq/exband.fmt'
      OPEN(UNIT=iuexphg, FILE = filename, ACTION = 'write', IOSTAT = ierr)
      DO iexq=1,nktotf
        DO imode=1,negnv_explrn-1
          WRITE(iuexphg, '(e20.10)',advance='no') eigval_ex(imode, iexq)
        ENDDO
        WRITE(iuexphg, '(e20.10)') eigval_ex(negnv_explrn, iexq)
      ENDDO
      CLOSE(iuexphg)
    ENDIF
    ! Write exciton-phonon matrix elements 
    WRITE(tp,"(I10)") iq-1  
    IF(my_pool_id .EQ. ionode_id) THEN  
      filename = './G_full_epmatq/G_full_epmatq_' // trim(adjustl(tp)) // '.dat'
      OPEN(UNIT=iuexphg, FILE = filename, ACTION = 'write', IOSTAT = ierr)
      DO iexq=1,nktotf
        ! WRITE(iuexphg, '(a, 10i)') 'iQ = ', iexq
        DO sbnd=1,negnv_explrn
          DO s1bnd=1,negnv_explrn
            DO imode=1,nmodes
              WRITE(iuexphg, '(4i9, f20.8, 2e20.10)') sbnd, s1bnd, iexq, imode, ryd2ev*wf(imode, iq), &
              ryd2ev*REAL(G_full_epmatq(sbnd, s1bnd, iexq, imode)), &
              ryd2ev*AIMAG(G_full_epmatq(sbnd, s1bnd, iexq, imode))                
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iuexphg)  
    ENDIF 
    !
    ! Write phonon eigenvectors 
    IF(my_pool_id .EQ. ionode_id) THEN  
      filename = './G_full_epmatq/ph_eigvec_' // trim(adjustl(tp)) // '.dat'
      OPEN(UNIT=iuexphg, FILE = filename, ACTION = 'write', IOSTAT = ierr)
      !
      DO imode=1,nmodes
        DO imode1=1,nmodes
          WRITE(iuexphg, '(2e20.10)') REAL(ph_eigvec(imode1, imode)), AIMAG(ph_eigvec(imode1, imode))               
        ENDDO
      ENDDO
      !
      CLOSE(iuexphg)  
    ENDIF  
    !
    DEALLOCATE(temp_eigval)
    DEALLOCATE(G_full_epmatq)
    DEALLOCATE(wkf,STAT=ierr)
    !
    CALL mp_barrier(world_comm)
    !
    CALL CPU_TIME(finish)
    WRITE(stdout,'(5x,a, I10, a/)') 'Finish computing exciton-phonon coupling. Costing ', &
                                      CEILING(finish - start), ' seconds'
    !                                      
  !-----------------------------------------------------------------------
  END SUBROUTINE calc_G_epmat_subbnd
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_epmat_oneexq(iexq,iexqplusq)
  !----------------------------------------------------------------------- 
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, meta_ionode_id
    USE io_var,        ONLY : iuindabs, iuexphg
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fermi_energy,                            &
                                degaussw, fsthick,filkf,nbndv_explrn,nbndc_explrn,negnv_explrn,     &
                                only_c_explrn, only_v_explrn, nqc1, nqc2, nqc3, only_pos_modes_explrn     
    USE global_var,    ONLY : etf, ibndmin, nkf, wf,                           &
                              nbndfst, nktotf,nkqtotf,epf17,nkqf,E_mesh,       &
                              xkf,xqf,wkf,gtemp,totf_pool,totcv_pool,tot_pool, &
                              index_buf_pool,H_quad,stop_qdabs,Energy,maxdim,  & 
                              G_epmatq,A_cvq, A_cvqplusq, G_full_epmatq,nbndep
    USE mp_global,     ONLY : my_pool_id,npool,ionode_id, world_comm
    USE mp,            ONLY : mp_barrier,mp_bcast,mp_sum
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, eps40, czero, eps8
    USE parallelism,   ONLY : fkbounds
    !
    INTEGER, INTENT(IN)                :: iexq
    !! Index for exciton momentum Q
    INTEGER                            :: lower_bnd, upper_bnd
    !! Lower and upper bound for k parallelization
    INTEGER, OPTIONAL                  :: iexqplusq 
    !! Index for Q+q. 
    INTEGER                            ::  imode
    !! Index for phonon mode
    INTEGER                            ::  sbnd, s1bnd
    !! Index for exciton band
    INTEGER                            ::  ik
    !! Index for crystal momentum k
    INTEGER                            ::  cbnd, vbnd, c1bnd, v1bnd
    !! Index for valence band and conduction band
    INTEGER                            ::  ierr
    !! Error status
    INTEGER                            ::  iphononq
    !! Index for phonon momentum q
    INTEGER                            ::  ikplusiphq, ikminusiexq
    !! Index for k+q and k-Q
    COMPLEX(KIND=DP), ALLOCATABLE      :: temp_G(:,:,:) 
    !! Temporary array for storing ex-ph matrix elements
    REAL(KIND=DP)                      :: abs_wf 
    !! Abs value of the phonon frequency
    !
    ! Check bound
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    CALL indexing_q2(iexq, iexqplusq, nqc1, nqc2, nqc3, -1, iphononq) 
    !
    ! WRITE(stdout,'(/5x,a, I10)')'inside epmatq for phonon q: ', iphononq
    ALLOCATE(temp_G(negnv_explrn, negnv_explrn, nmodes))
    temp_G = 0.d0
    !
    DO ik=1,nkf ! kpoint in pool
      !!ZD: Assume A is in our convention
      CALL indexing_q2(iexq, ik + lower_bnd - 1, nqc1, nqc2, nqc3, -1, ikminusiexq)
      ! k-Q
      CALL indexing_q2(ik + lower_bnd - 1, iphononq, nqc1, nqc2, nqc3, 1, ikplusiphq)
      ! k+q
      ! 
      DO sbnd=1,negnv_explrn  !! Exciton band s 
        DO s1bnd=1,negnv_explrn  !! Exciton band s'
          DO imode=1,nmodes  ! phonon mode imode    
            !       
            IF(only_pos_modes_explrn .EQV. .TRUE.) THEN
              abs_wf = wf(imode, iphononq)
            ELSE
              abs_wf = ABS(wf(imode, iphononq))
            ENDIF
            !
            DO cbnd=1,nbndc_explrn  ! conduction band
              DO vbnd=1,nbndv_explrn ! valence band
                !!Assume A is now in our convention
                IF((only_c_explrn .EQV. .TRUE.) .AND. (only_v_explrn .EQV. .FALSE.)) THEN
                  DO c1bnd=1,nbndc_explrn ! intermediate band conduction
                    !IF(ABS(wf(imode, iphononq)) .LE. eps6) THEN
                    IF(abs_wf .LE. eps6) THEN
                      temp_G(sbnd, s1bnd, imode)=temp_G(sbnd, s1bnd, imode)
                    ELSE
                      temp_G(sbnd, s1bnd, imode)= &
                                              temp_G(sbnd, s1bnd, imode) + &
                                              CONJG(A_cvqplusq(vbnd, cbnd, ikminusiexq, sbnd)) * &
                                              epf17(cbnd+nbndv_explrn, c1bnd+nbndv_explrn, ik, imode)* & 
                                              A_cvq(vbnd, c1bnd, ikminusiexq, s1bnd)* &
                                              DSQRT(1.0/( 2.0 * ABS(wf(imode, iphononq)) ))
                      ! ZD: ABS(wf) for now. But should it really be ABS(wf) or i*(ABS(wf))?   
                    ENDIF
                  ENDDO
                ELSEIF((only_c_explrn .EQV. .FALSE.) .AND. (only_v_explrn .EQV. .TRUE.)) THEN
                  DO v1bnd=1,nbndv_explrn ! intermediate band valence
                    !IF(ABS(wf(imode, iphononq)) .LE. eps6) THEN
                    IF(abs_wf <= eps6) THEN
                      temp_G(sbnd, s1bnd, imode) = temp_G(sbnd, s1bnd, imode)
                    ELSE
                      temp_G(sbnd, s1bnd, imode)= &
                                              temp_G(sbnd, s1bnd, imode)- &
                                              CONJG(A_cvqplusq(vbnd, cbnd, ik+lower_bnd-1, sbnd)) * &
                                              epf17(v1bnd, vbnd, ik, imode)* & 
                                              A_cvq(v1bnd, cbnd, ikplusiphq, s1bnd) * & 
                                              DSQRT(1.0/(2.0 * ABS(wf(imode, iphononq)) ))                       
                    ENDIF
                  ENDDO
                ELSE
                  DO c1bnd=1,nbndc_explrn ! intermediate band conduction
                    !IF(ABS(wf(imode, iphononq)) .LE. eps6) THEN
                    IF(abs_wf .LE. eps6) THEN
                      temp_G(sbnd,s1bnd,imode)=temp_G(sbnd,s1bnd,imode)
                    ELSE
                      temp_G(sbnd,s1bnd,imode)= &
                                              temp_G(sbnd,s1bnd,imode) +   &
                                              CONJG(A_cvqplusq(vbnd,cbnd,ikminusiexq,sbnd))*         &
                                              epf17(cbnd+nbndv_explrn,c1bnd+nbndv_explrn,ik,imode)*   & 
                                              A_cvq(vbnd,c1bnd,ikminusiexq,s1bnd)*                &
                                              DSQRT(1.0/(2.0 * ABS(wf(imode,iphononq)) ))
                    ENDIF
                  ENDDO
                  !
                  DO v1bnd=1,nbndv_explrn ! intermediate band valence
                    !IF(ABS(wf(imode, iphononq)) .LE. eps6) THEN
                    IF(abs_wf .LE. eps6) THEN
                      temp_G(sbnd,s1bnd,imode) = temp_G(sbnd,s1bnd,imode)
                    ELSE
                      temp_G(sbnd,s1bnd,imode)= &
                                              temp_G(sbnd,s1bnd,imode)- &
                                              CONJG(A_cvqplusq(vbnd,cbnd,ik+lower_bnd-1,sbnd))*      &
                                              epf17(v1bnd,vbnd,ik,imode)*                & 
                                              A_cvq(v1bnd,cbnd,ikplusiphq,s1bnd)*             & 
                                              DSQRT(1.0/(2.0 * ABS(wf(imode,iphononq)) ))
                    ENDIF
                  ENDDO                                        
                ENDIF  
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CALL mp_barrier(world_comm)
    CALL mp_sum(temp_G, world_comm)
    G_full_epmatq(:,:,iexq,:) = temp_G(:,:,:)
    DEALLOCATE(temp_G)
    CALL mp_barrier(world_comm)       
  !-----------------------------------------------------------------------
  END SUBROUTINE calculate_epmat_oneexq
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE indexing_q2(iq1, iq2, nqc1, nqc2, nqc3, mode, idq)
  !-----------------------------------------------------------------------  
    !! This finds the index of global iq index of iq2 + mode * iq1, 
    !! NOT assuming nqc1 = nqc2 = nqc3
    IMPLICIT NONE
    !
    INTEGER, INTENT(in)  :: iq1, iq2
    !! Index of the first q and second q
    INTEGER, INTENT(in)  :: nqc1, nqc2, nqc3
    !! Dimensions of the uniform q grid
    INTEGER, INTENT(in)  :: mode
    !! mode=1 means iq2+iq1, mode=-1 means iq2-iq1
    INTEGER, INTENT(out) :: idq
    !! Index of iq2+iq1 or iq2-iq1
    INTEGER              :: qx1, qy1, qz1
    !! The x, y, z index of iq1
    INTEGER              :: qx2, qy2, qz2
    !! The x, y, z index of iq2
    INTEGER              :: dqx, dqy, dqz
    !! The x, y, z index of idq
    INTEGER              :: shift_x, shift_y, shift_z
    !! The shift of x, y, z if (iq1 mode iq2) is beyond the 1BZ
    !
    qz1 = mod(iq1 - 1, nqc3)
    qy1 = mod( (iq1 - 1 - qz1) / nqc3, nqc2 )
    qx1 = (iq1 - 1) / (nqc3 * nqc2)
    !
    qz2 = mod(iq2 - 1, nqc3)
    qy2 = mod( (iq2 - 1 - qz2) / nqc3, nqc2 )
    qx2 = (iq2 - 1) / (nqc3 * nqc2)
    !
    if(mode == -1) then
        dqx = qx2 - qx1
        dqy = qy2 - qy1
        dqz = qz2 - qz1
    else
        dqx = qx2 + qx1
        dqy = qy2 + qy1
        dqz = qz2 + qz1
    end if
    !
    shift_x = dqx / nqc1
    if(dqx < 0) shift_x = shift_x - 1
    !
    shift_y = dqy / nqc2
    if(dqy < 0) shift_y = shift_y - 1
    !
    shift_z = dqz / nqc3
    if(dqz < 0) shift_z = shift_z - 1
    !
    dqx = dqx - shift_x * nqc1
    dqy = dqy - shift_y * nqc2
    dqz = dqz - shift_z * nqc3
    !
    idq = dqx * nqc2 * nqc3 + dqy * nqc3 + dqz + 1
  !-----------------------------------------------------------------------
  END SUBROUTINE indexing_q2
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE find_iq(at, nqc1, nqc2, nqc3, iq)
  !-----------------------------------------------------------------------  
    !! This finds the global index of a given q point in crystal coordinate
    USE kinds,      ONLY  : DP
    IMPLICIT NONE
    !
    INTEGER              :: qx, qy, qz
    !! The x, y, z index of iq corresponding to iq
    INTEGER              :: shift_x, shift_y, shift_z
    !! The shift of x, y, z if iq is beyond the 1BZ
    REAL(DP), INTENT(in) :: at(3)
    !! The crystal coordinate of the input q point
    INTEGER, INTENT(in)  :: nqc1, nqc2, nqc3
    !! Dimensions of the uniform q grid
    INTEGER, INTENT(out) :: iq
    !! The index of q corresponding to input q coordinate
    !
    qx = NINT(at(1) * nqc1)
    qy = NINT(at(2) * nqc2)
    qz = NINT(at(3) * nqc3)
    !
    shift_x = qx / nqc1
    IF(qx < 0) shift_x = shift_x - 1
    !
    shift_y = qy / nqc2
    IF(qy < 0) shift_y = shift_y - 1
    !
    shift_z = qz / nqc3
    IF(qz < 0) shift_z = shift_z - 1
    !
    qx = qx - shift_x * nqc1
    qy = qy - shift_y * nqc2
    qz = qz - shift_z * nqc3
    !
    iq = qx * nqc2 * nqc3 + qy * nqc3 + qz + 1
  !-----------------------------------------------------------------------
  END SUBROUTINE find_iq  
  !-----------------------------------------------------------------------
  !
#if defined(__HDF5)
  !-----------------------------------------------------------------------
  SUBROUTINE READ_HDF5(filename,grpname,dsetname,data_out,N)
  !-----------------------------------------------------------------------
    USE hdf5
    !
    CHARACTER(LEN=100),intent(in) :: filename     
    !! File name
    CHARACTER(LEN=100),intent(in) :: grpname      
    !! Dataset name
    CHARACTER(LEN=100),intent(in) :: dsetname     
    !! Dataset name
    REAL(kind=8),intent(inout)    :: data_out(:,:,:,:,:,:,:)
    !! Output data
    INTEGER, intent(in)           :: N
    !! Data dimension
    INTEGER(HID_T)                :: file_id       
    !! File identifier
    INTEGER(HID_T)                :: gset_id       
    !! Data
    INTEGER(HID_T)                :: dset_id       
    !! Data
    INTEGER(HSIZE_T), ALLOCATABLE :: data_dims(:)
    !! Data dimension
    INTEGER                       :: error 
    !! Error status
    ALLOCATE(data_dims(N))
    !
    CALL h5open_f(error)
    IF(error /= 0) WRITE(*, '(a)') 'error1'
    !
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    IF(error /= 0) WRITE(*, '(a)') 'error2'  
    !
    CALL h5gopen_f(file_id, grpname, gset_id, error)
    IF(error /= 0) WRITE(*, '(a)') 'error3'
    !
    CALL h5dopen_f(gset_id, dsetname, dset_id, error)
    IF(error /= 0) WRITE(*, '(a)') 'error4'  
    !
    data_out=0.d0
    !
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
    DEALLOCATE(data_dims)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(gset_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
  !-----------------------------------------------------------------------
  END SUBROUTINE READ_HDF5
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE READ_HDF5_ENE(filename,grpname,dsetname,data_out,N)
  !-----------------------------------------------------------------------  
    USE hdf5
    !
    CHARACTER(LEN=100),intent(in) :: filename     
    !! File name
    CHARACTER(LEN=100),intent(in) :: grpname      
    !! Dataset name
    CHARACTER(LEN=100),intent(in) :: dsetname     
    !! Dataset name
    REAL(kind=8),intent(inout)    :: data_out(:)
    !! Output data
    INTEGER, intent(in)           :: N
    !! Dimensions
    INTEGER(HID_T)                :: file_id       
    !! File identifier
    INTEGER(HID_T)                :: gset_id       
    !! Data
    INTEGER(HID_T)                :: dset_id       
    !! Data
    INTEGER(HSIZE_T), ALLOCATABLE :: data_dims(:)
    !! Data dimension
    INTEGER                       :: error 
    !! Error status
    ALLOCATE(data_dims(N))
    !
    CALL h5open_f(error)
    if(error /= 0) write(*, '(a)') 'error1'
    !
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    if(error /= 0) write(*, '(a)') 'error2'  
    !
    CALL h5gopen_f(file_id, grpname, gset_id, error)
    if(error /= 0) write(*, '(a)') 'error3'
    !
    CALL h5dopen_f(gset_id, dsetname, dset_id, error)
    if(error /= 0) write(*, '(a)') 'error4'  
    !
    data_out=0.d0
    !
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, data_dims, error)
    !
    DEALLOCATE(data_dims)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(gset_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
  !-----------------------------------------------------------------------
  END SUBROUTINE READ_HDF5_ENE
  !-----------------------------------------------------------------------
  SUBROUTINE READ_HDF5_SLAB(filename,grpname,dsetname,data_out, N, negnv_explrn, in_dims)
  !-----------------------------------------------------------------------
    USE hdf5
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=100),intent(in)   :: filename     
    !! File name
    CHARACTER(LEN=100),intent(in)   :: grpname      
    !! Dataset name
    CHARACTER(LEN=100),intent(in)   :: dsetname     
    !! Dataset name
    REAL(kind=8),intent(inout)      :: data_out(:,:,:,:,:,:,:)
    !! Output data
    INTEGER, intent(in)             :: N, negnv_explrn
    !! Data dimension and number of BSE eigenvectors
    INTEGER, intent(in)             :: in_dims(7)
    !! Data dimension
    INTEGER(HID_T)                  :: file_id       
    !! File identifier
    INTEGER(HID_T)                  :: gset_id       
    !! Data
    INTEGER(HID_T)                  :: dset_id       
    !! Data
    INTEGER(HID_T)                  :: dataspace     
    !! Dataspace identifier 
    INTEGER(HID_T)                  :: memspace      
    !! memspace identifier 
    INTEGER(HSIZE_T), ALLOCATABLE   :: data_dims(:)
    !! Data dimension
    INTEGER(HSIZE_T), ALLOCATABLE   :: sub_data_dims(:)
    !! Data dimension of subset of data
    INTEGER(HSIZE_T), ALLOCATABLE   :: count(:)   
    !! Size of hyperslab
    INTEGER(HSIZE_T), ALLOCATABLE   :: offset(:) 
    !! Hyperslab offset
    INTEGER(HSIZE_T), ALLOCATABLE   :: stride(:) 
    !! Hyperslab stride 
    INTEGER(HSIZE_T), ALLOCATABLE   :: block(:)  
    !! Hyperslab block size 
    INTEGER                         :: i
    !! Counter
    INTEGER                         :: error
    !! Error status
    !
    ALLOCATE(data_dims(N))
    ALLOCATE(sub_data_dims(N))
    ALLOCATE(count(N))
    ALLOCATE(offset(N))
    ALLOCATE(stride(N))
    ALLOCATE(block(N))
    !
    offset = 0
    block = 1
    stride = 1
    !
    CALL h5open_f(error)
    if(error /= 0) write(*, '(a)') 'error1'
    !
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    if(error /= 0) write(*, '(a)') 'error2'  
    !
    CALL h5gopen_f(file_id, grpname, gset_id, error)
    if(error /= 0) write(*, '(a)') 'error3'
    !
    CALL h5dopen_f(gset_id, dsetname, dset_id, error)
    if(error /= 0) write(*, '(a)') 'error4'  
    ! Get dataset's dataspace identifier and select subset.
    !
    CALL h5dget_space_f(dset_id, dataspace, error)
    if(error /= 0) write(*, '(a)') 'Error in getting the dataspace'
    !
    DO i=1,N
      data_dims(i) = in_dims(i)
      sub_data_dims(i) = in_dims(i)
      count(i) = in_dims(i)
    ENDDO
    !
    count(6) = negnv_explrn
    sub_data_dims(6) = negnv_explrn
    !
    CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
         offset, count, error, stride, BLOCK) 
    if(error /= 0) write(*, '(a)') 'Error in getting the hyperslab'
    !
    ! Create memory dataspace.
    !
    CALL h5screate_simple_f(N, sub_data_dims, memspace, error)
    !
    data_out=0.d0
    !
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, sub_data_dims, error, &
      memspace, dataspace)
    !
    DEALLOCATE(data_dims)
    DEALLOCATE(sub_data_dims)
    DEALLOCATE(count)
    DEALLOCATE(offset)
    DEALLOCATE(stride)
    DEALLOCATE(block)
    !
    ! Close everything opened.
    !
    CALL h5sclose_f(dataspace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(gset_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
  !-----------------------------------------------------------------------    
  END SUBROUTINE READ_HDF5_SLAB
  !-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------
END MODULE exphonon
!-----------------------------------------------------------------------
