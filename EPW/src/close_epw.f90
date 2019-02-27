  !                                                                            
  ! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !----------------------------------------------------------------------
  MODULE close_epw
  !----------------------------------------------------------------------
  !! 
  !! This module contains routines related to deallocation and closing of files
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    ! 
    !----------------------------------------------------------------------------
    SUBROUTINE iter_close
    !----------------------------------------------------------------------------
    ! 
    ! This subroutine opens all the files needed to save scattering rates for the IBTE.
    ! 
    USE kinds,         ONLY : DP
    USE io_files,      ONLY : tmp_dir, prefix
    USE io_epw,        ONLY : iufilibtev_sup, iunepmat, iunsparseq, iunsparsek, &
                              iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                              iunsparsekcb, iunrestart, iunsparseicb, iunsparsejcb,&
                              iunsparsetcb, iunepmatcb, iunepmatwp2
    USE epwcom,        ONLY : iterative_bte, mp_mesh_k, int_mob, carrier, etf_mem, &
                              epmatkqread
    USE elph2,         ONLY : inv_tau_all, zi_allvb, inv_tau_allcb, zi_allcb,   &
                              map_rebal, map_rebal_inv    
    USE transportcom,  ONLY : s_BZtoIBZ_full, ixkqf_tr
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL
#endif
    ! 
    IMPLICIT NONE
    ! 
    ! Local variables
    INTEGER :: ierr
    !! Error status    
    ! 
    DEALLOCATE (inv_tau_all)
    DEALLOCATE (zi_allvb)
    IF (mp_mesh_k .AND. iterative_bte .AND. epmatkqread) DEALLOCATE (s_BZtoIBZ_full)
    IF (mp_mesh_k .AND. iterative_bte .AND. epmatkqread) DEALLOCATE (ixkqf_tr)
    IF (int_mob .AND. carrier) DEALLOCATE (inv_tau_allcb)
    IF (int_mob .AND. carrier) DEALLOCATE (zi_allcb)
    ! 
#if defined(__MPI)
    IF (etf_mem == 1) then
      CALL MPI_FILE_CLOSE(iunepmatwp2,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1 )
    ENDIF
    ! 
    IF (iterative_bte) THEN
      CALL MPI_FILE_CLOSE(iunepmat,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseq,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsek,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsei,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsej,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparset,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
  
      CALL MPI_FILE_CLOSE(iunepmatcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseqcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsekcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseicb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsejcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsetcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_close', 'error in MPI_FILE_CLOSE',1)
    ENDIF
#endif
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_close
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------
    SUBROUTINE deallocate_epw
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_epw
    !!  this routine is unchanged as of 3/9/09 and should be cleaned and fixed
    !!  09/2009 Cleanup still necessary
    !!  12/2009 Added variables from elph.f90
    !!
    !!  RM - Nov/Dec 2014
    !!  Imported the noncolinear case implemented by xlzhang
    !!
    !----------------------------------------------------------------------
    USE phcom,             ONLY : alphap, dmuxc, drc, dyn, evq, dvpsi, &
                                  int5, vlocq, int2_so, int5_so
    USE lrus,              ONLY : becp1, int3, int3_nc
    USE phus,              ONLY : int1, int1_nc, int2, int4, int4_nc
    USE lr_symm_base,      ONLY : rtau
    USE noncollin_module,  ONLY : m_loc
    USE control_lr,        ONLY : nbnd_occ
    USE becmod,            ONLY : becp, deallocate_bec_type
    USE elph2,             ONLY : el_ph_mat, epf17, epsi, etf,&
                                  etq, et_all, wf, wkf, wqf, &
                                  xkq, xk_all, zstar, xkf, xqf, epmatwp, eps_rpa
    USE epwcom,            ONLY : epbread, epwread
    USE modes,             ONLY : npert, u, name_rap_mode, num_rap_mode
    USE qpoint,            ONLY : eigqts, igkq 
    USE klist,             ONLY : nks
    !
    IMPLICIT NONE
    ! 
    INTEGER :: ik
    !! k-point number
    INTEGER :: ipol
    !! Polarization number
    !
    IF ( epwread .and. .not. epbread ) THEN
      !  EPW variables ONLY
      !
      IF(ALLOCATED(el_ph_mat)) DEALLOCATE (el_ph_mat)
      IF(ALLOCATED(epmatwp))   DEALLOCATE (epmatwp)
      IF(ALLOCATED(epf17))     DEALLOCATE (epf17)
      IF(ALLOCATED(etq))       DEALLOCATE (etq)
      IF(ALLOCATED(etf))       DEALLOCATE (etf)
      IF(ALLOCATED(wf))        DEALLOCATE (wf)
      IF(ALLOCATED(xkq))       DEALLOCATE (xkq)
      IF(ALLOCATED(xkf))       DEALLOCATE (xkf)
      IF(ALLOCATED(wkf))       DEALLOCATE (wkf)
      IF(ALLOCATED(xqf))       DEALLOCATE (xqf)
      IF(ALLOCATED(wqf))       DEALLOCATE (wqf)
      IF(ALLOCATED(xk_all))    DEALLOCATE (xk_all)
      IF(ALLOCATED(et_all))    DEALLOCATE (et_all)
      IF(ALLOCATED(eps_rpa))   DEALLOCATE (eps_rpa)
      IF(ALLOCATED(eps_rpa))   DEALLOCATE (eps_rpa) 
      ! 
    ELSE
      !   
      IF(ASSOCIATED(evq)) DEALLOCATE(evq)
      IF(ASSOCIATED(igkq)) DEALLOCATE(igkq)
      !
      IF(ALLOCATED(dvpsi)) DEALLOCATE (dvpsi)    
      !
      IF(ALLOCATED(vlocq)) DEALLOCATE (vlocq)
      IF(ALLOCATED(dmuxc)) DEALLOCATE (dmuxc)
      !
      IF(ALLOCATED(eigqts)) DEALLOCATE (eigqts)
      IF(ALLOCATED(rtau)) DEALLOCATE (rtau)
      IF(ASSOCIATED(u)) DEALLOCATE (u)
      if(allocated(name_rap_mode)) deallocate (name_rap_mode)
      if(allocated(num_rap_mode)) deallocate (num_rap_mode)
      IF(ALLOCATED(dyn)) DEALLOCATE (dyn)
      IF(ALLOCATED(epsi)) DEALLOCATE (epsi)
      IF(ALLOCATED(zstar)) DEALLOCATE (zstar)
      !
      IF(ALLOCATED(npert)) DEALLOCATE (npert)    
      !
      IF(ALLOCATED(int1)) DEALLOCATE (int1)    
      IF(ALLOCATED(int2)) DEALLOCATE (int2)
      IF(ALLOCATED(int3)) DEALLOCATE (int3)
      IF(ALLOCATED(int4)) DEALLOCATE (int4)
      IF(ALLOCATED(int5)) DEALLOCATE (int5)
      IF(ALLOCATED(int1_nc)) DEALLOCATE(int1_nc)
      IF(ALLOCATED(int3_nc)) DEALLOCATE(int3_nc)
      IF(ALLOCATED(int4_nc)) DEALLOCATE(int4_nc)
      IF(ALLOCATED(int2_so)) DEALLOCATE(int2_so)
      IF(ALLOCATED(int5_so)) DEALLOCATE(int5_so)
      ! 
      IF (allocated(alphap)) THEN
        DO ik = 1, nks
          DO ipol = 1, 3
            CALL deallocate_bec_type( alphap(ipol,ik) )
          ENDDO
        ENDDO
        DEALLOCATE(alphap)
      ENDIF
      IF (allocated(becp1)) THEN
        DO ik = 1, size(becp1)
          CALL deallocate_bec_type( becp1(ik) )
        ENDDO
        DEALLOCATE(becp1)
      ENDIF
      CALL deallocate_bec_type ( becp )
      !
      IF(ALLOCATED(nbnd_occ))  DEALLOCATE(nbnd_occ)
      IF(ALLOCATED(m_loc))     DEALLOCATE(m_loc)
      !
      IF(ALLOCATED(drc)) DEALLOCATE(drc)
      !
      !  EPW variables
      !
      IF(ALLOCATED(el_ph_mat)) DEALLOCATE (el_ph_mat)    
      IF(ALLOCATED(epmatwp))   DEALLOCATE (epmatwp)
      IF(ALLOCATED(epf17))     DEALLOCATE (epf17)    
      IF(ALLOCATED(etq))       DEALLOCATE (etq)    
      IF(ALLOCATED(etf))       DEALLOCATE (etf)    
      IF(ALLOCATED(wf))        DEALLOCATE (wf)    
      IF(ALLOCATED(xkq))       DEALLOCATE (xkq)    
      IF(ALLOCATED(xkf))       DEALLOCATE (xkf)    
      IF(ALLOCATED(wkf))       DEALLOCATE (wkf)    
      IF(ALLOCATED(xqf))       DEALLOCATE (xqf)    
      IF(ALLOCATED(wqf))       DEALLOCATE (wqf)    
      IF(ALLOCATED(xk_all))    DEALLOCATE (xk_all)    
      IF(ALLOCATED(et_all))    DEALLOCATE (et_all)    
      IF(ALLOCATED(eps_rpa))   DEALLOCATE (eps_rpa)
    ENDIF ! epwread .and. .not. epbread 
    !
    END SUBROUTINE deallocate_epw
    ! ---------------------------------------------------------------
    ! 
    !------------------------------------------------------------------
    SUBROUTINE close_final
    !------------------------------------------------------------------
    !
    USE units_lr,  ONLY : iuwfc
    USE units_ph,  ONLY : iudwf, iudrho
    USE phcom,     ONLY : fildrho
    USE mp_global, ONLY : me_pool,root_pool
    USE io_epw,    ONLY : iunepmatwe
    USE epwcom,    ONLY : etf_mem
    !
    implicit none
    !
    IF (etf_mem == 1 .OR. etf_mem == 2) THEN
      CLOSE (unit = iunepmatwe, status = 'delete')
    ENDIF
    !
    CLOSE (unit = iuwfc, status = 'keep')
    CLOSE (unit = iudwf, status = 'keep')
    IF (me_pool == root_pool ) THEN
      IF (fildrho.ne.' ') CLOSE (unit = iudrho, status = 'keep')
    ENDIF
    !
    END SUBROUTINE close_final
    ! ------------------------------------------------------------------
    ! 
  END MODULE close_epw
