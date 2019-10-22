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
    !! 
    !! This subroutine opens all the files needed to save scattering rates for the IBTE.
    !! 
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iunepmat, iunsparseq,  &
                              iunsparseqcb, iunepmatcb, iunepmatwp2
    USE epwcom,        ONLY : iterative_bte, int_mob, carrier, etf_mem, assume_metal
    USE epwcom,        ONLY : int_mob, carrier, ncarrier
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
#if defined(__MPI)
    IF (etf_mem == 1) then
      CALL MPI_FILE_CLOSE(iunepmatwp2, ierr)
      IF (ierr /= 0) CALL errore('iter_close', 'error in MPI_FILE_CLOSE', 1)
    ENDIF
#endif
    ! 
    IF (iterative_bte) THEN
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier < 0.0)) .OR. assume_metal) THEN
        CLOSE(iunepmat)
        CLOSE(iunsparseq)
      ENDIF
      IF ((int_mob .AND. carrier) .OR. ((.NOT. int_mob .AND. carrier) .AND. (ncarrier > 0.0)) .AND. .NOT. assume_metal) THEN
        CLOSE(iunepmatcb)
        CLOSE(iunsparseqcb)
      ENDIF  ! in all other cases it is still to decide which files to open
    ENDIF
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_close
    !----------------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_epw()
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
    USE phcom,             ONLY : drc, dyn
    USE noncollin_module,  ONLY : m_loc
    USE control_lr,        ONLY : nbnd_occ
    USE elph2,             ONLY : epsi, etf, wkf, wqf, &
                                  zstar, xkf, xqf
    USE klist_epw,         ONLY : xk_all, xk_loc, xk_cryst, et_all, et_loc, & 
                                  isk_loc, isk_all
    USE epwcom,            ONLY : epbread, epwread
    USE qpoint,            ONLY : igkq 
    !
    IMPLICIT NONE
    ! 
    IF (epwread .AND. .NOT. epbread) THEN
      !  EPW variables only
      !
      IF(ALLOCATED(etf))       DEALLOCATE(etf)
      IF(ALLOCATED(xkf))       DEALLOCATE(xkf)
      IF(ALLOCATED(wkf))       DEALLOCATE(wkf)
      IF(ALLOCATED(xqf))       DEALLOCATE(xqf)
      IF(ALLOCATED(wqf))       DEALLOCATE(wqf)
      IF(ALLOCATED(et_all))    DEALLOCATE(et_all)
      ! 
    ELSE
      !   
      IF(ASSOCIATED(igkq))    DEALLOCATE(igkq)
      IF(ALLOCATED(dyn))      DEALLOCATE(dyn)
      IF(ALLOCATED(epsi))     DEALLOCATE(epsi)
      IF(ALLOCATED(zstar))    DEALLOCATE(zstar)
      IF(ALLOCATED(nbnd_occ)) DEALLOCATE(nbnd_occ)
      IF(ALLOCATED(m_loc))    DEALLOCATE(m_loc)
      IF(ALLOCATED(drc))      DEALLOCATE(drc)
      !
      !  EPW variables
      !
      IF(ALLOCATED(etf))       DEALLOCATE(etf)    
      IF(ALLOCATED(xkf))       DEALLOCATE(xkf)    
      IF(ALLOCATED(wkf))       DEALLOCATE(wkf)    
      IF(ALLOCATED(xqf))       DEALLOCATE(xqf)    
      IF(ALLOCATED(wqf))       DEALLOCATE(wqf)    
      IF(ALLOCATED(xk_all))    DEALLOCATE(xk_all)
      IF(ALLOCATED(xk_loc))    DEALLOCATE(xk_loc)
      IF(ALLOCATED(xk_cryst))  DEALLOCATE(xk_cryst)
      IF(ALLOCATED(et_all))    DEALLOCATE(et_all)    
      IF(ALLOCATED(et_loc))    DEALLOCATE(et_loc)    
      IF(ALLOCATED(isk_loc))   DEALLOCATE(isk_loc)    
      IF(ALLOCATED(isk_all))   DEALLOCATE(isk_all)    
    ENDIF ! epwread .AND. .NOT. epbread 
    !
    !---------------------------------------------------------------
    END SUBROUTINE deallocate_epw
    !---------------------------------------------------------------
    ! 
    !------------------------------------------------------------------
    SUBROUTINE close_final
    !------------------------------------------------------------------
    !
    USE units_lr,  ONLY : iuwfc
    USE units_ph,  ONLY : iudwf, iudrho
    USE phcom,     ONLY : fildrho
    USE mp_global, ONLY : me_pool,root_pool
    USE io_var,    ONLY : iunepmatwe
    USE epwcom,    ONLY : etf_mem
    !
    IMPLICIT NONE
    !
    IF (etf_mem == 1 .OR. etf_mem == 2) THEN
      CLOSE(UNIT = iunepmatwe, STATUS = 'delete')
    ENDIF
    !
    CLOSE(UNIT = iuwfc, STATUS = 'keep')
    CLOSE(UNIT = iudwf, STATUS = 'keep')
    IF (me_pool == root_pool) THEN
      IF (fildrho/=' ') CLOSE(UNIT = iudrho, STATUS = 'keep')
    ENDIF
    !
    !------------------------------------------------------------------
    END SUBROUTINE close_final
    !------------------------------------------------------------------
    ! 
  !------------------------------------------------------------------
  END MODULE close_epw
  !------------------------------------------------------------------
