  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/deallocate_phq - Quantum-ESPRESSO group
  !
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
  USE phcom,             ONLY : alphap, alphasum, alphasum_nc, &
                                becsum_nc, dmuxc, dpsi,&
                                drc, dpsi, dyn, evq, dvpsi,&
                                int5, vlocq, int2_so, int5_so
  USE lrus,              ONLY : becp1, int3, int3_nc
  USE phus,              ONLY : int1, int1_nc, int2, int4, int4_nc
  USE lr_symm_base,      ONLY : rtau
  USE noncollin_module,  ONLY : m_loc
  USE control_lr,        ONLY : lgamma, nbnd_occ
  USE becmod,            ONLY : becp, deallocate_bec_type
  USE elph2,             ONLY : el_ph_mat, epf17, epsi, etf,&
                                etq, et_all, wf, wkf, wqf, wslen,&
                                xkq, xk_all, zstar, xkf, xqf, epmatwp
  USE epwcom,            ONLY : epbread, epwread
  USE modes,             ONLY : t, npert, u, name_rap_mode, num_rap_mode
  USE qpoint,            ONLY : eigqts, igkq 
  USE klist,             ONLY : nks
  !
  IMPLICIT NONE
  INTEGER :: ik, ipol
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
    IF(ALLOCATED(wslen))     DEALLOCATE (wslen)
    ! 
  ELSE
    !   
    IF (lgamma) THEN
       IF(ASSOCIATED(evq)) NULLIFY(evq)
       IF(ASSOCIATED(igkq)) NULLIFY(igkq)
    ELSE
       IF(ASSOCIATED(evq)) DEALLOCATE(evq)
       IF(ASSOCIATED(igkq)) DEALLOCATE(igkq)
    END IF
    !
    IF(ALLOCATED(dvpsi)) DEALLOCATE (dvpsi)    
    IF(ALLOCATED(dpsi)) DEALLOCATE ( dpsi)    
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
    !IF(ASSOCIATED(t)) DEALLOCATE (t)
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
    IF(ALLOCATED(becsum_nc)) DEALLOCATE(becsum_nc)
    IF(ALLOCATED(alphasum_nc)) DEALLOCATE(alphasum_nc)
    IF(ALLOCATED(int2_so)) DEALLOCATE(int2_so)
    IF(ALLOCATED(int5_so)) DEALLOCATE(int5_so)
    IF(ALLOCATED(alphasum)) DEALLOCATE (alphasum)
    ! 
    if(allocated(alphap)) then
       do ik=1,nks
          do ipol=1,3
             call deallocate_bec_type ( alphap(ipol,ik) )
          enddo
       end do
       deallocate (alphap)
    endif
    if(allocated(becp1))  then
       do ik=1,size(becp1)
          call deallocate_bec_type ( becp1(ik) )
       end do
       deallocate(becp1)
    end if
    call deallocate_bec_type ( becp )
  
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
    IF(ALLOCATED(wslen))     DEALLOCATE (wslen)    
  ENDIF ! epwread .and. .not. epbread 
  !
  END SUBROUTINE deallocate_epw
