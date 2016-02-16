!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_evc()
  !-----------------------------------------------------------------------
  !
  ! This routine is used only for isolated atoms in combination with
  ! the flag one_atom_occupations.
  ! It makes linear combinations of the degenerate bands so 
  ! that they have maximum overlap with the atomic states, and order 
  ! the bands in the same order as the atomic states.
  ! Weights "wg" must have been set to fixed values (as read in input)
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE basis,                ONLY : natomwfc, swfcatom
  USE klist,                ONLY : nks, ngk
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE control_flags,        ONLY : gamma_only, iverbosity
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer, save_buffer
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum

  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER :: ik, ibnd, jbnd, igroup, npw
  ! counter on k points
  !    "    "  bands
  !    "    "  groups of bands

  REAL(DP), EXTERNAL :: ddot
  COMPLEX(DP), EXTERNAL :: zdotc
  COMPLEX(DP), ALLOCATABLE :: proj(:,:), aux(:,:), aux_proj(:,:), a(:,:), v(:,:)

  REAL(DP) :: max_value, save_value, aux_et, maxproj
  INTEGER :: select_ibnd, iatwfc, first_available_band, info, nsize, &
             current_band, ngroups 
  INTEGER, ALLOCATABLE :: ind(:), group_size(:), start_band(:), used_atwfc(:)
  REAL(DP), ALLOCATABLE :: wband(:)

  IF (natomwfc > nbnd) THEN
     WRITE(6,'(5x,"natomwfc=", i5, " nbnd=",i5)') natomwfc, nbnd
     CALL errore('new_evc','increase nbnd',1)
  ENDIF

  ALLOCATE(proj(natomwfc,nbnd))  
  ALLOCATE(wband(nbnd))
  ALLOCATE(group_size(nbnd))
  ALLOCATE(start_band(nbnd))
  ALLOCATE(used_atwfc(nbnd))
  ALLOCATE(ind(nbnd))
  !
  !  we start a loop over k points
  !
  DO ik = 1, nks
     IF (lsda) current_spin = isk(ik)
     npw = ngk (ik)
     IF (nks > 1) &
        CALL get_buffer  (evc, nwordwfc, iunwfc, ik)

     CALL get_buffer (swfcatom, nwordatwfc, iunsat, ik)
     !
     ! make the projection on the atomic wavefunctions,
     !
     DO ibnd = 1, nbnd
        DO iatwfc = 1, natomwfc
           IF ( gamma_only ) THEN 
              proj (iatwfc, ibnd) = 2.d0 * &
                    ddot(2*npw, swfcatom (1, iatwfc), 1, evc (1, ibnd), 1) 
              IF (gstart.EQ.2) proj (iatwfc, ibnd) = proj (iatwfc, ibnd) - &
                    swfcatom (1, iatwfc) * evc (1, ibnd)
           ELSE 
              proj(iatwfc, ibnd) = zdotc(npw,swfcatom(1,iatwfc),1,evc(1,ibnd),1)
              IF (noncolin) &
                 proj (iatwfc, ibnd) = proj(iatwfc, ibnd) + &
                   zdotc (npw, swfcatom(npwx+1,iatwfc), 1, evc(npwx+1,ibnd), 1)
           ENDIF
        ENDDO
     ENDDO

     CALL mp_sum ( proj, intra_bgrp_comm )

     IF ( iverbosity > 0 ) THEN
        DO ibnd=1,nbnd
           WRITE(6,*) 'bands ', ibnd, et(ibnd,ik)*rytoev
           WRITE(6,'(8f9.4)') (ABS(proj(iatwfc,ibnd)), iatwfc=1,natomwfc)
        END DO
     END IF
!
!  We have to select natomwfc bands that have the largest overlap with
!  the atomic states. The other bands are empty and will be put above the
!  natomwfc with large projections.
!
     IF (natomwfc < nbnd) THEN
        DO ibnd=1,nbnd
           wband(ibnd) =0.0_DP
           DO iatwfc=1,natomwfc
              wband(ibnd) = wband(ibnd) + ABS(proj(iatwfc,ibnd))
           ENDDO
           ind(ibnd)=ibnd
        ENDDO
!
!   order from the largest to the smaller overlap
!
        wband=-wband
        CALL hpsort(nbnd, wband, ind)
!
!    now put the bands with smaller overlap above the others, change also
!    the eigenvalues and the projectors
!
        ALLOCATE(aux(npwx*npol,1))
        ALLOCATE(aux_proj(natomwfc,1))
        current_band=natomwfc+1
        DO ibnd =1, natomwfc
           IF (ind(ibnd) > natomwfc) THEN
              DO jbnd=current_band,nbnd
                 IF (ind(jbnd)<=natomwfc) THEN
                    aux(:,1)=evc(:,ind(ibnd))
                    evc(:,ind(ibnd))=evc(:,ind(jbnd))
                    evc(:,ind(jbnd))=aux(:,1)
                    aux_proj(:,1)=proj(:,ind(ibnd))
                    proj(:,ind(ibnd))=proj(:,ind(jbnd))
                    proj(:,ind(jbnd))=aux_proj(:,1)
                    aux_et = et(ind(ibnd),ik)
                    et(ind(ibnd),ik)=et(ind(jbnd),ik)
                    et(ind(jbnd),ik)=aux_et
                    current_band=jbnd+1
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        DEALLOCATE(aux)
        DEALLOCATE(aux_proj)
     ENDIF
!
!  Here we partition the bands in groups of degenerate bands.
!
     ngroups=1
     group_size=1
     start_band(1)=1
     DO iatwfc=1,natomwfc-1
        IF ( ABS(et(iatwfc,ik)-et(iatwfc+1,ik))>1.d-4) THEN
           ngroups=ngroups+1
           start_band(ngroups)=iatwfc+1
        ELSE
           group_size(ngroups) = group_size(ngroups) + 1
        ENDIF
     ENDDO 
!
!  For each group of bands we decide which are the atomic states
!  with the largest projection on the group of bands
!
     used_atwfc=0
     DO igroup = 1, ngroups
        DO iatwfc=1, natomwfc
           wband(iatwfc) = 0.0_DP
           DO ibnd = start_band(igroup), start_band(igroup)+group_size(igroup)-1
              wband(iatwfc) = wband(iatwfc) + ABS(proj(iatwfc,ibnd))
           ENDDO
           ind(iatwfc) = iatwfc
           IF (used_atwfc(iatwfc)==1)  wband(iatwfc)=0.0_DP
        ENDDO
!
!   order the atomic states from the largest to the smaller projection
!
        wband=-wband
        
        CALL hpsort(natomwfc, wband, ind)
        nsize = group_size(igroup)
!
        DO iatwfc=1,nsize
           IF (used_atwfc(ind(iatwfc))==1) THEN
              CALL errore('new_evc','this atomic wfc already used',ind(iatwfc))
           ELSE
               used_atwfc(ind(iatwfc))=1
           ENDIF
        ENDDO
!
!  At this point we solve a linear system of size group_size x group_size
!  and find the linear combination of degenerate wavefunctions which has 
!  projection one on each atomic state.
!
        IF (nsize>1) THEN
           ALLOCATE(aux(npwx*npol,nsize))
           ALLOCATE(aux_proj(natomwfc,nsize))
           ALLOCATE(a(nsize,nsize))
           ALLOCATE(v(nsize,nsize))
           v=(0.0_DP,0.0_DP)
           DO ibnd = 1, nsize
              DO jbnd = 1, nsize
                 a(ibnd,jbnd) = proj(ind(ibnd),start_band(igroup)+jbnd-1)
              ENDDO
              v(ibnd,ibnd)=(1.0_DP,0.0_DP)
           ENDDO
           CALL ZGESV(nsize, nsize, a, nsize, ind, v, nsize, info)
!
!  We cannot use the vectors v to make the linear combinations
!  because they are not orthonormal. We orthonormalize them, so the
!  projection will not be exactly one, but quite close.
!
           CALL orthogonalize_vects(nsize, v)
!
!  And now make the linear combination. Update also the projections on
!  the atomic states.
!
           aux=(0.0_DP, 0.0_DP)
           aux_proj=(0.0_DP, 0.0_DP)
           DO ibnd=1, nsize 
              DO jbnd=1,nsize
                 aux(:,ibnd)=aux(:,ibnd)+ v(jbnd,ibnd)* &
                                     evc(:,start_band(igroup)+jbnd-1)
                 aux_proj(:,ibnd)=aux_proj(:,ibnd)+ v(jbnd,ibnd)* &
                                     proj(:,start_band(igroup)+jbnd-1)
              ENDDO
           ENDDO
           evc(:,start_band(igroup):start_band(igroup)+nsize-1)= aux(:,:)
           proj(:,start_band(igroup):start_band(igroup)+nsize-1)= aux_proj(:,:)
           DEALLOCATE(aux)
           DEALLOCATE(aux_proj)
           DEALLOCATE(a)
           DEALLOCATE(v)
        ENDIF
     ENDDO   ! loop over the groups of bands
!
!  Finally, we order the new bands as the atomic states
!
     ALLOCATE(aux(npwx*npol,natomwfc))
     used_atwfc=0
     DO ibnd=1,natomwfc
        current_band=1
        maxproj=0.0_DP
        DO iatwfc=1,natomwfc
           IF (ABS(proj(iatwfc,ibnd))>maxproj.AND.used_atwfc(iatwfc)==0) THEN
              current_band=iatwfc
              maxproj=ABS(proj(iatwfc,ibnd))
           ENDIF
        ENDDO
        used_atwfc(current_band)=1
        aux(:,current_band)=evc(:,ibnd) 
        wband(current_band)=et(ibnd,ik)
     ENDDO
     evc(:,1:natomwfc)=aux(:,:)
     et(1:natomwfc,ik)=wband(1:natomwfc)
     DEALLOCATE(aux)
!
!  If needed save the new bands on disk
!
     IF (nks > 1) THEN
        CALL save_buffer  (evc, nwordwfc, iunwfc, ik)
     END IF
  ENDDO

  DEALLOCATE(group_size)
  DEALLOCATE(start_band)
  DEALLOCATE(ind)
  DEALLOCATE(wband)
  DEALLOCATE(used_atwfc)
  DEALLOCATE(proj)

  RETURN

END SUBROUTINE new_evc

SUBROUTINE orthogonalize_vects(n,v)
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
COMPLEX(DP), INTENT(INOUT) :: v(n,n)
COMPLEX(DP) :: sca
REAL(DP) :: norm
INTEGER :: i,k
COMPLEX(DP), EXTERNAL :: zdotc
REAL(DP), EXTERNAL :: ddot

norm=ddot(2*n,v(:,1),1,v(:,1),1)
v(:,1)=v(:,1)/SQRT(norm)

DO i=2,n
   DO k=i-1, 1, -1
      sca=zdotc(n, v(:,k),1, v(:,i),1 )
      v(:,i)=v(:,i) - sca * v(:,k)
   ENDDO
   norm=ddot(2*n,v(:,i),1,v(:,i),1)
   v(:,i)=v(:,i)/SQRT(norm)
ENDDO
   
RETURN
END SUBROUTINE orthogonalize_vects
