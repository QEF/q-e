!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  this initializes som variables for DFT+U calculation,
!  espescially atomic wfc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


subroutine init_xanes_ldau
  USE ldaU,             ONLY : lda_plus_u, Hubbard_U, Hubbard_l, &
                                 Hubbard_alpha, Hubbard_lmax, U_projection 
  USE basis,            only : natomwfc  
  USE xspectra,            only : xread_wf, U_projection_type
  use wvfct,            ONLY : npwx,nbndx,nbnd,npw
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE lsda_mod,          ONLY : lsda, nspin
  USE uspp_param,           ONLY : upf
  implicit none

  logical :: exst, opnd
  integer :: ldim, nt
  INTEGER :: set_Hubbard_l


  call init_at_1


     Hubbard_lmax = -1
     !
     DO nt = 1, ntyp
        !
        IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
           !
            Hubbard_l(nt) = set_Hubbard_l( upf(nt)%psd )
           !
           Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
           !
           WRITE( UNIT = stdout, &
                  FMT = * ) ' HUBBARD L FOR TYPE ',upf(nt)%psd,' IS ', Hubbard_l(nt)
           !
        END IF
        !
     END DO
     !
     WRITE( UNIT = stdout, &
            FMT = * ) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
     !
     IF ( Hubbard_lmax == -1 ) &
        CALL errore( 'setup', &
                   & 'lda_plus_u calculation but Hubbard_l not set', 1 )

  ldim = 2 * Hubbard_lmax + 1
  U_projection=U_projection_type

 end subroutine init_xanes_ldau




! Here is a slightly modified version of orthoatwfc, for a 
! selected k-point and without writing the result
! modified by CG


!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_xanes_ldau_2(ik)
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunat, iunsat, nwordatwfc, iunigk, diropn
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc
  USE klist,      ONLY : nks, xk, ngk
  USE ldaU,       ONLY : swfcatom, U_projection
  USE wvfct,      ONLY : npwx, npw, igk
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  ! 
  IMPLICIT NONE
  !
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol
  ! the k point under consideration
  ! counter on bands
  REAL(DP) :: t0, scnds
  ! cpu time spent
  LOGICAL :: orthogonalize_wfc
     
  COMPLEX(DP) :: temp, t (5)
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:), work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)
  LOGICAL :: opnd, exst

  t0 = scnds ()
  
  IF (noncolin) THEN
     ALLOCATE (wfcatom( npwx*npol, natomwfc))    
  ELSE
     ALLOCATE (wfcatom( npwx, natomwfc))    
  END IF
  ALLOCATE (overlap( natomwfc , natomwfc))    
  ALLOCATE (work   ( natomwfc , natomwfc))    
  ALLOCATE (e      ( natomwfc))    

  IF (U_projection=="file") THEN
     WRITE( stdout,*) 'LDA+U Projector read from file '

       nwordatwfc=2*npwx*natomwfc*npol

       CALL diropn( iunsat, 'satwfc', nwordatwfc, exst )
       CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )

       INQUIRE( UNIT = iunsat, OPENED = opnd )
       IF ( opnd ) CLOSE( UNIT = iunsat, STATUS = 'KEEP' )

     RETURN
  END IF

  IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are orthogonalized'
     IF (gamma_only) THEN
        WRITE( stdout,*) 'Gamma-only calculation for this case not implemented'
        STOP
     END IF
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) THEN
        WRITE( stdout,*) 'Gamma-only calculation for this case not implemented'
        STOP
     END IF
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("orthoatwfc"," this U_projection_type is not valid",1)
  END IF

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  
     
     overlap(:,:) = (0.d0,0.d0)
     work(:,:) = (0.d0,0.d0)
     
     CALL atomic_wfc (ik, wfcatom)
     
     CALL calbec (npw, vkb, wfcatom, becp)

     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

   IF (orthogonalize_wfc) THEN
     !
     ! calculate overlap matrix
     !
     IF (noncolin) THEN
        CALL zgemm ('c', 'n', natomwfc, natomwfc, npwx*npol, (1.d0, 0.d0), &
             wfcatom, npwx, swfcatom, npwx, (0.d0,0.d0), overlap, natomwfc)
     ELSE
         CALL zgemm ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0), &
             wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0), overlap, natomwfc)
     END IF
#ifdef __MPI
     CALL mp_sum(  overlap, intra_pool_comm )
#endif
     IF (U_projection=="norm-atomic") THEN
        DO i = 1, natomwfc
           DO j = i+1, natomwfc
              overlap(i,j) = (0.d0,0.d0)
              overlap(j,i) = (0.d0,0.d0)
           ENDDO
        ENDDO
     END IF
     !
     ! find O^-.5
     !
     CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
     DO i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     ENDDO
     DO i = 1, natomwfc
        DO j = i, natomwfc
           temp = (0.d0, 0.d0)
           DO k = 1, natomwfc
              temp = temp + e (k) * work (j, k) * CONJG (work (i, k) )
           ENDDO
           overlap (i, j) = temp
           IF (j.NE.i) overlap (j, i) = CONJG (temp)
        ENDDO
     ENDDO
     !
     ! transform atomic orbitals O^-.5 psi
     !
     DO i = 1, npw
        work(:,1) = (0.d0,0.d0)
        IF (noncolin) THEN
           DO ipol=1,npol
              j = i + (ipol-1)*npwx
              CALL zgemv ('n',natomwfc,natomwfc,(1.d0,0.d0),overlap, &
                   natomwfc,swfcatom(j,1),npwx*npol, &
                                       (0.d0,0.d0),work,1)
              CALL zcopy (natomwfc,work,1,swfcatom(j,1),npwx*npol)
           END DO
        ELSE
           CALL zgemv ('n', natomwfc, natomwfc, (1.d0, 0.d0) , overlap, &
                natomwfc, swfcatom (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
           CALL zcopy (natomwfc, work, 1, swfcatom (i, 1), npwx)
        END IF
     ENDDO
        
   END IF ! orthogonalize_wfc

  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
  DEALLOCATE (wfcatom)
  CALL deallocate_bec_type (becp )
  !
  RETURN
     
END SUBROUTINE init_xanes_ldau_2


