!
! Copyright (C) 2006-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE write_p_avg(filp, spin_component, firstk, lastk)
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba2, at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : nrx1, nrx2, nrx3, nrxx, nr1, nr2, &
                                   nr3, ngm, nl, g, ecutwfc
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef
  USE wvfct,                ONLY : et, nbnd, npwx, npw, igk, g2kin
  USE klist,                ONLY : xk, nks, nkstot
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE becmod,               ONLY : becp, becp_nc, calbec, &
                                   allocate_bec, deallocate_bec
  USE noncollin_module,     ONLY : noncolin, npol
  USE ldaU,                 ONLY : lda_plus_u
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: spin_component, nks1, nks2, firstk, lastk
  INTEGER :: iunout, ios, ik, ibnd, jbnd, ipol, nbnd_occ
  COMPLEX(DP) :: ZDOTC
  COMPLEX(DP), ALLOCATABLE :: ppsi(:,:), ppsi_us(:,:), matp(:,:,:) 
  CHARACTER (LEN=256) :: filp, namefile
  !
  IF (lda_plus_u) CALL errore('write_p_avg', &
                       'write_p_avg not working with LDA+U',1)
  ALLOCATE(matp(nbnd,nbnd,3))
  CALL allocate_bec (nkb,nbnd)

  IF (nspin==1.OR.nspin==4) THEN
     nks1=MAX(1,firstk)
     nks2=MIN(nkstot, lastk)
     IF (spin_component .ne. 1)  &
        CALL errore('punch_bands','uncorrect spin_component',1)
  ELSE IF (nspin.eq.2) THEN
     IF (spin_component == 1) THEN
        nks1=MAX(1,firstk)
        nks2=MIN(nks/2,lastk)
     ELSE IF (spin_component==2) THEN
        nks1=nks/2 + MAX(1,firstk)
        nks2=nks/2 + MIN(nks/2,lastk)
     ELSE
        CALL errore('punch_bands','uncorrect spin_component',1)
     END IF
  END IF

  ios = 0
  IF ( ionode ) THEN
     iunout=58
     namefile=TRIM(filp)
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  END IF

  CALL mp_bcast (ios, ionode_id)
  IF ( ios/=0 ) CALL errore ('sym_band', 'Opening filband file', ABS (ios) )

  DO ik = nks1, nks2
     !
     !   Compute the number of occupated bands at this k point
     !
     DO ibnd = 1, nbnd
        IF (et (ibnd, ik).LE.ef) nbnd_occ = ibnd
     END DO
     IF (nbnd_occ==nbnd) WRITE( stdout, '(5x,/,&
             &"No empty band at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     !
     ALLOCATE(ppsi(npwx*npol,nbnd_occ))
     IF (okvan) ALLOCATE(ppsi_us(npwx*npol,nbnd_occ))
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)

     IF (noncolin) THEN
        CALL calbec ( npw, vkb, evc, becp_nc, nbnd_occ )
     ELSE
        CALL calbec ( npw, vkb, evc, becp, nbnd_occ )
     END IF
 
     DO ipol=1,3
        CALL compute_ppsi(ppsi, ppsi_us, ik, ipol, nbnd_occ, spin_component)
        DO ibnd=nbnd_occ+1,nbnd
           DO jbnd=1,nbnd_occ
              IF (noncolin) THEN
                 matp(ibnd-nbnd_occ,jbnd,ipol)=  &
                         ZDOTC(npwx*npol,evc(1,ibnd),1,ppsi(1,jbnd),1) 
                 IF (okvan) THEN
                    matp(ibnd-nbnd_occ,jbnd,ipol)=                  &
                         matp(ibnd-nbnd_occ,jbnd,ipol)+             &
                           (0.d0,0.5d0)*(et(ibnd,ik)-et(jbnd,ik))*  &
                         (ZDOTC(npwx*npol,evc(1,ibnd),1,ppsi_us(1,jbnd),1) )
                 END IF
              ELSE
                 matp(ibnd-nbnd_occ,jbnd,ipol)=  &
                            ZDOTC(npw,evc(1,ibnd),1,ppsi(1,jbnd),1)
                 IF (okvan) THEN
                    matp(ibnd-nbnd_occ,jbnd,ipol)= &
                               matp(ibnd-nbnd_occ,jbnd,ipol) +  &
                   (0.d0,0.5d0)*ZDOTC(npw,evc(1,ibnd),1,ppsi_us(1,jbnd),1)* &
                   (et(ibnd,ik)-et(jbnd,ik)) 
                    
                 ENDIF
              END IF
           END DO
        END DO
     END DO
     DEALLOCATE(ppsi)
     IF (okvan) DEALLOCATE(ppsi_us)
#ifdef __PARA
     CALL reduce(nbnd*nbnd*3*2,matp)
#endif

     IF (ionode) THEN
        IF (ik == nks1) &
           WRITE (iunout, '(" &p_mat nbnd=",i4,", nks=",i4," /")') &
                 nbnd, nks2-nks1+1
        WRITE (iunout, '(10x,3f10.6,i7)') xk(1,ik),xk(2,ik),xk(3,ik), &
                                          nbnd_occ
        
        DO ipol=1,3
           WRITE (iunout, '(i3)') ipol
           DO ibnd=nbnd_occ+1,nbnd
              WRITE (iunout, '(5f15.8)') &
                    (ABS(matp(ibnd-nbnd_occ,jbnd,ipol))**2, jbnd=1,nbnd_occ)
           
           END DO
        END DO
     END IF
  END DO

  IF (ionode) THEN
     CLOSE(iunout)
  END IF

  DEALLOCATE(matp)
  !
  RETURN
END SUBROUTINE write_p_avg
