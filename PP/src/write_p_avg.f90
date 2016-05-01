!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE write_p_avg(filp, spin_component, firstk, lastk)
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : ngm, nl, g
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef
  USE wvfct,                ONLY : et, nbnd, npwx
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE noncollin_module,     ONLY : noncolin, npol
  USE ldaU,                 ONLY : lda_plus_u
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_world,             ONLY : world_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: spin_component, nks1, nks2, firstk, lastk, npw
  INTEGER :: iunout, ios, ik, ibnd, jbnd, ipol, nbnd_occ
  COMPLEX(DP) :: zdotc
  COMPLEX(DP), ALLOCATABLE :: ppsi(:,:), ppsi_us(:,:), matp(:,:,:)
  CHARACTER (len=256) :: filp, namefile
  !
  IF (lda_plus_u) CALL errore('write_p_avg', &
                       'write_p_avg not working with LDA+U',1)
  ALLOCATE(matp(nbnd,nbnd,3))
  CALL allocate_bec_type ( nkb, nbnd, becp)

  IF (nspin==1.or.nspin==4) THEN
     nks1=max(1,firstk)
     nks2=min(nkstot, lastk)
     IF (spin_component /= 1)  &
        CALL errore('write_p_avg','incorrect spin_component',1)
  ELSEIF (nspin==2) THEN
     IF (spin_component == 1) THEN
        nks1=max(1,firstk)
        nks2=min(nks/2,lastk)
     ELSEIF (spin_component==2) THEN
        nks1=nks/2 + max(1,firstk)
        nks2=nks/2 + min(nks/2,lastk)
     ELSE
        CALL errore('write_p_avg','incorrect spin_component',1)
     ENDIF
  ENDIF

  ios = 0
  IF ( ionode ) THEN
     iunout=58
     namefile=trim(filp)
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast (ios, ionode_id, world_comm)
  IF ( ios/=0 ) CALL errore ('write_p_avg', 'Opening filband file', abs (ios) )

  DO ik = nks1, nks2
     !
     !   Compute the number of occupated bands at this k point
     !
     DO ibnd = 1, nbnd
        IF (et (ibnd, ik)<=ef) nbnd_occ = ibnd
     ENDDO
     IF (nbnd_occ==nbnd) WRITE( stdout, '(5x,/,&
             &"No empty band at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     !
     ALLOCATE(ppsi(npwx*npol,nbnd_occ))
     IF (okvan) ALLOCATE(ppsi_us(npwx*npol,nbnd_occ))
     !
     npw = ngk(ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     CALL calbec ( npw, vkb, evc, becp, nbnd_occ )

     DO ipol=1,3
        CALL compute_ppsi(ppsi, ppsi_us, ik, ipol, nbnd_occ, spin_component)
        DO ibnd=nbnd_occ+1,nbnd
           DO jbnd=1,nbnd_occ
              IF (noncolin) THEN
                 matp(ibnd-nbnd_occ,jbnd,ipol)=  &
                         zdotc(npwx*npol,evc(1,ibnd),1,ppsi(1,jbnd),1)
                 IF (okvan) THEN
                    matp(ibnd-nbnd_occ,jbnd,ipol)=                  &
                         matp(ibnd-nbnd_occ,jbnd,ipol)+             &
                           (0.d0,0.5d0)*(et(ibnd,ik)-et(jbnd,ik))*  &
                         (zdotc(npwx*npol,evc(1,ibnd),1,ppsi_us(1,jbnd),1) )
                 ENDIF
              ELSE
                 matp(ibnd-nbnd_occ,jbnd,ipol)=  &
                            zdotc(npw,evc(1,ibnd),1,ppsi(1,jbnd),1)
                 IF (okvan) THEN
                    matp(ibnd-nbnd_occ,jbnd,ipol)= &
                               matp(ibnd-nbnd_occ,jbnd,ipol) +  &
                   (0.d0,0.5d0)*zdotc(npw,evc(1,ibnd),1,ppsi_us(1,jbnd),1)* &
                   (et(ibnd,ik)-et(jbnd,ik))

                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(ppsi)
     IF (okvan) DEALLOCATE(ppsi_us)
     CALL mp_sum(matp, intra_bgrp_comm)

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
                    (abs(matp(ibnd-nbnd_occ,jbnd,ipol))**2, jbnd=1,nbnd_occ)

           ENDDO
        ENDDO
     ENDIF
  ENDDO

  IF (ionode) THEN
     CLOSE(iunout)
  ENDIF

  DEALLOCATE(matp)
  !
  RETURN
END SUBROUTINE write_p_avg
