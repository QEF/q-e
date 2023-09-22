!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE  write_xml_proj (filename, projs, lwrite_ovp, ovps )
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iunpun, restart_dir
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE lsda_mod,         ONLY : nspin
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE xmltools,         ONLY : xml_open_file, xml_closefile, xmlw_writetag,&
                               xmlw_opentag, xmlw_closetag, add_attr
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  COMPLEX(DP),   INTENT(IN) :: ovps(natomwfc,natomwfc,nkstot)
  LOGICAL,       INTENT(IN) :: lwrite_ovp
  !
  COMPLEX(DP) :: proj(nbnd)
  INTEGER :: ik, ik_eff, is, nwfc, ibnd, nspin_lsda, num_k_points
  !
  !
  iunpun = xml_open_file ( TRIM( restart_dir() )//TRIM(filename) )
  IF ( iunpun == -1 ) RETURN
  !
  nspin_lsda = 1
  IF ( nspin == 2 ) nspin_lsda = 2
  num_k_points = nkstot / nspin_lsda
  !
  ! <PROJECTIONS>
  !
  CALL xmlw_opentag ("PROJECTIONS")
  !
  ! <HEADER>
  !
  CALL add_attr ("NUMBER_OF_BANDS", nbnd)
  CALL add_attr ("NUMBER_OF_K-POINTS", num_k_points)
  CALL add_attr ("NUMBER_OF_SPIN_COMPONENTS", nspin_lsda)
  CALL add_attr ("NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL add_attr ("NUMBER_OF_ELECTRONS", nelec)
  CALL add_attr ("FERMI_ENERGY", ef)
  CALL xmlw_writetag ("HEADER", "")
  !
  ! <EIGENSTATES>
  !
  CALL xmlw_opentag ("EIGENSTATES")
  DO is = 1, nspin_lsda
     DO ik = 1, num_k_points
        ik_eff = ik + (is-1)*num_k_points
        CALL add_attr ( "Weight", wk(ik_eff) )
        CALL xmlw_writetag("K-POINT", xk(:,ik_eff) )
        !
        CALL xmlw_writetag ( "E", et(:,ik_eff) )
        !
        CALL xmlw_opentag ("PROJS")
        DO nwfc = 1, natomwfc
           CALL add_attr ( "index", nwfc )
           CALL add_attr ( "spin", is )
           ! NOTE: the complex to real conversion done inside xmlw_writetag
           !       using C pointer does not work  on intel compilers
           !       with an array section (non-contiguous memory)
           ! CALL xmlw_writetag ("ATOMIC_WFC", projs(nwfc,:,ik_eff) )
           proj(:) =  projs(nwfc,:,ik_eff)
           CALL xmlw_writetag ("ATOMIC_WFC", proj )
        ENDDO
        CALL xmlw_closetag ( )
        !
     END DO
     !
  ENDDO
  CALL xmlw_closetag ( )
  !
  ! </EIGENSTATE>
  ! <OVERLAPS> if required
  IF ( lwrite_ovp ) THEN
     CALL xmlw_opentag("OVERLAPS" )
     DO ik = 1, num_k_points
        DO is = 1, nspin_lsda
           ik_eff = ik + (is-1)*num_k_points
           CALL add_attr ("dim", natomwfc )
           CALL add_attr ("spin", is )
           CALL xmlw_writetag ("OVPS", ovps(:,:,ik_eff) )
        END DO
     END DO
     CALL xmlw_closetag ( )
  END IF
  !
  CALL xmlw_closetag ( )
  !
  ! </PROJECTIONS>
  !
  CALL xml_closefile ( )
  !
END SUBROUTINE write_xml_proj
!
!-----------------------------------------------------------------------
SUBROUTINE write_proj_file ( filproj, proj )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module, ONLY : noncolin, lspinorb
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : xk, nkstot
  USE run_info,  ONLY : title
  USE cell_base, ONLY : at, ibrav, celldm
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, atm
  USE wvfct,     ONLY : nbnd
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : gcutm 
  USE gvecs,     ONLY : dual
  USE gvecw,     ONLY : ecutwfc
  USE projections, ONLY : nlmchi, compute_mj
  !
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: proj(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256) :: filename
  INTEGER :: is, ik, nwfc, ibnd, nk_, nksinit, nkslast, iunproj=33
  !
  IF ( TRIM(filproj) == ' ' ) RETURN
  !
  DO is=1,nspin
     IF (nspin==2) THEN
        IF (is==1) filename=trim(filproj)//'.projwfc_up'
        IF (is==2) filename=trim(filproj)//'.projwfc_down'
        nksinit=(nkstot/2)*(is-1)+1
        nkslast=(nkstot/2)*is
        nk_ = nkstot/2
     ELSE
        filename=trim(filproj)//'.projwfc_up'
        nksinit=1
        nkslast=nkstot
        nk_ = nkstot
     ENDIF
     CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
          dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
          celldm, at, gcutm, dual, ecutwfc, nk_, nbnd, natomwfc)
     DO nwfc = 1, natomwfc
        IF (lspinorb) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,2i5,f5.1,1x,f5.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%jj, &
                compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l, nlmchi(nwfc)%m)
        ELSE IF (noncolin) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5,1x,f4.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%m, &
                0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
        ELSE
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l,  &
                nlmchi(nwfc)%m
        END IF
        DO ik=nksinit,nkslast
           DO ibnd=1,nbnd
              WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, abs(proj(nwfc,ibnd,ik))
           ENDDO
        ENDDO
     ENDDO
     CLOSE(iunproj)
  ENDDO
  !
END SUBROUTINE write_proj_file
