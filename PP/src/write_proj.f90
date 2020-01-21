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
  USE FoX_wxml,         ONLY : xmlf_t, xml_openfile, xml_close, &
       xml_newElement, xml_addAttribute, xml_addCharacters, xml_endElement
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  LOGICAL,       INTENT(IN) :: lwrite_ovp
  REAL(DP),      INTENT(IN) :: ovps(2*natomwfc,natomwfc,nkstot)
  ! slightly dirty trick to make ovps a real array
  !
  INTEGER :: ik, ik_eff, is, nwfc, ibnd, ierr, nspin_lsda, num_k_points
  REAL(DP):: proj_tmp(2*nbnd) 
  TYPE(xmlf_t) :: xf
  !
  !
  CALL xml_openfile(FILENAME = TRIM( restart_dir() )//TRIM(filename), XF = xf, &
       UNIT = iunpun, PRETTY_PRINT = .TRUE., REPLACE  = .TRUE., &
       NAMESPACE = .FALSE., IOSTAT = ierr ) 
  !
  IF ( ierr /= 0 ) RETURN
  !
  nspin_lsda = 1
  IF ( nspin == 2 ) nspin_lsda = 2
  num_k_points = nkstot / nspin_lsda
  !
  ! <PROJECTIONS>
  !
  CALL xml_newElement (xf, "PROJECTIONS")
  !
  ! <HEADER>
  !
  CALL xml_newElement (xf, "HEADER")
  CALL xml_addAttribute (xf, "NUMBER_OF_BANDS", nbnd)
  CALL xml_addAttribute (xf, "NUMBER_OF_K-POINTS", num_k_points)
  CALL xml_addAttribute (xf, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL xml_addAttribute (xf, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL xml_addAttribute (xf, "NUMBER_OF_ELECTRONS", nelec)
  CALL xml_addAttribute (xf, "FERMI_ENERGY", ef)
  CALL xml_endElement (xf, "HEADER")
  !
  ! </HEADER>
  ! <EIGENSTATES>
  !
  CALL xml_newElement (xf, "EIGENSTATES")
  DO ik = 1, num_k_points
     DO is = 1, nspin_lsda
        ik_eff = ik + (is-1)*num_k_points
        CALL xml_newElement(xf, "K-POINT")
        CALL xml_addAttribute (xf, "Weight", wk(ik_eff) )
        CALL xml_addCharacters (xf, xk(:,ik_eff) )
        CALL xml_endElement(xf, "K-POINT")
        !
        CALL xml_newElement (xf, "E")
        CALL xml_addCharacters (xf, et(:,ik_eff) )
        CALL xml_endElement (xf, "E" )
        !
        CALL xml_newElement (xf, "PROJS")
        DO nwfc = 1, natomwfc
           CALL xml_newElement (xf, "ATOMIC_WFC")
           CALL xml_addAttribute (xf, "index", nwfc )
           CALL xml_addAttribute (xf, "spin", is )
           ! not-so-smart way of copying complex into double
           DO ibnd = 1, nbnd
              proj_tmp(2*ibnd-1) = DBLE( projs(nwfc,ibnd,ik_eff) )
              proj_tmp(2*ibnd  ) =AIMAG( projs(nwfc,ibnd,ik_eff) )
           END DO
           CALL xml_addCharacters (xf, proj_tmp )
           CALL xml_endElement (xf, "ATOMIC_WFC" )
        ENDDO
        CALL xml_endElement (xf, "PROJS" )
        !
     END DO
     !
  ENDDO
  CALL xml_endElement (xf, "EIGENSTATES")
  !
  ! </EIGENSTATE>
  ! <OVERLAPS> if required
  IF ( lwrite_ovp ) THEN
     CALL xml_newElement (xf, "OVERLAPS" )
     DO ik = 1, num_k_points
        DO is = 1, nspin_lsda
           ik_eff = ik + (is-1)*num_k_points
           CALL xml_newElement (xf, "OVPS")
           CALL xml_addAttribute (xf, "dim", natomwfc )
           CALL xml_addAttribute (xf, "spin", is )
           CALL xml_addCharacters (xf, ovps(:,:,ik_eff) )
           CALL xml_endElement (xf, "OVPS" )
        END DO
     END DO
     CALL xml_endElement (xf, "OVERLAPS" )
  END IF
  !
  CALL xml_endElement (xf, "PROJECTIONS")
  !
  ! </PROJECTIONS>
  !
  CALL xml_close ( xf)
  !
END SUBROUTINE write_xml_proj
!
!-----------------------------------------------------------------------
SUBROUTINE write_proj_file ( filproj, proj )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,  ONLY : nspin
  USE spin_orb,  ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin
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
  USE projections, ONLY : nlmchi
  !
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: proj(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256) :: filename
  REAL (DP), EXTERNAL :: compute_mj
  INTEGER :: is, ik, nwfc, ibnd, nksinit, nkslast, iunproj=33
  !
  IF ( TRIM(filproj) == ' ' ) RETURN
  !
  DO is=1,nspin
     IF (nspin==2) THEN
        IF (is==1) filename=trim(filproj)//'.projwfc_up'
        IF (is==2) filename=trim(filproj)//'.projwfc_down'
        nksinit=(nkstot/2)*(is-1)+1
        nkslast=(nkstot/2)*is
     ELSE
        filename=trim(filproj)//'.projwfc_up'
        nksinit=1
        nkslast=nkstot
     ENDIF
     CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
          dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
          celldm, at, gcutm, dual, ecutwfc, nkstot/nspin, nbnd, natomwfc)
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
