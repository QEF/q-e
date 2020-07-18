!
! Copyright (C) 2020 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE read_proj
  !-----------------------------------------------------------------------
CONTAINS
  !
  SUBROUTINE read_xml_proj ( filename, ierr, natomwfc, nbnd, nkstot, &
     nspin, nelec, ef, xk, wk, et, projs, ovps )
  !-----------------------------------------------------------------------
  USE kinds,    ONLY : dp
  USE xmltools, ONLY : xml_openfile, xml_closefile, xmlr_readtag,&
                       xmlr_opentag, xmlr_closetag, get_attr
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(in) :: filename
  INTEGER, INTENT(OUT) :: ierr, natomwfc, nbnd, nkstot, nspin
  REAL(dp), INTENT(OUT) :: nelec, ef
  REAL(dp), INTENT(OUT), ALLOCATABLE :: et(:,:), xk(:,:), wk(:) 
  COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: projs(:,:,:)
  COMPLEX(dp), INTENT(OUT), OPTIONAL, ALLOCATABLE :: ovps(:,:,:)
  !
  COMPLEX(dp), ALLOCATABLE :: proj(:)
  INTEGER :: num_k_points
  INTEGER :: iun, nw, nw_, ik, ik_eff, is, is_
  LOGICAL :: found
  CHARACTER(len=1) :: dummy
  !
  !
  INQUIRE ( FILE = filename, exist=found )
  IF (.NOT. found ) THEN
     ierr = 1
     CALL infomsg ('read_xml_proj', 'xml data file not found')
  END IF
  !
  iun = xml_openfile (filename)
  IF ( iun == -1 ) THEN
     ierr = 2
     CALL infomsg ('read_xml_proj', 'xml data file not readable')
  END IF
  !
  CALL xmlr_opentag ("PROJECTIONS")
  !
  CALL xmlr_readtag ("HEADER", dummy)
  CALL get_attr ("NUMBER_OF_BANDS", nbnd)
  CALL get_attr ("NUMBER_OF_K-POINTS", num_k_points)
  CALL get_attr ("NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL get_attr ("NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL get_attr ("NUMBER_OF_ELECTRONS", nelec)
  CALL get_attr ("FERMI_ENERGY", ef)
  !
  CALL xmlr_opentag ("EIGENSTATES")
  !
  nkstot = num_k_points
  IF ( nspin == 2 ) nkstot = 2*nkstot
  ALLOCATE ( xk(3,nkstot) , wk(nkstot) )
  ALLOCATE ( et(nbnd, nkstot) )
  ALLOCATE ( projs(natomwfc, nbnd, nkstot) )
  ALLOCATE ( proj(nbnd) )
  !
  DO is = 1, nspin
     DO ik = 1, num_k_points
        ik_eff = ik + (is-1)*num_k_points
        CALL xmlr_readtag("K-POINT", xk(:,ik_eff) )
        CALL get_attr ( "Weight", wk(ik_eff) )
        CALL xmlr_readtag ( "E", et(:,ik_eff) )
        !
        CALL xmlr_opentag ("PROJS")
        DO nw = 1, natomwfc
           ! see comment in write_proj when writing this tag
           ! CALL xmlr_readtag ("ATOMIC_WFC", projs(nw,:,ik_eff) )
           CALL xmlr_readtag ("ATOMIC_WFC", proj)
           projs(nw,:,ik_eff) = proj(:)
           CALL get_attr ( "index", nw_ )
           CALL get_attr ( "spin", is_ )
        ENDDO
        CALL xmlr_closetag ( )
        !
     END DO
     !
  ENDDO
  CALL xmlr_closetag ( )
  DEALLOCATE (proj)
  !
  ! </EIGENSTATES>
  !
  IF ( PRESENT(ovps) ) THEN
     ALLOCATE ( ovps(natomwfc,natomwfc,nkstot) )
     CALL xmlr_opentag("OVERLAPS" )
     DO ik = 1, num_k_points
        DO is = 1, nspin
           ik_eff = ik + (is-1)*num_k_points
           CALL xmlr_readtag ("OVPS", ovps(:,:,ik_eff) )
           CALL get_attr ("dim", nw_ )
           CALL get_attr ("spin",is_ )
        END DO
     END DO
     CALL xmlr_closetag ( )
  END IF
  !
  CALL xmlr_closetag ( )
  !
  ! </PROJECTIONS>
  !
  CALL xml_closefile ( )
  !
  ierr = 0
  !
END subroutine read_xml_proj

END MODULE read_proj
