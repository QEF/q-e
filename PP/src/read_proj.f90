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
  SUBROUTINE read_xml_proj ( filename, ierr, natomwfc, nbnd, num_k_points, &
     nspin, nelec, ef, xk, wk, et, projs, ovps )
  !-----------------------------------------------------------------------
  USE kinds,    ONLY : dp
  USE FoX_dom
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(in) :: filename
  INTEGER, INTENT(OUT) :: ierr, natomwfc, nbnd, num_k_points, nspin
  REAL(dp), INTENT(OUT) :: nelec, ef
  REAL(dp), INTENT(OUT), ALLOCATABLE :: et(:,:), xk(:,:), wk(:) 
  COMPLEX(dp), INTENT(OUT), ALLOCATABLE :: projs(:,:,:)
  COMPLEX(dp), INTENT(OUT), OPTIONAL, ALLOCATABLE :: ovps(:,:,:)
  !
  TYPE (Node), POINTER :: root
  TYPE (Node), POINTER :: node1
  TYPE (Node), POINTER :: node2
  TYPE (Node), POINTER :: node3
  TYPE (Node), POINTER :: node4
  TYPE (NodeList), POINTER :: node3list, node4list
  !
  INTEGER :: nw, nwfc, ik, is, nks
  LOGICAL :: found
  !
  !
  INQUIRE ( FILE = filename, exist=found )
  IF (.NOT. found ) &
       CALL errore ('read_xml_proj', 'xml data file not found', 1)
  !
  ! read XML file into "root" object
  !
  root => parseFile ( filename )
  !
  ierr = 1
  node1 => item( getElementsByTagname( root, "PROJECTIONS"), 0)
  IF ( .NOT.ASSOCIATED( node1)) RETURN
  !
  ierr = 10
  node2 => item( getElementsByTagname( node1, "HEADER"), 0)
  IF ( .NOT.ASSOCIATED( node2)) RETURN
  !
  IF ( .NOT. hasAttribute( node2, "NUMBER_OF_BANDS") ) RETURN
  CALL extractDataAttribute( node2, "NUMBER_OF_BANDS", nbnd)
  !
  IF ( .NOT. hasAttribute( node2, "NUMBER_OF_K-POINTS") ) RETURN
  CALL extractDataAttribute( node2, "NUMBER_OF_K-POINTS", num_k_points)
  !
  IF ( .NOT. hasAttribute( node2, "NUMBER_OF_SPIN_COMPONENTS") ) RETURN
  CALL extractDataAttribute( node2, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  !
  IF ( .NOT. hasAttribute( node2, "NUMBER_OF_ATOMIC_WFC") ) RETURN
  CALL extractDataAttribute( node2, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  !
  IF ( .NOT. hasAttribute( node2, "NUMBER_OF_ELECTRONS") ) RETURN
  CALL extractDataAttribute( node2, "NUMBER_OF_ELECTRONS", nelec)
  !
  IF ( .NOT. hasAttribute( node2, "FERMI_ENERGY") ) RETURN
  CALL extractDataAttribute( node2, "FERMI_ENERGY", ef)
  !
  ierr = 2
  node2 => item( getElementsByTagname( node1, "EIGENSTATES"), 0)
  IF ( .NOT. ASSOCIATED( node2) ) RETURN
  !
  IF ( nspin == 2 ) num_k_Points = num_k_points*nspin
  ALLOCATE ( xk(3,num_k_points) , wk(num_k_points) )
  ALLOCATE ( et(nbnd, num_k_points) )
  ALLOCATE ( projs(natomwfc, nbnd, num_k_points) )
  !
  node3list => getElementsByTagname( node2, "K-POINT")
  IF ( getLength (node3list) /= num_k_points  ) RETURN
  !
  DO ik=1,num_k_points
     !
     ierr = 20+ik
     node3 => item( node3list, ik-1)
     IF ( .NOT. ASSOCIATED( node3) ) RETURN
     CALL extractDataContent( node3, xk(:,ik))
     IF ( .NOT. hasAttribute( node3, "Weight") ) RETURN
     !
     CALL extractDataAttribute( node3, "Weight", wk(ik))
     !
  END DO
  !
  node3list => getElementsByTagname( node2, "E")
  IF ( getLength (node3list) /= num_k_points  ) RETURN
  !
  DO ik=1,num_k_points
     !
     ierr = 20+ik
     node3 => item( node3list, ik-1)
     IF ( .NOT. ASSOCIATED( node3) ) RETURN
     !
     CALL extractDataContent( node3, et(:,ik))
     !
  END DO
  !
  node3list => getElementsByTagname( node2, "PROJS")
  IF ( getLength (node3list) /= num_k_points ) RETURN
  !
  DO ik=1,num_k_points
     !
     node3 => item( node3list, ik-1)
     IF ( .NOT. ASSOCIATED( node3) ) RETURN
     node4list => getElementsByTagname( node3, "ATOMIC_WFC")
     IF ( getLength (node4list) /= natomwfc ) RETURN
     !
     DO nw = 1, natomwfc
        node4 => item( node4list, nw-1)
        IF ( .NOT. hasAttribute( node4, "index") ) RETURN
        CALL extractDataAttribute( node4, "index", nwfc)
        IF ( .NOT. hasAttribute( node4, "spin") ) RETURN
        CALL extractDataAttribute( node4, "spin", is)
        IF ( nwfc /= nw ) RETURN
        CALL extractDataContent( node4, projs(nw,:,ik))
     END DO
     !
  END DO
  !
  IF ( PRESENT(ovps) ) THEN
     !
     node2 => item( getElementsByTagname( node1, "OVERLAPS"), 0)
     IF ( ASSOCIATED( node2) ) THEN
        !
        ALLOCATE ( ovps(natomwfc,natomwfc,num_k_points) )
        !
        ierr = 3
        node3list => getElementsByTagname( node2, "OVPS")
        IF ( getLength (node3list) /= num_k_points ) RETURN
        !
        DO ik=1,num_k_points
           !
           node3 => item( node3list, ik-1)
           IF ( .NOT. ASSOCIATED( node3) ) RETURN
           !
           IF ( .NOT. hasAttribute( node3, "dim") ) RETURN
           CALL extractDataAttribute( node3, "dim", nwfc)
           IF ( nwfc /= natomwfc ) RETURN
           !
           IF ( .NOT. hasAttribute( node3, "spin") ) RETURN
           CALL extractDataAttribute( node3, "spin", is)
           !
           CALL extractDataContent( node3, ovps(:,:,ik))
           !
        END DO
        !
     END IF
     !
  END IF
  ierr = 0 
  CALL destroy (root) 
  !
END subroutine read_xml_proj

END MODULE read_proj
