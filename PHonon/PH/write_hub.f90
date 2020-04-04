!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE write_hub
!-----------------------------------------------------------------------
!
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE write_dnsscf_ph
  !---------------------------------------------------------------------
  !
  ! DFPT+U: This routine transforms dnsscf_all_modes 
  ! (which is the response of occupation matrices due to
  ! atomic displacements) from the pattern basis u to cartesian 
  ! coordinates and then writes it to the standard output.
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat, ityp
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, is_hubbard
  USE ldaU_ph,         ONLY : dnsscf_all_modes
  USE lsda_mod,        ONLY : nspin
  USE modes,           ONLY : u, nmodes
  !
  IMPLICIT none
  !
  INTEGER :: na_icart, nah, is, m1, m2, na, icart, nt, &
             nap_jcar, nap, na_icar, jcar, icar, nb, imode  
  COMPLEX(DP), ALLOCATABLE :: dnsscf_all_modes_cart(:,:,:,:,:)
  !
  ALLOCATE (dnsscf_all_modes_cart(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,nmodes))
  dnsscf_all_modes_cart = (0.d0, 0.d0)
  !
  ! Transform dnsscf_all_modes from the pattern basis u to cartesian coordinates
  !    
  DO na_icart = 1, 3*nat
     DO imode = 1, nmodes
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 DO m1 = 1, 2 * Hubbard_l(nt) + 1
                    DO m2 = 1, 2 * Hubbard_l(nt) + 1
                       dnsscf_all_modes_cart(m1,m2,is,nah,na_icart) = &
                       dnsscf_all_modes_cart(m1,m2,is,nah,na_icart) + &
                       dnsscf_all_modes(m1,m2,is,nah,imode) * CONJG(u(na_icart,imode))
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  ! Write dnsscf_all_modes_cart to the standard output
  !
  WRITE(stdout,*)
  WRITE(stdout,*) 'DNS_SCF SYMMETRIZED IN CARTESIAN COORDINATES'
  DO na = 1, nat
     DO icart = 1, 3
        WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
        na_icart = 3 * (na - 1) + icart
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                 DO m1 = 1, 2 * Hubbard_l(nt) + 1
                    WRITE( stdout,'(14(f15.10,1x))') dnsscf_all_modes_cart (m1,:,is,nah,na_icart)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  WRITE(stdout,*)
  !
  DEALLOCATE(dnsscf_all_modes_cart)
  !
  RETURN
  !
END SUBROUTINE write_dnsscf_ph
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE write_dnsscf_e 
  !---------------------------------------------------------------------
  ! 
  ! DFPT+U: This routine transforms dnsscf_all_modes 
  ! (which is the response of occupation matrices due to
  ! the electric field perturbation) from crystal to cartesian 
  ! coordinates and then writes it to the standard output.
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat, ityp
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, is_hubbard   
  USE ldaU_ph,         ONLY : dnsscf_all_modes 
  USE cell_base,       ONLY : at
  USE lsda_mod,        ONLY : nspin
  !
  IMPLICIT none
  !  
  INTEGER :: icart, ipol, nt, nah, m1, m2, is
  COMPLEX(DP), ALLOCATABLE :: dnsscf_all_modes_cart(:,:,:,:,:)
  !
  ALLOCATE (dnsscf_all_modes_cart(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3))
  dnsscf_all_modes_cart = (0.d0, 0.d0)
  !
  ! Transform dnsscf_all_modes from crystal to cartesian coordinates
  !
  DO icart = 1, 3 
     DO ipol = 1, 3
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       dnsscf_all_modes_cart(m1,m2,is,nah,icart) = &
                       dnsscf_all_modes_cart(m1,m2,is,nah,icart) + &
                           dnsscf_all_modes(m1,m2,is,nah,ipol)* at(icart,ipol) 
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  ! Write dnsscf_all_modes_cart to the standard output
  !
  WRITE(stdout,*)
  WRITE(stdout,*) 'DNS_SCF SYMMETRIZED IN ELECTRIC FIELD IN CARTESIAN COORDINATES'
  DO icart = 1, 3
     WRITE(stdout,'(a,1x,i2)') 'icart=', icart
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           DO is = 1, nspin
              WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
              DO m1 = 1, 2*Hubbard_l(nt)+1
                 WRITE(stdout,'(14(f15.10,1x))') dnsscf_all_modes_cart (m1,:,is,nah,icart) 
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  WRITE(stdout,*)
  !
  DEALLOCATE(dnsscf_all_modes_cart)
  !
  RETURN
  !
END SUBROUTINE write_dnsscf_e
!-------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE write_dynmat_hub 
  !---------------------------------------------------------------------
  ! 
  ! DFPT+U: This routine writes the scf and total hubbard dynamical 
  ! matrix to the standard output. 
  !
  ! Written  by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  ! 
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode
  USE dynmat,        ONLY : dyn_hub_scf, dyn_hub_bare
  USE ldaU_ph,       ONLY : dnsscf_all_modes      
  USE ldaU,          ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, is_hubbard
  USE control_flags, ONLY : iverbosity 
  USE lsda_mod,      ONLY : nspin
  USE ions_base,     ONLY : nat, ityp
  USE modes,         ONLY : u, nmodes
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: dyn_hub_tot(:,:)
  !
  ALLOCATE(dyn_hub_tot(3*nat,3*nat))
  dyn_hub_tot = (0.d0, 0.d0)
  !
  ! Write the UNSYMMETRIZED SCF Hubbard dynamical matrix 
  ! in the pattern basis
  ! 
  CALL tra_write_matrix_no_sym ('dyn_hub_scf NOT SYMMETRIZED',dyn_hub_scf,nat)
  !
  ! Write the SYMMETRIZED SCF Hubbard dynamical matrix 
  ! in carthesian coordinates
  !
  CALL tra_write_matrix ('dyn_hub_scf SYMMETRIZED',dyn_hub_scf,u,nat)
  !
  ! The total Hubbard dynamical matrix
  !
  IF (ALLOCATED(dyn_hub_bare)) THEN
     dyn_hub_tot = dyn_hub_scf + dyn_hub_bare 
  ELSE
     WRITE(stdout,'("Warning! dyn_hub_bare is not allocated.")')
  ENDIF
  !
  ! Write the UNSYMMETRIZED total Hubbard dynamical matrix 
  ! in the pattern basis
  !
  CALL tra_write_matrix_no_sym('dyn_hub_tot NOT SYMMETRIZED',dyn_hub_tot,nat)
  !
  ! Write the SYMMETRIZED total Hubbard dynamical matrix 
  ! in carthesian coordinates
  !
  CALL tra_write_matrix('dyn_hub_tot SYMMETRIZED',dyn_hub_tot,u,nat)
  !
  DEALLOCATE (dyn_hub_tot)
  !
  RETURN 
  ! 
END SUBROUTINE write_dynmat_hub
!---------------------------------------------------------------------------

END MODULE write_hub
!-----------------------------------------------------------------------
