!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   module non_local_core_correction
!=----------------------------------------------------------------------------=!

     USE kinds

     IMPLICIT NONE
     SAVE

     PRIVATE

     PUBLIC :: add_core_charge,  core_charge_forces

!=----------------------------------------------------------------------------=!
   contains
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   subroutine add_core_charge( rhoetg, rhoetr, sfac, ps, gv, nsp)
!=----------------------------------------------------------------------------=!

     USE fft, ONLY: pinvfft
     USE cp_types, ONLY: pseudo, recvecs

     implicit none

     integer :: nsp
     COMPLEX(dbl) :: rhoetg(:)
     REAL(dbl)    :: rhoetr(:,:,:)
     type (pseudo),  intent(in) :: ps
     type (recvecs), intent(in) :: gv
     COMPLEX(dbl), INTENT(IN) :: sfac(:,:)
          
     COMPLEX(dbl), ALLOCATABLE :: vtemp(:)
     integer :: is, ig, ng

     ng = gv%ng_l

     ALLOCATE( vtemp( ng ) )
     vtemp = CMPLX( 0.0d0 )

     DO is = 1, nsp
       if( ps%tnlcc( is ) ) then
         do ig = 1, ng
           vtemp(ig) = vtemp(ig) + sfac( is, ig ) * CMPLX(ps%rhoc1(ig,is),0.0d0)
         end do
       endif
     end do

     rhoetg( 1:ng ) = rhoetg( 1:ng ) + vtemp( 1:ng )

     !  rhoetr = 1.0d0 * rhoetr + INVFFT( vtemp )
     CALL pinvfft( rhoetr, vtemp, 1.0d0 )

     DEALLOCATE( vtemp )

     RETURN
!=----------------------------------------------------------------------------=!
   end subroutine add_core_charge
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   SUBROUTINE core_charge_forces( fion, vxc, rhoc1, tnlcc, atoms, ht, eigr, gv, kp )
!=----------------------------------------------------------------------------=!

     !   This subroutine computes the non local core correction
     !   contribution to the atomic forces

     USE constants
     USE cell_base, ONLY: tpiba
     USE cell_module, ONLY: boxdimensions
     USE brillouin, ONLY: kpoints
     USE atoms_type_module, ONLY: atoms_type
     USE cp_types, ONLY: phase_factors, recvecs

     IMPLICIT NONE

     TYPE (atoms_type), INTENT(IN) :: atoms    !   atomic positions
     TYPE (recvecs), INTENT(IN) :: gv          !   reciprocal space vectors
     TYPE (boxdimensions), INTENT(IN) :: ht    !   cell parameters
     TYPE (kpoints), INTENT(IN) :: kp          !   k points
     TYPE (phase_factors), INTENT(IN) :: eigr  !   phase factors
     LOGICAL      :: tnlcc(:)                  !   NLCC flags
     REAL(dbl)    :: fion(:,:)                 !   ionic forces
     REAL(dbl)    :: rhoc1(:,:)                !   derivative of the core charge
     COMPLEX(dbl) :: vxc(:,:)                  !   XC potential

     INTEGER :: ig, ig1, ig2, ig3, isa, ia, is, ispin, nspin
     COMPLEX(dbl) :: gx, gy, gz, tx, ty, tz, teigr, cxc
     COMPLEX(dbl), ALLOCATABLE :: ftmp(:,:)
     REAL(dbl) :: cost

     IF( ANY( tnlcc ) ) then

       nspin = SIZE( vxc, 2)
       ALLOCATE( ftmp( 3, atoms%nat ) )
       ftmp = CMPLX( 0.0d0, 0.0d0 )

       DO ig = gv%gstart, gv%ng_l
         ig1  = gv%mill(1,ig)
         ig2  = gv%mill(2,ig)
         ig3  = gv%mill(3,ig)
         GX   = CMPLX(0.D0, gv%gx_l(1,ig))
         GY   = CMPLX(0.D0, gv%gx_l(2,ig))
         GZ   = CMPLX(0.D0, gv%gx_l(3,ig))
         isa = 1
         DO is = 1, atoms%nsp
           IF ( tnlcc(is) ) THEN
             CXC = 0.0_dbl
             DO ispin = 1, nspin
               CXC = CXC + rhoc1( ig, is ) * CONJG( vxc( ig, ispin ) )
             END DO
             TX = CXC * GX
             TY = CXC * GY
             TZ = CXC * GZ
             DO IA = 1, atoms%na(is)
               teigr = eigr%x( ig1, isa ) * eigr%y( ig2, isa ) * eigr%z( ig3, isa )
               ftmp( 1, isa ) = ftmp( 1, isa ) + TEIGR * TX
               ftmp( 2, isa ) = ftmp( 2, isa ) + TEIGR * TY
               ftmp( 3, isa ) = ftmp( 3, isa ) + TEIGR * TZ
               isa = isa + 1
             END DO
           ELSE
             isa = isa + atoms%na(is)
           END IF
         END DO
       END DO

       ! ...  each processor add its own contribution to the array FION
       IF(kp%gamma_only) THEN
         cost = 2.D0 * ht%deth * tpiba
       ELSE
         cost =        ht%deth * tpiba
       END IF

       DO isa = 1, atoms%nat
         FION(1,ISA) = FION(1,ISA) + REAL(ftmp(1,ISA)) * cost
         FION(2,ISA) = FION(2,ISA) + REAL(ftmp(2,ISA)) * cost
         FION(3,ISA) = FION(3,ISA) + REAL(ftmp(3,ISA)) * cost
       END DO

       DEALLOCATE( ftmp )

     END IF

     RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   end module non_local_core_correction
!=----------------------------------------------------------------------------=!
