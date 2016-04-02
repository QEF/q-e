!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Calculation of Grimme D2 contribution to the dyamical matrix
! See module "london_module" in Modules/mm_dispersion.f90
! Written by Fabrizio Masullo for his M.Sc. in Mathematic at UniUD
! under the supervision of Paolo Giannozzi
!------------------------------------------------------------------------------
!
!
SUBROUTINE d2ionq_mm ( alat , nat , ityp , at , bg , tau, q, deriv2_london )
  !
  USE london_module
  USE constants,    ONLY : tpi, eps8
  USE mp_images,    ONLY : me_image , nproc_image , intra_image_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER , INTENT ( in ) :: nat , ityp ( nat )
  ! input:
  ! nat  : number of atoms
  ! ityp : type of each atom
  !
  REAL ( DP ) , INTENT ( in ) :: alat, tau(3, nat),  at(3 , 3), bg(3 , 3), q(3)
  ! input:
  ! alat : the cell parameter
  ! tau  : atomic positions in alat units
  ! at   : direct lattice vectors
  ! bg   : reciprocal lattice vectors
  ! q    : wave-vector (in 2pi/alat units)
  !
  COMPLEX ( DP ), INTENT(OUT) :: deriv2_london ( 3, nat, 3, nat )
  !
  INTEGER :: ata , atb , nrm , nr , ipol, jpol
  ! locals :
  ! ata , atb : atom counters
  ! nrm       : actual number of vectors computed by rgen
  ! nr        : counter on neighbours shells
  ! ipol, jpol: counters on coords
  !
  INTEGER :: first , last , resto, divid
  ! locals :
  ! first  : lower bound on processor
  ! last   : upper
  !
  REAL ( DP ) :: dist , f_damp , dtau ( 3 ) , &
       exparg , expval, par , par2, fac , facF, add, addF, auxr
  COMPLEX ( DP ) :: eiqr
  ! locals :
  ! dist         : distance R_ij between the current pair of atoms
  ! f_damp       :  damping function
  ! dtau         :  \vec R_ij
  ! ... and many other temporary variables, plus some buffers:
  !
  REAL ( DP ) ::    aux (3, 3, nat)
  COMPLEX ( DP ) :: aux2(3, 3, nat)
  !
  !
  deriv2_london ( : , : , : , :) = 0.d0
  !
#if defined __MPI
  !
  ! parallelization: divide atoms across processors of this image
  ! (different images have different atomic positions)
  !
  resto = mod ( nat , nproc_image )
  divid = nat / nproc_image
  !
  IF ( me_image + 1 <= resto ) THEN
     !
     first = ( divid  + 1 ) * me_image + 1
     last  = ( divid  + 1 ) * ( me_image + 1 )
     !
  ELSE
     !
     first = ( ( divid + 1 ) * resto ) + ( divid ) * ( me_image-resto ) + 1
     last  = ( divid  + 1 ) * resto + ( divid ) * ( me_image - resto + 1 )
     !
  ENDIF
  !
#else
  !
  first = 1
  last  = nat
#endif
  !
  DO ata = first , last
     !
     aux(:,:,:) = 0.d0
     aux2(:,:,:) = 0.d0
     !
     DO atb = 1 , nat
        !
        dtau ( : ) = tau ( : , ata ) - tau ( : , atb )
        !
        ! generate neighbours shells
        
        CALL rgen ( dtau, r_cut, mxr, at, bg, r, dist2, nrm )
        !
        ! compute forces
        !
        par = beta / ( R_sum ( ityp ( atb ) , ityp ( ata ) ) )
        !
        par2 = par**2
        !
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1600)
!$omp parallel do private(nr,dist,exparg,expval,fac,add,eiqr,facF,addF,auxr,ipol,jpol) default(shared), reduction(+:aux), reduction(+:aux2)
#endif
        DO nr = 1 , nrm
           !
           dist  = alat * sqrt ( dist2 ( nr ) )
           IF ( dist > eps8 ) THEN
              !
              exparg = - beta * ( dist / ( R_sum ( ityp(atb) , ityp(ata) ) ) - 1.0_dp )
              expval = exp ( exparg )
              !
              fac = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist**8
              add = 48.d0 / dist**2
              !
              eiqr = exp ((0_dp,1_dp)*(q(1)*(r(1,nr)+dtau(1))+&
                                       q(2)*(r(2,nr)+dtau(2))+&
                                       q(3)*(r(3,nr)+dtau(3)) ) * tpi )
              !
              facF = C6_ij ( ityp ( atb ) , ityp ( ata ) ) / dist**6
              addF = 6.d0 / dist
              !
              DO ipol = 1 , 3
                 DO jpol = 1 , 3
                    IF (ipol /= jpol) THEN
                       auxr = ( scal6 / ( 1.d0 + expval ) * fac * &
                            ( add - par*(13.d0*expval/((1.d0+expval)*dist))-&
                            par2 *expval*(1.0_dp - expval)/ &
                            ( (1.d0 + expval) * (1.d0 + expval) ) ) * &
                            r ( ipol, nr ) * alat * r(jpol,nr) * alat )
                       
                    ELSE
                       auxr = ( scal6 / ( 1.d0 + expval ) * fac * &
                            ( add - par*(13.d0*expval/((1.d0+expval)*dist))-&
                            par2 *expval*(1.0_dp - expval)/ &
                            ( (1.d0 + expval) * (1.d0 + expval) ) ) * &
                            r ( ipol, nr ) * alat * r(jpol,nr) * alat )- &
                            ( scal6 / ( 1.0_dp + expval ) * facF * &
                            ( - par * expval / ( 1.d0 + expval ) + addF ) * &
                            1.d0/ dist )
                    ENDIF
                    !
                    aux (ipol,jpol,atb) = aux (ipol,jpol,atb) + auxr
                    aux2(ipol,jpol,atb) = aux2(ipol,jpol,atb) + auxr*eiqr
                    !
                 ENDDO
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER < 1600)
!$omp end parallel do
#endif
        DO ipol =1,3
           DO jpol = 1,3
              deriv2_london (ipol, ata, jpol, atb) = aux2(ipol,jpol,atb)
           ENDDO
        ENDDO
     ENDDO
     !
     DO atb = 1, nat
        DO ipol =1,3
           DO jpol = 1,3
              deriv2_london (ipol, ata, jpol, ata) = &
                   deriv2_london (ipol, ata, jpol, ata) - aux(ipol,jpol,atb)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO
  
  CALL mp_sum ( deriv2_london , intra_image_comm )
  !
  RETURN
  !
END SUBROUTINE d2ionq_mm
