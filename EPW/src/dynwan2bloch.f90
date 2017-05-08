  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynwan2bloch ( nmodes, nrr, irvec, ndegen, xxq, cuf, eig)
  !--------------------------------------------------------------------------
  !!
  !!
  !!  WARNING: this SUBROUTINE is identical to hamwan2bloch.f90, except
  !!           that here rdw is a real array, not a complex one. This is
  !!           required to obtain proper phonon dispersion interpolation
  !!           and corresponds to the reality of the interatomic force
  !!           constants
  !!
  !! -------------------------------------------------------------------------
  !!
  !!  From the Hamiltonian in Wannier representation, find the corresponding
  !!  Hamiltonian in Bloch representation for a given k point
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE phcom,     ONLY : nq1, nq2, nq3
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : rdw, epsi, zstar
  USE epwcom,    ONLY : lpolar
  USE constants_epw, ONLY : twopi, ci, czero
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nmodes
  !! number of modes (possibly of the optimal subspace)
  INTEGER, INTENT (in) :: nrr
  !! kpoint number for the interpolation
  INTEGER, INTENT (in) :: irvec (3, nrr)
  !! record length and unit for direct write of rotation matrix
  INTEGER, INTENT (in) :: ndegen (nrr)
  !! number of WS points, crystal coordinates, degeneracy
  !
  REAL(kind=DP), INTENT (in) :: xxq (3)
  !! kpoint coordinates for the interpolation
  REAL(kind=DP), INTENT (out) :: eig (nmodes)
  !! interpolated hamiltonian eigenvalues for this kpoint
  COMPLEX(kind=DP), INTENT (out) :: cuf(nmodes, nmodes)
  !! Rotation matrix, fine mesh 
  !
  ! work variables
  ! variables for lapack ZHPEVX
  integer :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
  real(kind=DP) :: w( nmodes ), rwork( 7*nmodes )
  complex(kind=DP) :: champ( nmodes*(nmodes+1)/2 ), &
    cwork( 2*nmodes ), cz( nmodes, nmodes)
  !
  real(kind=DP) :: xq (3)
  complex(kind=DP) :: chf(nmodes, nmodes)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: imode, jmode, ir, na, nb
  real(kind=DP) :: rdotk, massfac
  complex(kind=DP) :: cfac
  !
  CALL start_clock ( 'DynW2B' )
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform to fine k and k+q meshes
  !----------------------------------------------------------
  !
  !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
  !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
  !
  !  H~_k   is chf ( nmodes, nmodes, 2*ik-1 )
  !  H~_k+q is chf ( nmodes, nmodes, 2*ik   )
  !
  xq = xxq
  chf (:,:) = czero
  !
  DO ir = 1, nrr
    !
    rdotk = twopi * dot_product( xq, dble(irvec( :, ir) ))
    cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
    chf = chf + cfac * rdw (:,:, ir )
    !
  ENDDO
  !
  ! bring xq in cart. coordinates (needed for rgd_blk call)
  CALL cryst_to_cart (1, xq, bg, 1)
  !
  !  add the long-range term to D(q)
  IF (lpolar) THEN
    ! xq has to be in 2pi/a     
    CALL rgd_blk (nq1,nq2,nq3,nat,chf,xq,tau,epsi,zstar,+1.d0)
    !
  ENDIF
  !
  !  divide by the square root of masses 
  !
  DO na = 1, nat
    DO nb = 1, nat
      massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
      !
      chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
         chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
      ! 
    ENDDO
  ENDDO
  !
  ! bring xq back to crystal coordinates
  CALL cryst_to_cart (1, xq, at, -1)
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jmode = 1, nmodes
   DO imode = 1, jmode
      champ (imode + (jmode - 1) * jmode/2 ) = &
      ( chf ( imode, jmode) + conjg ( chf ( jmode, imode) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nmodes, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nmodes, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues
  ! [in Ry, mind when comparing with wannier code]
  !
  cuf = cz
  eig = w
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynwan2bloch
  !--------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------- 
  SUBROUTINE dynifc2blochf ( nmodes, rws, nrws, xxq, cuf, eig)
  !--------------------------------------------------------------------------
  !!
  !!  From the IFCs in the format of q2r, find the corresponding
  !!  dynamical matrix for a given q point (as in matdyn.x) on the fine grid
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE phcom,     ONLY : nq1, nq2, nq3
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : ifc, epsi, zstar
  USE epwcom,    ONLY : lpolar
  USE constants_epw, ONLY : twopi, ci, czero
  USE io_global, ONLY : stdout
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nmodes
  !! number of modes 
  INTEGER, INTENT (in) :: nrws
  !! Number of Wigner-Size real space vectors
  INTEGER, INTENT (in) :: rws(0:3,nrws)
  !! 
  REAL(kind=DP), INTENT (in) :: xxq (3)
  !! qpoint coordinates for the interpolation
  REAL(kind=DP), INTENT (out) :: eig (nmodes)
  !! interpolated phonon eigenvalues for this qpoint
  COMPLEX(kind=DP), INTENT (out) :: cuf(nmodes, nmodes)
  !! Rotation matrix, fine mesh 
  !
  ! work variables
  ! variables for lapack ZHPEVX
  INTEGER :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
  REAL(kind=DP) :: w( nmodes ), rwork( 7*nmodes )
  COMPLEX(kind=DP) :: champ( nmodes*(nmodes+1)/2 )
  COMPLEX(kind=DP) :: cwork( 2*nmodes ), cz( nmodes, nmodes)
  !
  LOGICAL,SAVE :: first=.true.
  INTEGER :: imode, jmode, ir, na, nb, n1,n2,n3, m1,m2,m3, ipol,jpol, i
  REAL(kind=DP) :: xq (3)
  REAL(kind=DP) :: rdotk, massfac
  REAL(kind=DP), EXTERNAL :: wsweight
  REAL(kind=DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
  REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)
  COMPLEX(kind=DP) :: chf(nmodes, nmodes), dyn(3,3,nat,nat)
  ! Dyn mat in Bloch basis, fine mesh
  COMPLEX(kind=DP) :: cfac
  !
  xq = xxq
  ! bring xq in cart. coordinates
  CALL cryst_to_cart (1, xq, bg, 1)
  !
  FIRST_TIME : IF (first) THEN
    first=.false.
    ALLOCATE( wscache(-2*nq3:2*nq3, -2*nq2:2*nq2, -2*nq1:2*nq1, nat,nat) )
    DO na=1, nat
       DO nb=1, nat
          total_weight=0.0d0
          !
          DO n1=-2*nq1,2*nq1
             DO n2=-2*nq2,2*nq2
                DO n3=-2*nq3,2*nq3
                   DO i=1, 3
                      r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                   END DO
                   wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
                ENDDO
             ENDDO
          ENDDO
      ENDDO
    ENDDO
  ENDIF FIRST_TIME
  !
  CALL start_clock ( 'DynW2B' )
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform to fine k and k+q meshes
  !----------------------------------------------------------
  !
  !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
  !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
  !
  !  H~_k   is chf ( nmodes, nmodes, 2*ik-1 )
  !  H~_k+q is chf ( nmodes, nmodes, 2*ik   )
  !
  chf = czero
  dyn = czero
  !
  DO na=1, nat
     DO nb=1, nat
        total_weight=0.0d0
        DO n1=-2*nq1,2*nq1
           DO n2=-2*nq2,2*nq2
              DO n3=-2*nq3,2*nq3
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                 END DO
                 !
                 weight = wscache(n3,n2,n1,nb,na)
                 IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nq1)
                    IF(m1.LE.0) m1=m1+nq1
                    m2 = MOD(n2+1,nq2)
                    IF(m2.LE.0) m2=m2+nq2
                    m3 = MOD(n3+1,nq3)
                    IF(m3.LE.0) m3=m3+nq3
                    !
                    arg = twopi*(xq(1)*r(1) + xq(2)*r(2) + xq(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               ifc(m1,m2,m3,ipol,jpol,na,nb)*CMPLX(COS(arg),-SIN(arg),kind=DP)*weight
                       END DO
                    END DO
                 END IF
                 total_weight=total_weight + weight
              END DO
           END DO
        END DO
        IF (ABS(total_weight-nq1*nq2*nq3).GT.1.0d-8) THEN
           WRITE(stdout,*) total_weight
           CALL errore ('dynifc2bloch','wrong total_weight',1)
        END IF
     END DO
  END DO
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              chf((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  IF (lpolar) THEN
    ! xq has to be in 2pi/a     
    CALL rgd_blk (nq1,nq2,nq3,nat,chf,xq,tau,epsi,zstar,+1.d0)
    !
  ENDIF
  !
  !  divide by the square root of masses 
  !
  DO na = 1, nat
    DO nb = 1, nat
      massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )
      !
      chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
         chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
      ! 
    ENDDO
  ENDDO
  !
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jmode = 1, nmodes
   DO imode = 1, jmode
      champ (imode + (jmode - 1) * jmode/2 ) = &
      ( chf ( imode, jmode) + conjg ( chf ( jmode, imode) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nmodes, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nmodes, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues
  ! [in Ry, mind when comparing with wannier code]
  !
  cuf = cz
  eig = w
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynifc2blochf
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynifc2blochc ( nmodes, rws, nrws, xq, chf)
  !--------------------------------------------------------------------------
  !!
  !!  From the IFCs in the format of q2r, find the corresponding
  !!  dynamical matrix for a given q point (as in matdyn.x) on the coarse grid
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE phcom,     ONLY : nq1, nq2, nq3
  USE ions_base, ONLY : amass, tau, nat, ityp
  USE elph2,     ONLY : ifc, epsi, zstar
  USE epwcom,    ONLY : lpolar
  USE constants_epw, ONLY : twopi, ci, czero
  USE io_global, ONLY : stdout
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nmodes
  !! number of modes 
  INTEGER, INTENT (in) :: nrws
  !! Number of Wigner-Size real space vectors
  INTEGER, INTENT (in) :: rws(0:3,nrws)
  !! Wigner-Seitz radius 
  REAL(kind=DP), INTENT (in) :: xq (3)
  !! qpoint coordinates for the interpolation
  COMPLEX(kind=DP), INTENT (out) :: chf(nmodes, nmodes)
  !! dyn mat (not divided by the masses)
  !
  ! work variables
  !
  ! Dyn mat in Bloch basis, fine mesh
  LOGICAL,SAVE :: first=.true.
  INTEGER :: imode, jmode, ir, na, nb, n1,n2,n3, m1,m2,m3, ipol,jpol, i
  REAL(kind=DP) :: rdotk, massfac
  REAL(kind=DP), EXTERNAL :: wsweight
  REAL(kind=DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
  REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)  
  COMPLEX(kind=DP) :: cfac
  COMPLEX(kind=DP) :: dyn(3,3,nat,nat)
  !
  ! bring xq in cart. coordinates
  !CALL cryst_to_cart (1, xq, bg, 1)
  !
  FIRST_TIME : IF (first) THEN
    first=.false.
    ALLOCATE( wscache(-2*nq3:2*nq3, -2*nq2:2*nq2, -2*nq1:2*nq1, nat,nat) )
    DO na=1, nat
       DO nb=1, nat
          total_weight=0.0d0
          !
          DO n1=-2*nq1,2*nq1
             DO n2=-2*nq2,2*nq2
                DO n3=-2*nq3,2*nq3
                   DO i=1, 3
                      r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                   END DO
                   wscache(n3,n2,n1,nb,na) = wsweight(r_ws,rws,nrws)
                ENDDO
             ENDDO
          ENDDO
      ENDDO
    ENDDO
  ENDIF FIRST_TIME
  !
  chf = czero
  dyn = czero
  !
  DO na=1, nat
     DO nb=1, nat
        total_weight=0.0d0
        DO n1=-2*nq1,2*nq1
           DO n2=-2*nq2,2*nq2
              DO n3=-2*nq3,2*nq3
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                 END DO
                 !
                 weight = wscache(n3,n2,n1,nb,na)
                 IF (weight .GT. 0.0d0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nq1)
                    IF(m1.LE.0) m1=m1+nq1
                    m2 = MOD(n2+1,nq2)
                    IF(m2.LE.0) m2=m2+nq2
                    m3 = MOD(n3+1,nq3)
                    IF(m3.LE.0) m3=m3+nq3
                    !
                    arg = twopi*(xq(1)*r(1) + xq(2)*r(2) + xq(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               ifc(m1,m2,m3,ipol,jpol,na,nb)*CMPLX(COS(arg),-SIN(arg),kind=DP)*weight
                       END DO
                    END DO
                 END IF
                 total_weight=total_weight + weight
              END DO
           END DO
        END DO
        IF (ABS(total_weight-nq1*nq2*nq3).GT.1.0d-8) THEN
           WRITE(stdout,*) total_weight
           CALL errore ('dynifc2bloch','wrong total_weight',1)
        END IF
     END DO
  END DO
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              chf((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  IF (lpolar) THEN
    ! xq has to be in 2pi/a     
    CALL rgd_blk (nq1,nq2,nq3,nat,chf,xq,tau,epsi,zstar,+1.d0)
    !
  ENDIF
  !
  !
  END SUBROUTINE dynifc2blochc
  !-----------------------------------------------------------------------------------------
