  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE wan2bloch
  !----------------------------------------------------------------------
  !! Contains all the routines that transforms quantities from Wannier 
  !! space to Bloch space. 
  !! 
  !!  
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE hamwan2bloch ( nbnd, nrr, cuf, eig, chw, cfac)
    !--------------------------------------------------------------------------
    !
    !  From the Hamiltonian in Wannier representation, find the corresponding
    !  Hamiltonian in Bloch representation for a given k point
    !  
    !  input  : number of bands nbnd
    !           number of WS vectors, coordinates and degeneracy 
    !           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
    !           kpoint coordinate xk(3)
    !
    !  output : rotation matrix cuf (nbnd, nbnd)
    !           interpolated hamiltonian eigenvalues eig(nbnd)
    !
    !  SP [optimization]
    !
    !  Feliciano Giustino, UCB
    !
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, cone
    !
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT (in) :: nrr
    !! number of WS points
    ! 
    REAL(kind=DP), INTENT (out) :: eig (nbnd)
    !! interpolated hamiltonian eigenvalues for this kpoint 
    ! 
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(kind=DP), INTENT (in) :: chw ( nbnd, nbnd, nrr)
    !! Hamiltonian in Bloch basis, fine mesh
    COMPLEX(kind=DP), INTENT (out) :: cuf(nbnd, nbnd)
    !! Rotation matrix, fine mesh
    !
    ! variables for lapack ZHPEVX 
    !
    INTEGER :: neig, info, ifail( nbnd ), iwork( 5*nbnd )
    REAL(kind=DP) :: w( nbnd )
    REAL(kind=DP) :: rwork( 7*nbnd )
    COMPLEX(kind=DP) :: champ( nbnd*(nbnd+1)/2 )
    COMPLEX(kind=DP) :: cwork( 2*nbnd )
    COMPLEX(kind=DP) :: cz( nbnd, nbnd)
    !
    ! work variables 
    !
    INTEGER :: ibnd, jbnd
    COMPLEX(kind=DP) :: chf(nbnd, nbnd) 
    !
    CALL start_clock('HamW2B')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
    !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
    !
    !  H~_k   is chf ( nbnd, nbnd, 2*ik-1 )
    !  H~_k+q is chf ( nbnd, nbnd, 2*ik   )
    !
    chf (:,:) = czero
    !
    ! Previous implementation
    !DO ir = 1, nrr
    !   !
    !   rdotk = twopi * ( xk(1)*irvec(1,ir) + xk(2)*irvec(2,ir) + xk(3)*irvec(3,ir))
    !   cfac = exp( ci*rdotk ) / ndegen(ir)
    !   !
    !   ! SP : Significantly faster !
    !   !chf (:,:) = chf (:,:) + cfac * chw (:,:, ir )
    !   CALL zaxpy(nbnd**2,cfac, chw (1,1, ir ),1, chf (1,1),1)
    !   !
    !ENDDO
    ! New one
    CALL zgemv('n', nbnd**2, nrr, cone, chw, nbnd**2, cfac, 1, cone, chf, 1 )
    !
    !---------------------------------------------------------------------
    !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
    !---------------------------------------------------------------------
    !
    ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
    ! after hermitian-ization
    !
    DO jbnd = 1, nbnd
     DO ibnd = 1, jbnd
        champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
        ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) * 0.5d0
     ENDDO
    ENDDO
    !
    CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
                 0, 0,-1.0, neig, w, cz, nbnd, cwork, &
                 rwork, iwork, ifail, info)
    !
    ! rotation matrix and Ham eigenvalues 
    ! [in Ry, mind when comparing with wannier code]
    ! 
    cuf = conjg( transpose ( cz ) )
    eig = w 
    !
    CALL stop_clock('HamW2B')
    !
    END SUBROUTINE hamwan2bloch
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
    USE constants_epw, ONLY : twopi, czero
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
    INTEGER :: imode, jmode, na, nb, n1,n2,n3, m1,m2,m3, ipol,jpol, i
    REAL(kind=DP) :: xq (3)
    REAL(kind=DP) :: massfac
    REAL(kind=DP), EXTERNAL :: wsweight
    REAL(kind=DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
    REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)
    COMPLEX(kind=DP) :: chf(nmodes, nmodes), dyn(3,3,nat,nat)
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
    USE cell_base, ONLY : at 
    USE phcom,     ONLY : nq1, nq2, nq3
    USE ions_base, ONLY : tau, nat
    USE elph2,     ONLY : ifc, epsi, zstar
    USE epwcom,    ONLY : lpolar
    USE constants_epw, ONLY : twopi, czero
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
    INTEGER :: na, nb, n1,n2,n3, m1,m2,m3, ipol,jpol, i
    REAL(kind=DP), EXTERNAL :: wsweight
    REAL(kind=DP),SAVE,ALLOCATABLE :: wscache(:,:,:,:,:)
    REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)  
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
    !--------------------------------------------------------------------------
    SUBROUTINE dmewan2bloch ( nbnd, nrr, cuf, dmef, etf, etf_ks, cfac)
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Wannier representation, find the corresponding
    !!  Dipole in Bloch representation for a given k point
    !!  
    !!  input  : number of bands nbnd
    !!           number of WS vectors, coordinates and degeneracy 
    !!           Dipole in Wannier representation  cdmew(3,nbnd,nbnd,nrr)
    !!
    !!  output : interpolated dipole matrix elements (dmef)
    !!
    !!  SP 09/2016: Optimization
    !!  JN, EK 09/2010
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : cdmew
    USE epwcom,        ONLY : eig_read
    USE constants_epw, ONLY : cone, czero
    USE mp_world,      ONLY : mpime
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nbnd 
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT (in) :: nrr 
    !! kpoint number for the interpolation
    !
    REAL(kind=DP), INTENT (in) :: etf (nbnd)
    !! Eigenenergies on the fine grid
    REAL(kind=DP), INTENT (in) :: etf_ks (nbnd) 
    !! 
    COMPLEX(kind=DP), INTENT (in) :: cuf (nbnd, nbnd)
    !! Rotation matrix, fine mesh 
    COMPLEX(kind=DP), INTENT (out) :: dmef (3, nbnd, nbnd)
    !! interpolated dipole matrix elements in Bloch basis, fine mesh
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    ! 
    ! local variables
    !
    INTEGER :: ibnd, jbnd, j 
    COMPLEX( kind=DP ) :: congj_cuf (nbnd, nbnd)
    COMPLEX( kind=DP ) :: cdmef (3, nbnd, nbnd)
    COMPLEX( kind=DP ) :: cdmef_tmp (3, nbnd, nbnd)
    ! dipole matrix elements in Bloch basis, fine mesh
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
    !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
    !
    !  H~_k   is chf ( nbnd, nbnd, 2*ik-1 )
    !  H~_k+q is chf ( nbnd, nbnd, 2*ik   )
    !
    ! Initialization
    cdmef_tmp (:,:,:) = czero
    dmef (:,:,: ) = czero
    cdmef (:,:,:) = czero
    !   
    ! SUM on ir of: cdmef (1,ibnd,jbnd) = cdmef (1,ibnd,jbnd) + cfac( ir) * cdmew (1, ibnd,jbnd, ir )
    CALL zgemv('n', 3*(nbnd**2), nrr, cone, cdmew(:,:,:,:), 3*(nbnd**2), cfac(:), 1, cone, cdmef(:,:,:), 1  )
    ! 
    !
    ! pmn(k) = U p~ U^dagger
    ! cuf,  passed from hamwan2bloch.
    !
    congj_cuf = Conjg(Transpose(cuf))
    DO j=1, 3
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cdmef(j,:,:), &
           nbnd, congj_cuf(:,:), nbnd, czero, cdmef_tmp(j,:,:), nbnd)
    ENDDO
    !
    DO j=1, 3
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:,:), &
                 nbnd, cdmef_tmp(j,:,:), nbnd, czero, dmef(j,:,:), nbnd)
    ENDDO
    !
    ! Satisfy
    ! Phys. Rev. B 62, 4927–4944 (2000) , Eq. (30)
    !
    IF (eig_read) THEN
       DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
             IF (abs(etf_ks(ibnd) - etf_ks(jbnd)) .gt. 1.d-4) THEN
                dmef(:, ibnd, jbnd) = dmef(:,ibnd, jbnd) * &
                     ( etf(ibnd)    - etf(jbnd) )/ &
                     ( etf_ks(ibnd) - etf_ks(jbnd) )
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    END SUBROUTINE dmewan2bloch
    !
    !--------------------------------------------------------------------------
    subroutine vmewan2bloch ( nbnd, nrr, irvec, ndegen, xk, cuf, vmef, et, et_ks, chw)
    !--------------------------------------------------------------------------
    !
    !  From the Velocity matrix elements in Wannier representation, find the corresponding
    !  MEs in Bloch representation for a given k point
    !  
    !  input  : nbnd, nrr, irvec, ndegen, xk, cuf, et
    !
    !  output : vmef; velocity matrix elements on the fine mesh
    !
    !  Adapted from hamwan2bloch by Jesse Noffsinger and Emmanouil Kioupakis
    !
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    use elph2,         ONLY : cvmew !, chw
    use cell_base,     ONLY : at, celldm
    USE epwcom,        ONLY : eig_read
    USE constants_epw, ONLY : twopi, ci, czero
    !
    implicit none
    !
    !  input variables
    !
    integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr), ipol
    ! number of bands (possibly of the optimal subspace)
    ! kpoint number for the interpolation
    ! record length and unit for direct write of rotation matrix
    ! number of WS points, crystal coordinates, degeneracy
    !
    ! Hamiltonian in wannier basis
    !
    real(kind=DP) :: xk (3), et(nbnd),et_ks(nbnd),irvec_tmp(3)
    ! kpoint coordinates for the interpolation
    !
    ! output variables
    !  
    complex(kind=DP) :: chf_a(3,nbnd, nbnd), chf_a_tmp(nbnd, nbnd)
    complex(kind=DP) :: vmef (3,nbnd,nbnd)
    ! interpolated hamiltonian eigenvalues for this kpoint 
    complex(kind=DP) :: cuf(nbnd, nbnd), chw(nbnd, nbnd, nrr)
    ! Rotation matrix, fine mesh 
    !
    !
    complex(kind=DP) :: cvmef(3,nbnd, nbnd), cvmef_tmp(nbnd, nbnd)
    ! Hamiltonian in Bloch basis, fine mesh
    integer :: ibnd, jbnd, ir
    real(kind=DP) :: rdotk
    complex(kind=DP) :: cfac
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
    !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
    !
    !  H~_k   is chf ( nbnd, nbnd, 2*ik-1 )
    !  H~_k+q is chf ( nbnd, nbnd, 2*ik   )
    !
    !  
    cvmef (:,:,:) = czero
    !
    DO ir = 1, nrr
       !
       rdotk = twopi * dot_product( xk, dble(irvec( :, ir) ))
       cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
       DO ipol = 1,3
          cvmef (ipol,:,:) = cvmef (ipol,:,:) + cfac * cvmew (ipol, :,:, ir )
       ENDDO
       !
    ENDDO
    !
    !
    ! vmn(k) = U v(amn)~ U^dagger
    !cuf,  passed from hamwan2bloch.
    DO ipol = 1,3
       cvmef_tmp(:,:) = matmul( cvmef(ipol,:,:) ,  conjg(transpose(cuf(:,:)))  )
       vmef(ipol, :,:) = matmul(cuf(:,:) , cvmef_tmp(:,:) )
    ENDDO
    !
    !
    !
    !  get k-derivative of the Hamiltonian in the Wannier gauge
    !
    chf_a (:,:,:) = czero
    !
    DO ir = 1, nrr
       !
       rdotk = twopi * dot_product( xk, dble(irvec( :, ir) ))
       cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
       irvec_tmp(:) = celldm(1) * matmul ( at, dble(irvec(:,ir)) )
       DO ipol = 1, 3
          chf_a (ipol,:,:) = chf_a (ipol,:,:) + &
               ci * irvec_tmp( ipol ) * cfac * chw (:,:, ir )
       ENDDO
       !
    ENDDO
    !
    ! H'mn(k) = U H'~ U^dagger
    ! cuf,  passed from hamwan2bloch.
    DO ipol = 1,3
       chf_a_tmp(:,:) = matmul( chf_a(ipol,:,:) ,  conjg(transpose(cuf(:,:)))  )
       chf_a(ipol, :,:) = matmul(cuf(:,:) , chf_a_tmp(:,:) )
    ENDDO
    !
    !
    DO ibnd = 1, nbnd
       DO jbnd = 1, nbnd
          vmef (:,ibnd,jbnd) = chf_a(:, ibnd, jbnd) - &
               ci * (et_ks(jbnd) - et_ks(ibnd) ) *  vmef(:,ibnd,jbnd)
       ENDDO
    ENDDO
    !
    IF (eig_read) THEN
       DO ibnd = 1, nbnd
          DO jbnd = 1, nbnd
             IF (abs(et_ks(ibnd) - et_ks(jbnd)) .gt. 1.d-4) THEN
                vmef(:, ibnd, jbnd) = vmef(:,ibnd, jbnd) * &
                     ( et(ibnd)    - et(jbnd) )/ &
                     ( et_ks(ibnd) - et_ks(jbnd) )
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    END SUBROUTINE vmewan2bloch
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp ( nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
    !---------------------------------------------------------------------------
    !!
    !! even though this is for phonons, I use the same notations
    !! adopted for the electronic case (nmodes->nmodes etc)
    !!
    USE kinds,         only : DP
    USE epwcom,        only : parallel_k, parallel_q, etf_mem
    USE elph2,         only : epmatwp
    USE constants_epw, ONLY : twopi, ci, czero
    USE io_epw,        ONLY : iunepmatwp, iunepmatwp2
    USE mp_global,     ONLY : mp_sum
    USE mp_world,      ONLY : world_comm
    USE parallel_include
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT (in) :: nrr_q
    !! Number of WS points
    INTEGER, INTENT (in) :: irvec ( 3, nrr_q)
    !! Coordinates of WS points
    INTEGER, INTENT (in) :: ndegen (nrr_q)
    !! Number of degeneracy of WS points
    INTEGER, INTENT (in) :: nbnd
    !! Number of bands
    INTEGER, INTENT (in) ::  nrr_k
    !! Number of electronic WS points
    REAL(kind=DP) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(kind=DP), INTENT (in) :: cuf (nmodes, nmodes)
    !! e-p matrix in Wanner representation
    COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nrr_k, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    ! 
    ! Local variables 
    !
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: ir_start
    !! Starting ir for this cores
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: ierr
    !! Return if there is an error
    !integer(kind=8) ::  lrepmatw,  lrepmatw2
#if defined(__MPI)
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw2
    !! Offset to tell where to start reading the file
#else
    INTEGER(kind=8) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER(kind=8) :: lrepmatw2
    !! Offset to tell where to start reading the file
#endif
    !
    REAL(kind=DP) :: rdotk
    !! Exponential for the FT
    !
    COMPLEX(kind=DP) :: eptmp( nbnd, nbnd, nrr_k, nmodes)
    !! Temporary matrix to store el-ph
    COMPLEX(kind=DP) :: cfac(nrr_q)
    !! Factor for the FT
    COMPLEX(kind=DP), ALLOCATABLE :: epmatw ( :,:,:,:)
    !! El-ph matrix elements
    !
    CALL start_clock('ephW2Bp')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
    !
    !  g~(k') is epmatf (nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    IF (parallel_k) THEN
       CALL para_bounds(ir_start, ir_stop, nrr_q)
    ELSEIF (parallel_q) THEN
       ir_start = 1
       ir_stop  = nrr_q
    ELSE 
       CALL errore ('ephwan2blochp', 'Problem with parallel_k/q scheme', nrr_q)
    ENDIF
    !
    eptmp = czero
    cfac(:) = czero
    !
    DO ir = ir_start, ir_stop
       !   
       ! note xxq is assumed to be already in cryst coord
       !
       rdotk = twopi * dot_product ( xxq, dble(irvec(:, ir)) )
       cfac(ir) = exp( ci*rdotk ) / dble( ndegen(ir) )
    ENDDO
    ! 
    IF (etf_mem == 0) then
      !      
      ! SP: This is faster by 20 % 
      Call zgemv( 'n',  nbnd * nbnd * nrr_k * nmodes, ir_stop - ir_start + 1, ( 1.d0, 0.d0 ),&
               epmatwp(1,1,1,1,ir_start), nbnd * nbnd * nrr_k * nmodes, cfac(ir_start),1,( 0.d0, 0.d0),eptmp, 1 )    
      !
    ELSE
      !
      ALLOCATE(epmatw ( nbnd, nbnd, nrr_k, nmodes))
      !
      !lrepmatw2   = 2 * nbnd * nbnd * nrr_k * nmodes
#if defined(__MPI)
      ! Although this should almost never be problematic (see explaination below)
      lrepmatw2 = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                      INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                      INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
                                      INT( nmodes, kind = MPI_OFFSET_KIND )
#else
      lrepmatw2 = INT( 2 * nbnd * nbnd * nrr_k * nmodes, kind = 8)
#endif
      ! 
      DO ir = ir_start, ir_stop
#if defined(__MPI)
        ! DEBUG: print*,'Process ',my_id,' do ',ir,'/ ',ir_stop
        !
        !  Direct read of epmatwp for this ir
        !lrepmatw   = 2 * nbnd * nbnd * nrr_k * nmodes * 8 * (ir-1)
        ! 
        ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or 
        !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
        !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
        lrepmatw = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                       INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                       INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
                                       INT( nmodes, kind = MPI_OFFSET_KIND ) * &
                 8_MPI_OFFSET_KIND * ( INT( ir    , kind = MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND )
        !
        ! SP: mpi seek is used to set the position at which we should start
        ! reading the file. It is given in bits. 
        ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
        !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ. 
        !        Here we want non blocking because not all the process have the same nb of ir. 
        !
        CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
        IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
        !CALL MPI_FILE_READ(iunepmatwp2, aux, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
        CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
        IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
        ! 
#else      
        call rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
#endif
        !
        CALL ZAXPY(nbnd * nbnd * nrr_k * nmodes, cfac(ir), epmatw, 1, eptmp, 1)
        ! 
      ENDDO
      DEALLOCATE(epmatw)
    ENDIF
    !
#if defined(__MPI)
    IF (parallel_k) CALL mp_sum(eptmp, world_comm)
#endif  
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
    !
    Call zgemm( 'n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, ( 1.d0, 0.d0 ),eptmp  , nbnd * nbnd * nrr_k, &
                                                                                      cuf, nmodes         , &
                                                               ( 0.d0, 0.d0 ),epmatf, nbnd * nbnd * nrr_k )
    !
    CALL stop_clock('ephW2Bp')
    !
    END SUBROUTINE ephwan2blochp
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch ( nbnd, nrr, irvec, ndegen, epmatw, &
           xk, cufkk, cufkq, epmatf, nmodes)
    !---------------------------------------------------------------------------
    !!
    !! Interpolation from Wannier to the fine Bloch grid of the electron-phonon 
    !! matrix elements
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : twopi, ci, czero, cone
    implicit none
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    INTEGER, INTENT (in) :: nrr
    !! Number of Wigner-Size points
    INTEGER, INTENT (in) :: irvec ( 3, nrr)
    !! Coordinates of WS points
    INTEGER, INTENT (in) :: ndegen (nrr)
    !! Degeneracy of WS points
    INTEGER, INTENT (in) :: nmodes
    !! number of phonon modes
    !
    REAL(kind=DP), INTENT (in) :: xk(3)
    !! kpoint for the interpolation (WARNING: this must be in crystal coord!)
    !
    COMPLEX(kind=DP), INTENT (in) :: epmatw ( nbnd, nbnd, nrr, nmodes)
    !! e-p matrix in Wannier representation
    COMPLEX(kind=DP), INTENT (in) :: cufkk (nbnd, nbnd)
    !! rotation matrix U(k)
    COMPLEX(kind=DP), INTENT (in) :: cufkq (nbnd, nbnd)
    !! rotation matrix U(k+q)
    COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! work variables 
    integer :: ir, imode
    real(kind=DP) :: rdotk
    complex(kind=DP) :: cfac, eptmp( nbnd, nbnd)
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
    !
    !  g~(k') is epmatf (nbnd, nbnd, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    epmatf = czero
    !
    DO ir = 1, nrr
       !
       ! note xk is assumed to be already in cryst coord
       !
       rdotk = twopi * dot_product ( xk, dble(irvec(:, ir)) )
       cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
       !
       DO imode = 1, nmodes
         epmatf (:, :, imode) = epmatf (:, :, imode) + cfac * epmatw ( :, :, ir, imode)
       ENDDO
       !
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  g (k') = U_q^\dagger (k') g~ (k') U_k (k')
    !
    ! the two zgemm calls perform the following ops:
    !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
    !
    DO imode = 1, nmodes
      !
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
           nbnd, epmatf (:,:,imode), nbnd, czero, eptmp, nbnd)
      CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, eptmp, &
           nbnd, cufkk, nbnd, czero, epmatf(:,:,imode), nbnd)
      !
    ENDDO
    !
    !
    END SUBROUTINE ephwan2bloch
    ! 
  END MODULE wan2bloch
