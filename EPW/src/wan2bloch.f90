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
  !! 
  !! Modules that contains all the routines that transforms quantities from Wannier 
  !! space to Bloch space. 
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
    !  output : rotation matrix cuf(nbnd, nbnd)
    !           interpolated hamiltonian eigenvalues eig(nbnd)
    !
    !  SP [optimization]
    !
    !  Feliciano Giustino, UCB
    !
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, cone, zero, one, eps12
    USE epwcom,        ONLY : lphase
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
    COMPLEX(kind=DP), INTENT (in) :: chw( nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(kind=DP), INTENT (out) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh
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
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    !
    COMPLEX(kind=DP) :: chf(nbnd, nbnd)
    !! Hamiltonian in Bloch basis, fine mesh
    COMPLEX(KIND=DP) :: zdotu
    !! Dot product between the two phonon eigenvectors.
    !
    CALL start_clock('HamW2B')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 31 of PRB 76, 165108 (2007)]
    !  H~(k')    = 1/ndegen(R) sum_R e^{ik'R     } H(R)
    !  H~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} H(R)
    !  Note H~(k') is H^(W)(k') in PRB 74, 195118 (2006) notations
    !
    !  H~(k')    is chf( nbnd, nbnd, 2*ik-1 )
    !  H~(k'+q') is chf( nbnd, nbnd, 2*ik   )
    !
    chf(:,:) = czero
    !
    ! Previous implementation
    !DO ir = 1, nrr
    !   !
    !   rdotk = twopi * ( xk(1)*irvec(1,ir) + xk(2)*irvec(2,ir) + xk(3)*irvec(3,ir))
    !   cfac = exp( ci*rdotk ) / ndegen(ir)
    !   !
    !   ! SP : Significantly faster !
    !   !chf(:,:) = chf(:,:) + cfac * chw(:,:,ir)
    !   CALL zaxpy(nbnd**2, cfac, chw(1,1,ir), 1, chf(1,1), 1)
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
        champ(ibnd + (jbnd - 1) * jbnd/2 ) = &
        ( chf( ibnd, jbnd) + conjg ( chf( jbnd, ibnd) ) ) * 0.5d0
     ENDDO
    ENDDO
    !
    CALL zhpevx ('V', 'A', 'U', nbnd, champ , zero, zero, &
                 0, 0, -one, neig, w, cz, nbnd, cwork, &
                 rwork, iwork, ifail, info)
    ! clean noise
    DO jbnd = 1, nbnd
      DO ibnd = 1, nbnd
        IF ( ABS( cz(ibnd,jbnd) ) < eps12 ) cz(ibnd,jbnd) = czero
      ENDDO
    ENDDO
    !  
    ! DS - Impose phase 
    IF (lphase) THEN
      DO jbnd = 1, nbnd
        INNER : DO ibnd = 1, nbnd
          IF ( ABS(cz(ibnd, jbnd)) > eps12 ) THEN
            cz(:, jbnd) = cz(:, jbnd) * conjg( cz(ibnd,jbnd) )
            cz(:, jbnd) = cz(:, jbnd)/sqrt(zdotu(nbnd,conjg(cz(:,jbnd)),1,cz(:, jbnd),1) )
            EXIT INNER
          ENDIF
        END DO INNER
      ENDDO
    ENDIF
    ! 
    ! rotation matrix and Ham eigenvalues in Ryd
    ! [mind when comparing with wannier code (eV units)]
    ! 
    ! U^\dagger is cuf(nbnd,nbnd)
    !
    cuf = conjg( transpose ( cz ) )
    eig = w 
    !
    CALL stop_clock('HamW2B')
    !
    END SUBROUTINE hamwan2bloch
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dynwan2bloch ( nmodes, nrr_q, irvec_q, ndegen_q, xxq, cuf, eig)
    !--------------------------------------------------------------------------
    !!
    !!
    !!  WARNING: this subroutine is identical to hamwan2bloch.f90, except
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
    USE epwcom,    ONLY : lpolar, lphase
    USE constants_epw, ONLY : twopi, ci, czero, zero, one, eps12
    USE rigid,     ONLY : cdiagh2
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nmodes
    !! number of modes (possibly of the optimal subspace)
    INTEGER, INTENT (in) :: nrr_q
    !! number of WS points
    INTEGER, INTENT (in) :: irvec_q(3, nrr_q)
    !! coordinates of phononic WS points
    INTEGER, INTENT (in) :: ndegen_q(nrr_q, nat, nat)
    !! degeneracy of WS points
    !
    REAL(kind=DP), INTENT (in) :: xxq(3)
    !! kpoint coordinates for the interpolation
    REAL(kind=DP), INTENT (out) :: eig(nmodes)
    !! interpolated dynamical matrix eigenvalues for this kpoint
    COMPLEX(kind=DP), INTENT (out) :: cuf(nmodes, nmodes)
    !! Rotation matrix, fine mesh 
    !
    ! work variables
    !
    ! variables for lapack ZHPEVX
    INTEGER :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
    REAL(kind=DP) :: w( nmodes ), rwork( 7*nmodes )
    COMPLEX(kind=DP) :: champ( nmodes*(nmodes+1)/2 )
    COMPLEX(kind=DP) :: cwork( 2*nmodes ), cz( nmodes, nmodes)
    !
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: jmode
    !! Counter on modes
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    !
    REAL(kind=DP) :: xq(3)
    !! Coordinates q-point
    REAL(kind=DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}   
    REAL(kind=DP) :: massfac
    !! inverse square root of masses
    !
    COMPLEX(kind=DP) :: chf(nmodes, nmodes)
    ! Dynamical matrix in Bloch basis, fine mesh
    COMPLEX(kind=DP) :: cfac
    !! Complex prefactor for Fourier transform. 
    COMPLEX(KIND=DP) :: zdotu
    !! Dot product between the two phonon eigenvectors. 
    !
    CALL start_clock ( 'DynW2B' )
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  [Eqn. 32 of PRB 76, 165108 (2007)]
    !  D~(k')    = 1/ndegen(R) sum_R e^{ik'R     } D(R)
    !  D~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} D(R)
    !
    !  D~(k')    is chf ( nmodes, nmodes, 2*ik-1 )
    !  D~(k'+q') is chf ( nmodes, nmodes, 2*ik   )
    !
    xq = xxq
    chf(:,:) = czero
    !
    DO ir = 1, nrr_q
      !
      rdotk = twopi * dot_product( xq, dble(irvec_q( :, ir) ))
      ! 
      DO na = 1, nat
        DO nb = 1, nat
          IF ( ndegen_q(ir, na, nb) > 0 ) THEN
            cfac = exp( ci*rdotk ) / dble( ndegen_q(ir,na,nb) )
            ! To map atom coordinate to mode basis. 
            chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
              chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) &
            + cfac * rdw(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb, ir )
          ENDIF
        ENDDO
      ENDDO
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
       champ(imode + (jmode - 1) * jmode/2 ) = &
       ( chf( imode, jmode) + conjg ( chf( jmode, imode) ) ) * 0.5d0
     ENDDO
    ENDDO
    !
    !CALL zhpevx ('V', 'A', 'U', nmodes, champ, zero, zero, &
    !             0, 0, -one, neig, w, cz, nmodes, cwork, &
    !             rwork, iwork, ifail, info)
    CALL cdiagh2(nmodes,chf,nmodes,w,cz)
    ! 
    ! clean noise
    DO jmode=1,nmodes
      DO imode=1,nmodes
        IF ( ABS( cz(imode,jmode) ) < eps12 ) cz(imode,jmode) = czero
      ENDDO
    ENDDO
    ! 
    ! DS - Impose phase 
    IF (lphase) THEN
      DO jmode = 1,nmodes
        INNER : DO imode = 1,nmodes
          IF ( ABS(cz(imode, jmode)) > eps12 ) THEN
            cz(:, jmode) = cz(:, jmode) * conjg( cz(imode,jmode) )
            cz(:, jmode) = cz(:, jmode)/sqrt( zdotu(nmodes,conjg(cz(:, jmode)),1,cz(:, jmode),1) )
            EXIT INNER
          ENDIF
        END DO INNER
      ENDDO
    ENDIF
    !
    ! cuf(nmodes,nmodes) is rotation matrix (eigenmodes e_k)
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
    USE constants_epw, ONLY : twopi, czero, zero, one
    USE io_global, ONLY : stdout
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nmodes
    !! number of modes 
    INTEGER, INTENT (in) :: nrws
    !! Number of Wigner-Size real space vectors
    ! DBSP
    !INTEGER, INTENT (in) :: rws(0:3,nrws)
    REAL(kind=DP), INTENT (in) :: rws(0:3,nrws)
    !! Real space Wigner-Seitz vector
    REAL(kind=DP), INTENT (in) :: xxq (3)
    !! qpoint coordinates for the interpolation
    REAL(kind=DP), INTENT (out) :: eig (nmodes)
    !! interpolated phonon eigenvalues for this qpoint
    COMPLEX(kind=DP), INTENT (out) :: cuf(nmodes, nmodes)
    !! Rotation matrix, fine mesh 
    !
    ! work variables
    !
    ! variables for lapack ZHPEVX
    INTEGER :: neig, info, ifail( nmodes ), iwork( 5*nmodes )
    REAL(kind=DP) :: w( nmodes ), rwork( 7*nmodes )
    COMPLEX(kind=DP) :: champ( nmodes*(nmodes+1)/2 )
    COMPLEX(kind=DP) :: cwork( 2*nmodes ), cz( nmodes, nmodes)
    !
    LOGICAL, SAVE :: first=.true.
    !
    INTEGER :: n1,n2,n3, m1,m2,m3, i
    !! 
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: jmode
    !! Counter on modes
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: jpol
    !! Counter on polarizations
    !
    REAL(kind=DP) :: xq(3)
    !! Coordinates q-point
    REAL(kind=DP) :: massfac
    !! inverse square root of masses
    !
    REAL(kind=DP), EXTERNAL :: wsweight
    REAL(kind=DP), SAVE, ALLOCATABLE :: wscache(:,:,:,:,:)
    REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)
    !
    COMPLEX(kind=DP) :: chf(nmodes, nmodes)
    !! Dynamical matrix in Bloch basis, fine mesh
    COMPLEX(kind=DP) :: dyn(3,3,nat,nat)
    !! Dynamical matrix
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
            total_weight = zero
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
    !  [Eqn. 32 of PRB 76, 165108 (2007)]
    !  D~(k')    = 1/ndegen(R) sum_R e^{ik'R     } D(R)
    !  D~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} D(R)
    !
    !  D~(k')    is chf ( nmodes, nmodes, 2*ik-1 )
    !  D~(k'+q') is chf ( nmodes, nmodes, 2*ik   )
    !
    chf = czero
    dyn = czero
    !
    DO na=1, nat
       DO nb=1, nat
          total_weight=zero
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
        ( chf ( imode, jmode) + conjg ( chf ( jmode, imode) ) ) * 0.5d0
     ENDDO
    ENDDO
    !
    CALL zhpevx ('V', 'A', 'U', nmodes, champ , zero, zero, &
                 0, 0, -one, neig, w, cz, nmodes, cwork, &
                 rwork, iwork, ifail, info)
    !
    ! cuf(nmodes,nmodes) is rotation matrix (eigenmodes e_k)
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
    USE constants_epw, ONLY : twopi, czero, zero
    USE io_global, ONLY : stdout
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nmodes
    !! number of modes 
    INTEGER, INTENT (in) :: nrws
    !! Number of Wigner-Size real space vectors
    ! DBSP
    !!INTEGER, INTENT (in) :: rws(0:3,nrws)
    REAL(kind=DP), INTENT (in) :: rws(0:3,nrws)
    !! Wigner-Seitz radius 
    REAL(kind=DP), INTENT (in) :: xq(3)
    !! qpoint coordinates for the interpolation
    COMPLEX(kind=DP), INTENT (out) :: chf(nmodes, nmodes)
    !! dyn mat (not divided by the masses)
    !
    ! work variables
    !
    ! Dyn mat in Bloch basis, fine mesh
    LOGICAL,SAVE :: first=.true.
    !
    INTEGER :: n1,n2,n3, m1,m2,m3, i
    !! 
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: jpol
    !! Counter on polarizations
    !
    REAL(kind=DP), EXTERNAL :: wsweight
    REAL(kind=DP), SAVE, ALLOCATABLE :: wscache(:,:,:,:,:)
    REAL(kind=DP) total_weight, weight, arg, r(3), r_ws(3)  
    !
    COMPLEX(kind=DP) :: dyn(3,3,nat,nat)
    !! Dynamical matrix
    !
    ! bring xq in cart. coordinates
    !CALL cryst_to_cart (1, xq, bg, 1)
    !
    FIRST_TIME : IF (first) THEN
      first=.false.
      ALLOCATE( wscache(-2*nq3:2*nq3, -2*nq2:2*nq2, -2*nq1:2*nq1, nat,nat) )
      DO na=1, nat
         DO nb=1, nat
            total_weight = zero
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
          total_weight = zero
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
    USE constants_epw, ONLY : cone, czero, eps4
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nbnd 
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT (in) :: nrr 
    !! kpoint number for the interpolation
    !
    REAL(kind=DP), INTENT (in) :: etf(nbnd)
    !! Eigenenergies on the fine grid
    REAL(kind=DP), INTENT (in) :: etf_ks(nbnd) 
    !! Kohn-Sham eigenvalues
    !
    COMPLEX(kind=DP), INTENT (in) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh 
    COMPLEX(kind=DP), INTENT (out) :: dmef(3, nbnd, nbnd)
    !! interpolated dipole matrix elements in Bloch basis, fine mesh
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    !! Exponential factor
    ! 
    ! local variables
    !
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    !
    COMPLEX( kind=DP ) :: cdmef(3, nbnd, nbnd)
    !! dipole matrix elements in Bloch basis, fine mesh
    COMPLEX( kind=DP ) :: cdmef_tmp(nbnd, nbnd)
    !! dipole matrix elements in Bloch basis, fine mesh
    !
    ! Initialization
    cdmef_tmp(:,:) = czero
    dmef(:,:,:) = czero
    cdmef(:,:,:) = czero
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    !  p~(k')    = 1/ndegen(R) sum_R e^{ik'R     } p(R)
    !  p~(k'+q') = 1/ndegen(R) sum_R e^{i(k'+q')R} p(R)
    !  Note p~(k') is p^(W)(k') in PRB 74, 195118 (2006) notations
    !
    !  p~(k' )   is chf( nbnd, nbnd, 2*ik-1 )
    !  p~(k'+q') is chf( nbnd, nbnd, 2*ik   )
    !
    ! SUM on ir of: cdmef(1,ibnd,jbnd) = cdmef(1,ibnd,jbnd) + cfac(ir) * cdmew(1,ibnd,jbnd,ir)
    CALL zgemv('n', 3*(nbnd**2), nrr, cone, cdmew(:,:,:,:), 3*(nbnd**2), cfac(:), 1, cone, cdmef(:,:,:), 1  )
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! p(k') = U(k')^\dagger * p~(k') * U(k')
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    ! Note p(k') is p^(H)(k') in PRB 74, 195118 (2006) notations
    !
    DO ipol = 1, 3
      CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, cdmef(ipol,:,:), &
                 nbnd, cuf(:,:), nbnd, czero, cdmef_tmp(:,:), nbnd)
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:,:), &
                 nbnd, cdmef_tmp(:,:), nbnd, czero, dmef(ipol,:,:), nbnd)
    ENDDO
    !
    ! Satisfy
    ! Phys. Rev. B 62, 4927-4944 (2000) , Eq. (30)
    !
    IF (eig_read) THEN
       DO ibnd = 1, nbnd
         DO jbnd = 1, nbnd
           IF (abs(etf_ks(ibnd) - etf_ks(jbnd)) .gt. eps4) THEN
              dmef(:,ibnd,jbnd) = dmef(:,ibnd,jbnd) * &
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
    SUBROUTINE vmewan2bloch ( nbnd, nrr, irvec, cuf, vmef, etf, etf_ks, chw, cfac)
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
    !  RM 04/2018: optimized
    !
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    use elph2,         ONLY : cvmew 
    use cell_base,     ONLY : at, alat
    USE epwcom,        ONLY : eig_read
    USE constants_epw, ONLY : twopi, ci, czero, cone, zero, eps4, bohr2ang
    !USE io_global, ONLY : ionode_id
    !USE mp_world,  ONLY : mpime
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands (possibly of the optimal subspace)
    INTEGER, INTENT (in) :: nrr
    !! number of WS points
    INTEGER :: irvec(3, nrr)
    !! coordinates of WS points
    !
    REAL(kind=DP), INTENT (in) :: etf(nbnd)
    !! Eigenenergies on the fine grid
    REAL(kind=DP), INTENT (in) :: etf_ks(nbnd)
    !! Kohn-Sham eigenvalues
    !
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(kind=DP), INTENT (in) :: chw(nbnd, nbnd, nrr)
    !! Hamiltonian in Wannier basis
    COMPLEX(kind=DP), INTENT (in) :: cuf(nbnd, nbnd)
    !! Rotation matrix U^\dagger, fine mesh
    !
    COMPLEX(kind=DP), INTENT (out) :: vmef(3,nbnd,nbnd)
    !! interpolated velocity matrix elements in Bloch basis, fine mesh
    !
    ! local variables
    !
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ipol
    !! Counter on polarization
    !
    REAL(kind=DP) :: irvec_tmp(3)
    !! coordinates of WS points for the interpolation, cartesian coordinates
    !
    COMPLEX(kind=DP) :: chf_a(3, nbnd, nbnd)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(kind=DP) :: chf_a_tmp(nbnd, nbnd)
    !! derivative of interpolated hamiltonian eigenvalues, fine mesh
    COMPLEX(kind=DP) :: cvmef(3, nbnd, nbnd)
    !! velocity matrix elements in Bloch basis, fine mesh
    COMPLEX(kind=DP) :: cvmef_tmp(nbnd, nbnd)
    !! velocity matrix elements in Bloch basis, fine mesh
    !
    ! Initialization
    !
    cvmef_tmp(:,:) = czero
    cvmef(:,:,:) = czero
    vmef(:,:,:) = czero
    chf_a_tmp(:,:) = czero
    chf_a(:,:,:) = czero
    irvec_tmp(:) = zero
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform to fine k and k+q meshes
    !----------------------------------------------------------
    !
    ! [Eqn. 39 of PRB 74, 195118 (2006)]
    ! A^(W)_{mn,\alpha}(k') = 1/ndegen(R) sum_R e^{ik'R} r_{\alpha}(R)
    !
    ! SUM on ir of: cvmef(1,ibnd,jbnd) = cvmef(1,ibnd,jbnd) + cfac(ir) * cvmew(1,ibnd,jbnd,ir)
    CALL zgemv('n', 3*(nbnd**2), nrr, cone, cvmew(:,:,:,:), 3*(nbnd**2), cfac(:), 1, cone, cvmef(:,:,:), 1  )
    !
    ! k-derivative of the Hamiltonian in the Wannier gauge
    ! [Eqn. 38 of PRB 74, 195118 (2006)] or [Eq. 29 of PRB 75, 195121 (2007)]
    !
    ! dH~_{\alpha}(k')    = 1/ndegen(R) sum_R i*R_{\alpha} e^{ik'R     } H(R)
    ! dH~_{\alpha}(k'+q') = 1/ndegen(R) sum_R i*R_{\alpha} e^{i(k'+q')R} H(R)
    !
    ! Note dH~_{\alpha}(k') is H^(W)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations
    !
    ! dH~(k')    is chf_a( nbnd, nbnd, 2*ik-1 )
    ! dH~(k'+q') is chf_a( nbnd, nbnd, 2*ik   )
    !
    DO ir = 1, nrr
      !
      ! convert irvec from reduce to cartesian coordinates
      ! multiply by alat since the crystal axis 'at' are in 
      ! cart. coords. in units of a_0
      irvec_tmp(:) = alat * matmul( at, dble(irvec(:,ir)) )
      DO ipol = 1, 3
        chf_a(ipol,:,:) = chf_a(ipol,:,:) + &
              ci * irvec_tmp(ipol) * cfac(ir) * chw(:,:,ir)
      ENDDO
      !
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! A^(H)_{mn,\alpha}(k') = U(k')^\dagger A^(W)_{mn,\alpha}(k') U(k')
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    DO ipol = 1, 3
      !
      ! cvmef_tmp(:,:) = matmul( cvmef(ipol,:,:), conjg(transpose(cuf(:,:))) )
      ! vmef(ipol,:,:) = matmul( cuf(:,:), cvmef_tmp(:,:) )
      !
      CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, cvmef(ipol,:,:), &
                 nbnd, cuf(:,:), nbnd, czero, cvmef_tmp(:,:), nbnd)
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:,:), &
                 nbnd, cvmef_tmp(:,:), nbnd, czero, vmef(ipol,:,:), nbnd)
    ENDDO
    !
    ! [Eqn. 21 of PRB 74, 195118 (2006)]
    ! dH_{\alpha}(k') = U(k')^\dagger dH~_{\alpha}(k') U(k')
    !
    ! Note dH_{\alpha}(k') is H^(H)_{mn,\alpha}(k') in PRB 74, 195118 (2006) notations
    !
    ! U^\dagger is cuf(nbnd,nbnd), passed from hamwan2bloch.
    !
    DO ipol = 1, 3
      !
      ! chf_a_tmp(:,:) = matmul( chf_a(ipol,:,:), conjg(transpose(cuf(:,:))) )
      ! chf_a(ipol,:,:) = matmul(cuf(:,:), chf_a_tmp(:,:) )
      !
      CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, chf_a(ipol,:,:), &
                 nbnd, cuf(:,:), nbnd, czero, chf_a_tmp(:,:), nbnd)
      CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cuf(:,:), &
                 nbnd, chf_a_tmp(:,:), nbnd, czero, chf_a(ipol,:,:), nbnd)
    ENDDO
    !DBSP
    !if (mpime==1) write(902,*)SUM(irvec_tmp), celldm(1), at, SUM(irvec(:,ir))
    !if (mpime==1) write(902,*) 'celldm',celldm
    !if (mpime==1) write(902,*) 'celldm',alat
    !if (mpime==1) write(902,*) 'at',at
    !if (mpime==1) write(902,*)SUM(chw)
    !if (mpime==1) write(902,*)SUM(cuf(:,:))
    !if (mpime==1) write(902,*)SUM(vmef(:,:,:))
    !if (mpime==1) write(902,*)SUM(chf_a_tmp(:,:))
    !if (mpime==1) write(902,*)SUM(chf_a(:,:,:))

    !
    ! velocity matrix elements
    ! [Eqn. 31 of PRB 74, 195118 (2006)]
    ! \hbar v_{mn,\alpha}(k') = H^(H)_{mn,\alpha}(k') &
    !                         - (E^(H)_nk'-E^(H)_mk') * A^(H)_{mn,\alpha}(k')
    !
    ! RM - use etf instead of etf_ks when eig_read=.false.
    IF (eig_read) THEN
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          vmef(:,ibnd,jbnd) = chf_a(:,ibnd,jbnd) - &
               ci * ( etf_ks(jbnd) - etf_ks(ibnd) ) * vmef(:,ibnd,jbnd)
        ENDDO
      ENDDO
    ELSE
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          vmef(:,ibnd,jbnd) = chf_a(:,ibnd,jbnd) - &
              ci * ( etf(jbnd) - etf(ibnd) ) * vmef(:,ibnd,jbnd)
        ENDDO
      ENDDO
    ENDIF
    !
    ! Satisfy
    ! Phys. Rev. B 62, 4927-4944 (2000) , Eq. (30)
    !
    IF (eig_read) THEN
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          IF (abs(etf_ks(ibnd) - etf_ks(jbnd)) .gt. eps4) THEN
            vmef(:,ibnd,jbnd) = vmef(:,ibnd,jbnd) * &
                   ( etf(ibnd)    - etf(jbnd) )/ &
                   ( etf_ks(ibnd) - etf_ks(jbnd) )
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    END SUBROUTINE vmewan2bloch
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp ( nmodes, xxq, irvec_g, ndegen_g, nrr_g, cuf, epmatf, nbnd, nrr_k )
    !---------------------------------------------------------------------------
    !!
    !! even though this is for phonons, I use the same notations
    !! adopted for the electronic case (nmodes->nmodes etc)
    !!
    USE kinds,            ONLY : DP
    USE epwcom,           ONLY : etf_mem
    USE elph2,            ONLY : epmatwp
    USE constants_epw,    ONLY : twopi, ci, czero, cone
    USE io_epw,           ONLY : iunepmatwp, iunepmatwp2
    USE mp_global,        ONLY : mp_sum
    USE mp_world,         ONLY : world_comm
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
#endif
    USE ions_base,        ONLY : nat
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT (in) :: nrr_g
    !! Number of phononic WS points
    INTEGER, INTENT (in) :: irvec_g ( 3, nrr_g)
    !! Coordinates of WS points
    INTEGER, INTENT (in) :: ndegen_g (nrr_g, nat)
    !! Number of degeneracy of WS points
    INTEGER, INTENT (in) :: nbnd
    !! Number of bands
    INTEGER, INTENT (in) ::  nrr_k
    !! Number of electronic WS points
    REAL(kind=DP) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(kind=DP), INTENT (in) :: cuf(nmodes, nmodes)
    !! e-p matrix in Wanner representation
    COMPLEX(kind=DP), INTENT (out) :: epmatf(nbnd, nbnd, nrr_k, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    ! 
    ! Local variables 
    !
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: irn
    !! Combined WS and atom index  
    INTEGER :: ir_start
    !! Starting ir for this cores
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: ierr
    !! Return if there is an error
    INTEGER :: na
    !! Atom index
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
    COMPLEX(kind=DP) :: cfac(nat,nrr_g)
    !! Factor for the FT
    COMPLEX(kind=DP), ALLOCATABLE :: epmatw( :,:,:,:)
    !! El-ph matrix elements
    !
    CALL start_clock('ephW2Bp')
    ! 
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(R_e,q') = 1/ndegen(R_p) sum_R_p e^{iq'R_p} g(R_e,R_p)
    !
    !  g~(R_e,q') is epmatf(nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    ! SP: Because nrr_g can be quite small, we do a combined parallelization
    !     on WS vector and atoms  
    CALL para_bounds(ir_start, ir_stop, nrr_g * nat)
    !
    eptmp(:,:,:,:) = czero
    cfac(:,:) = czero
    !
    DO irn = ir_start, ir_stop
      ir = (irn-1)/nat + 1
      na = MOD(irn-1,nat) +1   
      !   
      ! note xxq is assumed to be already in cryst coord
      !
      rdotk = twopi * dot_product ( xxq, dble(irvec_g(:, ir)) )
      IF (ndegen_g(ir,na) > 0) &
        cfac(na,ir) = exp( ci*rdotk ) / dble( ndegen_g(ir,na) )
    ENDDO
    ! 
    IF (etf_mem == 0) then
      !      
      DO irn = ir_start, ir_stop
        ir = (irn-1)/nat + 1
        na = MOD(irn-1,nat) + 1
        ! 
        !print*,'irn ',irn, shape(cfac), shape(epmatwp(:,:,:,:,:)), 3*(na-1)+1,3*na, ir
        CALL ZAXPY(nbnd * nbnd * nrr_k * 3, cfac(na,ir), epmatwp(:,:,:,3*(na-1)+1:3*na,ir), 1, &
             eptmp(:,:,:,3*(na-1)+1:3*na), 1)
        !CALL zgemv( 'n',  nbnd * nbnd * nrr_k * 3, ir_stop - ir_start + 1, cone, &
        !     epmatwp(:,:,:,3*(na-1)+1:3*na,ir_start:ir_stop), nbnd * nbnd * nrr_k * 3, &
        !     cfac(ir_start:ir_stop), 1, czero, eptmp(:,:,:,3*(na-1)+1:3*na), 1 )
      ENDDO
      !
    ELSE
      !
      !ALLOCATE(epmatw ( nbnd, nbnd, nrr_k, 3))
      !
      !lrepmatw2   = 2 * nbnd * nbnd * nrr_k * nmodes
#if defined(__MPI)
      ALLOCATE(epmatw ( nbnd, nbnd, nrr_k, 3))
      ! Although this should almost never be problematic (see explaination below)
      lrepmatw2 = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                      INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                      INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
                                      3_MPI_OFFSET_KIND
                                      !INT( nmodes, kind = MPI_OFFSET_KIND )
#else
      ALLOCATE(epmatw ( nbnd, nbnd, nrr_k, nmodes))
      lrepmatw2 = INT( 2 * nbnd * nbnd * nrr_k * 3, kind = 8)
#endif
      ! 
      DO irn = ir_start, ir_stop
        ir = (irn-1)/nat + 1
        na = MOD(irn-1,nat) + 1      
#if defined(__MPI)
        ! DEBUG: print*,'Process ',my_id,' do ',ir,'/ ',ir_stop
        !
        !  Direct read of epmatwp for this ir
        !lrepmatw   = 2 * nbnd * nbnd * nrr_k * nmodes * 8 * (ir-1)
        ! 
        ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or 
        !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
        !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
        lrepmatw = 2_MPI_OFFSET_KIND * 8_MPI_OFFSET_KIND * &
                                       INT( nbnd  ,kind=MPI_OFFSET_KIND ) * &
                                       INT( nbnd  ,kind=MPI_OFFSET_KIND ) * &
                                       INT( nrr_k ,kind=MPI_OFFSET_KIND ) * &
                                       ( INT( 3_MPI_OFFSET_KIND * ( na - 1_MPI_OFFSET_KIND) ,kind=MPI_OFFSET_KIND ) + &
        INT( 3_MPI_OFFSET_KIND * nat ,kind=MPI_OFFSET_KIND ) * ( INT( ir,kind=MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND ) )


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
        !print*,'irn ir na ',irn, ir, na
        !print*,'lrepmatw ',lrepmatw
        !print*,'shape epmatw ', shape(epmatw), sum(epmatw)
        CALL ZAXPY(nbnd * nbnd * nrr_k * 3, cfac(na,ir), epmatw(:,:,:,:), 1, &
                eptmp(:,:,:,3*(na-1)+1:3*na), 1)        
        ! 
#else      
        CALL rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
        !
        CALL ZAXPY(nbnd * nbnd * nrr_k * 3, cfac(na,ir), epmatw(:,:,:,3*(na-1)+1:3*na), 1, &
                eptmp(:,:,:,3*(na-1)+1:3*na), 1)         

#endif
        ! 
      ENDDO
      DEALLOCATE(epmatw)
    ENDIF
    !
#if defined(__MPI)
    CALL mp_sum(eptmp, world_comm)
#endif  
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    ! [Eqn. 22 of PRB 76, 165108 (2007)]
    ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
    !
    Call zgemm( 'n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, cone, eptmp, & 
                nbnd * nbnd * nrr_k, cuf, nmodes, czero, epmatf, nbnd * nbnd * nrr_k )

    !
    CALL stop_clock('ephW2Bp')
    !
    END SUBROUTINE ephwan2blochp
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch ( nbnd, nrr, epmatw, cufkk, cufkq, &
           epmatf, nmodes, cfac )
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
    INTEGER, INTENT (in) :: nmodes
    !! number of phonon modes
    !
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(kind=DP), INTENT (in) :: epmatw( nbnd, nbnd, nrr, nmodes)
    !! e-p matrix in Wannier representation
    COMPLEX(kind=DP), INTENT (in) :: cufkk(nbnd, nbnd)
    !! rotation matrix U(k)^\dagger, fine k mesh
    COMPLEX(kind=DP), INTENT (in) :: cufkq(nbnd, nbnd)
    !! rotation matrix U(k+q)^\dagger, fine q mesh
    COMPLEX(kind=DP), INTENT (out) :: epmatf(nbnd, nbnd, nmodes)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! work variables 
    !
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: imode
    !! Counter on  phonon modes
    !
    COMPLEX(kind=DP) :: eptmp( nbnd, nbnd)
    !! Temporary variable
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(k',q') = 1/ndegen(R_e) sum_R_e e^{ik'R_e} g(R_e,q')
    !
    !  g~(k',q') is epmatf(nmodes, nmodes, ik)
    !  every pool works with its own subset of k points on the fine grid
    !
    epmatf = czero
    !
    !DO ir = 1, nrr
    !   !
    !   ! rdotk = twopi * dot_product ( xk, dble(irvec(:,ir)) )
    !   ! cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
    !   !
    !   DO imode = 1, nmodes
    !     epmatf(:,:,imode) = epmatf(:,:,imode) + cfac(ir) * epmatw(:,:,ir,imode)
    !   ENDDO
    !   !
    !ENDDO
    !
    DO imode = 1, nmodes
      ! SUM on ir of: epmatf(ibnd,jbnd,imode) = epmatf(ibnd,jbnd,imode) + cfac(ir) * epmatw(ibnd,jbnd,ir,imode)
      CALL zgemv('n', nbnd**2, nrr, cone, epmatw(:,:,:,imode), nbnd**2, cfac(:), 1, cone, epmatf(:,:,imode), 1 )
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g(k',q') = U(k'+q') * g~(k',q') * U(k')^\dagger
    !
    !  RM - this is what is calculated
    !  g(k',q') = U(k'+q')^\dagger * g~(k',q') * U(k')
    !
    !  the two zgemm calls perform the following ops:
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
    END SUBROUTINE ephwan2bloch
    ! 
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2bloch_mem ( nbnd, nrr, epmatw, cufkk, cufkq, &
           epmatf, cfac )
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
    !
    COMPLEX(kind=DP), INTENT (in) :: cfac(nrr)
    !! Exponential factor
    COMPLEX(kind=DP), INTENT (in) :: epmatw( nbnd, nbnd, nrr)
    !! e-p matrix in Wannier representation
    COMPLEX(kind=DP), INTENT (in) :: cufkk(nbnd, nbnd)
    !! rotation matrix U(k)^\dagger, fine k mesh
    COMPLEX(kind=DP), INTENT (in) :: cufkq(nbnd, nbnd)
    !! rotation matrix U(k+q)^\dagger, fine k mesh
    COMPLEX(kind=DP), INTENT (out) :: epmatf(nbnd, nbnd)
    !! e-p matrix in Bloch representation, fine grid
    !
    ! work variables 
    !
    INTEGER :: ir
    !! Counter on real-space index
    !
    COMPLEX(kind=DP) :: eptmp( nbnd, nbnd)
    !! Temporary variable
    !
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(k',q') = 1/ndegen(R_e) sum_R_e e^{ik'R_e} g(R_e,q')
    !
    !  g~(k',q') is epmatf(nmodes, nmodes, ik)
    !  every pool works with its own subset of k points on the fine grid
    !
    epmatf = czero
    !
    !DO ir = 1, nrr
    !   !
    !   ! note xk is assumed to be already in cryst coord
    !   !
    !   ! rdotk = twopi * dot_product ( xk, dble(irvec(:,ir)) )
    !   ! cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
    !   !
    !   epmatf(:,:) = epmatf(:,:) + cfac(ir) * epmatw(:,:,ir)
    !   !
    !ENDDO
    !
    ! SUM on ir of: epmatf(ibnd,jbnd) = epmatf(ibnd,jbnd) + cfac(ir) * epmatw(ibnd,jbnd,ir)
    CALL zgemv('n', nbnd**2, nrr, cone, epmatw(:,:,:), nbnd**2, cfac(:), 1, cone, epmatf(:,:), 1 )
    !
    !----------------------------------------------------------
    !  STEP 4: un-rotate to Bloch space, fine grid
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g(k',q') = U(k'+q') * g~(k',q') * U(k')^\dagger
    !
    !  RM - this is what is calculated
    !  g(k',q') = U(k'+q')^\dagger * g~(k',q') * U(k')
    !
    !  the two zgemm calls perform the following ops:
    !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
    !
    !
    CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
               nbnd, epmatf (:,:), nbnd, czero, eptmp, nbnd)
    CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, eptmp, &
               nbnd, cufkk, nbnd, czero, epmatf(:,:), nbnd)
    !
    END SUBROUTINE ephwan2bloch_mem
    ! 
    !---------------------------------------------------------------------------
    SUBROUTINE ephwan2blochp_mem (imode, nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatf, nbnd, nrr_k )
    !---------------------------------------------------------------------------
    !!
    !! Even though this is for phonons, I use the same notations
    !! adopted for the electronic case (nmodes->nmodes etc)
    !!
    USE kinds,            ONLY : DP
    USE constants_epw,    ONLY : twopi, ci, czero
    USE io_files,         ONLY : prefix, tmp_dir
    USE mp_global,        ONLY : mp_sum
    USE mp_world,         ONLY : world_comm
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, &
                                 MPI_MODE_RDONLY,MPI_INFO_NULL
#endif
    USE ions_base,        ONLY : nat

    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: imode
    !! Current mode  
    INTEGER, INTENT (in) :: nmodes
    !! Total number of modes
    INTEGER, INTENT (in) :: nrr_g
    !! Number of phononic WS points
    INTEGER, INTENT (in) :: irvec_g( 3, nrr_g)
    !! Coordinates of WS points
    INTEGER, INTENT (in) :: ndegen_g(nrr_g, nat)
    !! Number of degeneracy of WS points
    INTEGER, INTENT (in) :: nbnd
    !! Number of bands
    INTEGER, INTENT (in) ::  nrr_k
    !! Number of electronic WS points
    REAL(kind=DP) :: xxq(3)
    !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
    COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nrr_k)
    !! e-p matrix in Bloch representation, fine grid
    ! 
    ! Local variables 
    !
    CHARACTER (len=256) :: filint
    !! File name
    !
    INTEGER :: ir
    !! Real space WS index
    INTEGER :: ir_start
    !! Starting ir for this pool
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: iunepmatwp2
    !! Return the file unit
    INTEGER :: ierr
    !! Return if there is an error
    INTEGER :: na
    !! Index on atom
#if defined(__MPI)  
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw2
    !! Offset to tell where to start reading the file
#endif  
    !
    REAL(kind=DP) :: rdotk
    !! Exponential for the FT
    !
    COMPLEX(kind=DP) :: cfac(nrr_g)
    !! Factor for the FT
    COMPLEX(kind=DP), ALLOCATABLE :: epmatw( :,:,:)
    !! El-ph matrix elements
    !
    CALL start_clock('ephW2Bp')
    !----------------------------------------------------------
    !  STEP 3: inverse Fourier transform of g to fine k mesh
    !----------------------------------------------------------
    !
    !  [Eqn. 22 of PRB 76, 165108 (2007)]
    !  g~(R_e,q') = 1/ndegen(R_p) sum_R_p e^{iq'R_p} g(R_e,R_p)
    !
    !  g~(R_e,q') is epmatf(nmodes, nmodes, ik )
    !  every pool works with its own subset of k points on the fine grid
    !
    CALL para_bounds(ir_start, ir_stop, nrr_g)
    !
#if defined(__MPI)  
    filint = trim(tmp_dir)//trim(prefix)//'.epmatwp1'
    CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp_mem', 'error in MPI_FILE_OPEN',1 )
#endif  
    !
    cfac(:) = czero
    !
    DO ir = ir_start, ir_stop
      !   
      ! note xxq is assumed to be already in cryst coord
      rdotk = twopi * dot_product ( xxq, dble(irvec_g(:, ir)) )
      na = (imode - 1) / 3 + 1
      IF (ndegen_g(ir, na) > 0) &
        cfac(ir) = exp( ci*rdotk ) / dble( ndegen_g(ir, na) )
    ENDDO
    ! 
    ALLOCATE(epmatw( nbnd, nbnd, nrr_k))
    !
#if defined(__MPI)  
    lrepmatw2 = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                    INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                    INT( nrr_k , kind = MPI_OFFSET_KIND )
#endif                          
    ! 
    DO ir = ir_start, ir_stop
      !
      ! SP: The following needs a small explaination: although lrepmatw is correctly defined as kind 8 bits or 
      !     kind=MPI_OFFSET_KIND, the number "2" and "8" are default kind 4. The other as well. Therefore
      !     if the product is too large, this will crash. The solution (kind help recieved from Ian Bush) is below:
#if defined(__MPI)    
      lrepmatw = 2_MPI_OFFSET_KIND * INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                     INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                     INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
                                     INT( nmodes, kind = MPI_OFFSET_KIND ) * &
               8_MPI_OFFSET_KIND * ( INT( ir    , kind = MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND ) + &
               2_MPI_OFFSET_KIND *   INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                     INT( nbnd  , kind = MPI_OFFSET_KIND ) * &
                                     INT( nrr_k , kind = MPI_OFFSET_KIND ) * &
               8_MPI_OFFSET_KIND * ( INT( imode , kind = MPI_OFFSET_KIND ) - 1_MPI_OFFSET_KIND )
      !
      ! SP: mpi seek is used to set the position at which we should start
      ! reading the file. It is given in bits. 
      ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
      !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ. 
      !        Here we want non blocking because not all the process have the same nb of ir. 
      !
      CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
#endif    
      ! 
      !
      CALL ZAXPY(nbnd * nbnd * nrr_k, cfac(ir), epmatw, 1, epmatf, 1)
      ! 
    ENDDO
    DEALLOCATE(epmatw)
    !
    CALL mp_sum(epmatf, world_comm)
    !
#if defined(__MPI)  
    CALL MPI_FILE_CLOSE(iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp_mem', 'error in MPI_FILE_CLOSE',1 )
#endif  
    !
    CALL stop_clock('ephW2Bp')
    !
    END SUBROUTINE ephwan2blochp_mem
    ! 
  END MODULE wan2bloch
