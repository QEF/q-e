MODULE adduscore

 USE kinds,          ONLY : DP
! change for US
!    USE paw_variables,  ONLY : paw_info
 IMPLICIT NONE

 PUBLIC :: US_make_ae_charge
 PRIVATE

 !maxZ is the atomic number of Rn, the last
 ! atom that has been tabulated with effective STOs so far.
 INTEGER, PARAMETER :: maxZ=86
 !maxSTO is the cumulative maximum number of core STO for the
 ! first 84 elements after He (Li-Rn, H and He have no core).
 !This number is in the minimal (-n-) assumption of pseudization
 ! (maximal number of core orbitals), i.e.
 ! with (n-2)f,(n-1)d,n(s,p) in valence.
 ! e.g. Pt val=4f,5d,6s
 INTEGER, PARAMETER :: maxSTO=392
 REAL (DP), PARAMETER :: rc_all=3.d0 !core radius for all atoms
 INTEGER, EXTERNAL :: atomic_number
 REAL (DP) :: Zcmax(maxZ)
 REAL (DP) :: nstar(maxSTO)
 REAL (DP) :: normc(maxSTO)
 REAL (DP) :: occ(maxSTO)
 REAL (DP) :: zeta(maxSTO)
 INTEGER   :: nSTO(maxZ)
 INTEGER   :: iSTO(maxZ)

 CONTAINS

! This is paw_postproc.f90 modified by Giovanni La Penna
! (HPC-Europa, Cork (IE), 28/01/2019).
! Reconstruct the all-electron density 
! using the effective Slater approximation for atomic cores.
! It assumes a core electron density associated to the isolated atom
! (frozen core) and described with effective Slater orbitals.
! Effective STO (Slater type orbital) data:
!  see Eq. 9-92 of Eyring-Walter-Kimball, Quantum chemistry (1944)
! https://archive.org/details/quantumchemistry033517mbp
! (Zener, Phys.Rev., 36,51(1930);Slater Phys.Rev.36,57(1930) */
!
! In practice: this creates on-the-fly what is passed as 
!  upf(i%t)%paw%ae_rho_atc(ir) / nspin * (2._DP * sqrtpi)
! in the original paw_postproc.f90
! for atoms, with the information in the UPF file.

SUBROUTINE US_make_ae_charge(rho)
   !
   USE ions_base,         ONLY : nat, ntyp => nsp, ityp, tau
   USE lsda_mod,          ONLY : nspin
   USE uspp,              ONLY : okvan
   USE uspp_param,        ONLY : nh, nhm, upf
   USE scf,               ONLY : scf_type
   USE fft_base,          ONLY : dfftp
   USE fft_types,         ONLY : fft_index_to_3d
   USE mp_global,         ONLY : me_pool
   USE cell_base,         ONLY : at, bg, alat, omega
   USE io_global,         ONLY : stdout, ionode

   TYPE(scf_type), INTENT(inout) :: rho
   !
   INTEGER, PARAMETER      :: Nc=16
   INTEGER                 :: ipol             ! counter on x,y,z
   INTEGER                 :: ir               ! counter on grid point
   INTEGER                 :: ir1,ir2          ! counters on atomic core mesh
   INTEGER                 :: is               ! counter on spin component
   INTEGER                 :: j,k,l
   INTEGER                 :: i, ia
   LOGICAL                 :: offrange
   REAL(DP)                :: posi0(3),posi1(3),posi2(3),posi3(3),posi4(3),&
                              posi5(3),posi6(3),posi7(3)
   REAL(DP)                :: inv_nr1, inv_nr2, inv_nr3, &
                              distsqmin,dist(8),distsq(8),rhoat
   REAL(DP)                :: rcsq,dr,sum1,sum2,dvol,X0(3),&
                              rc(3),distc,distsqc,inv_Nc,alatc,&
                              volc,inv_nrc1,inv_nrc2,inv_nrc3
   INTEGER                 :: imin,ic,jc,kc

   ! Some initialization
   !
   inv_nr1 = 1.D0 / dble(  dfftp%nr1 )
   inv_nr2 = 1.D0 / dble(  dfftp%nr2 )
   inv_nr3 = 1.D0 / dble(  dfftp%nr3 )
   inv_Nc=1.d0/REAL(Nc)
   volc=inv_Nc*inv_Nc*inv_Nc
   inv_nrc1 = inv_Nc*inv_nr1
   inv_nrc2 = inv_Nc*inv_nr2
   inv_nrc3 = inv_Nc*inv_nr3
   dvol = omega*inv_nr1*inv_nr2*inv_nr3
   !
   ! I cannot parallelize on atoms, because it is already parallelized
   ! on charge slabs
   !

   ifUS: IF (okvan) THEN
      !
      ! sets all cores up (effective STOs, Z<=Z(Rn)).
      ! All atoms are assumed described with US pseudopotentials.
      ! Compute rho of core contribution at point r
      ! This computes the data that are, for PAW, in:
      !  upf(i%t)%paw%ae_rho_atc(ir) / nspin * (2._DP * sqrtpi)
      !
      CALL init_cores(ntyp,rcsq)
      !
   atoms: DO ia = 1, nat
         sum2=0.d0
         i=ityp(ia)    ! type of atom ia
         !
         ir2=0
         rsp_point : DO ir = 1, dfftp%nr1x * dfftp%my_nr2p * dfftp%my_nr3p
            !
            ! three dimensional indices (l,j,k)
            !
            CALL fft_index_to_3d (ir, dfftp, l,j,k, offrange)
            !
            ! ... do not include points outside the physical range!
            IF ( offrange ) CYCLE rsp_point
            !
            ! using the 8 vertices to find when the grid element
            !  is out of atomic core
            DO ipol = 1, 3
               posi0(ipol) = dble( l )*inv_nr1*at(ipol,1) + &
                             dble( j )*inv_nr2*at(ipol,2) + &
                             dble( k )*inv_nr3*at(ipol,3)
               posi1(ipol) = dble( l+1 )*inv_nr1*at(ipol,1) + &
                             dble( j )*inv_nr2*at(ipol,2) + &
                             dble( k )*inv_nr3*at(ipol,3)
               posi2(ipol) = dble( l+1 )*inv_nr1*at(ipol,1) + &
                             dble( j+1 )*inv_nr2*at(ipol,2) + &
                             dble( k )*inv_nr3*at(ipol,3)
               posi3(ipol) = dble( l )*inv_nr1*at(ipol,1) + &
                             dble( j+1 )*inv_nr2*at(ipol,2) + &
                             dble( k )*inv_nr3*at(ipol,3)
               posi4(ipol) = dble( l )*inv_nr1*at(ipol,1) + &
                             dble( j )*inv_nr2*at(ipol,2) + &
                             dble( k+1 )*inv_nr3*at(ipol,3)
               posi5(ipol) = dble( l+1 )*inv_nr1*at(ipol,1) + &
                             dble( j )*inv_nr2*at(ipol,2) + &
                             dble( k+1 )*inv_nr3*at(ipol,3)
               posi6(ipol) = dble( l+1 )*inv_nr1*at(ipol,1) + &
                             dble( j+1 )*inv_nr2*at(ipol,2) + &
                             dble( k+1 )*inv_nr3*at(ipol,3)
               posi7(ipol) = dble( l )*inv_nr1*at(ipol,1) + &
                             dble( j+1 )*inv_nr2*at(ipol,2) + &
                             dble( k+1 )*inv_nr3*at(ipol,3)
            ENDDO
            X0(:) = posi0(:)
            !
            ! find the distance of real-space grid's point ir w.r.t
            ! closer periodic image of atom ia
            !
!            write(stdout,'(a,f,f,f)') "posi0= ",posi0(1),posi0(2),posi0(3)
            posi0(:) = posi0(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi0, bg, -1 )
            posi0(:) = posi0(:) - anint( posi0(:) )
            CALL cryst_to_cart( 1, posi0, at, 1 )
            posi1(:) = posi1(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi1, bg, -1 )
            posi1(:) = posi1(:) - anint( posi1(:) )
            CALL cryst_to_cart( 1, posi1, at, 1 )
            posi2(:) = posi2(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi2, bg, -1 )
            posi2(:) = posi2(:) - anint( posi2(:) )
            CALL cryst_to_cart( 1, posi2, at, 1 )
            posi3(:) = posi3(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi3, bg, -1 )
            posi3(:) = posi3(:) - anint( posi3(:) )
            CALL cryst_to_cart( 1, posi3, at, 1 )
            posi4(:) = posi4(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi4, bg, -1 )
            posi4(:) = posi4(:) - anint( posi4(:) )
            CALL cryst_to_cart( 1, posi4, at, 1 )
            posi5(:) = posi5(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi5, bg, -1 )
            posi5(:) = posi5(:) - anint( posi5(:) )
            CALL cryst_to_cart( 1, posi5, at, 1 )
            posi6(:) = posi6(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi6, bg, -1 )
            posi6(:) = posi6(:) - anint( posi6(:) )
            CALL cryst_to_cart( 1, posi6, at, 1 )
            posi7(:) = posi7(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi7, bg, -1 )
            posi7(:) = posi7(:) - anint( posi7(:) )
            CALL cryst_to_cart( 1, posi7, at, 1 )
            !
            posi0(:) = posi0(:) * alat
            posi1(:) = posi1(:) * alat
            posi2(:) = posi2(:) * alat
            posi3(:) = posi3(:) * alat
            posi4(:) = posi4(:) * alat
            posi5(:) = posi5(:) * alat
            posi6(:) = posi6(:) * alat
            posi7(:) = posi7(:) * alat
            distsq(1) = posi0(1)**2 + posi0(2)**2 + posi0(3)**2
            distsq(2) = posi1(1)**2 + posi1(2)**2 + posi1(3)**2
            distsq(3) = posi2(1)**2 + posi2(2)**2 + posi2(3)**2
            distsq(4) = posi3(1)**2 + posi3(2)**2 + posi3(3)**2
            distsq(5) = posi4(1)**2 + posi4(2)**2 + posi4(3)**2
            distsq(6) = posi5(1)**2 + posi5(2)**2 + posi5(3)**2
            distsq(7) = posi6(1)**2 + posi6(2)**2 + posi6(3)**2
            distsq(8) = posi7(1)**2 + posi7(2)**2 + posi7(3)**2

            distsqmin=MINVAL(distsq)
!            imin=MINLOC(distsq,1)
            ! don't consider grid points too far from the atom:
            IF ( distsqmin > rcsq ) CYCLE rsp_point
            ! LEGO-like integral within each grid point
            sum1=0.d0
            DO ic=0,Nc-1
             DO jc=0,Nc-1
              DO kc=0,Nc-1
               DO ipol = 1, 3
                rc(ipol) = dble( ic )*inv_nrc1*at(ipol,1) + &
                           dble( jc )*inv_nrc2*at(ipol,2) + &
                           dble( kc )*inv_nrc3*at(ipol,3)
                rc(ipol)=X0(ipol)+rc(ipol)
               ENDDO
               rc(:) = rc(:) - tau(:,ia)
               CALL cryst_to_cart( 1, rc, bg, -1 )
               rc(:) = rc(:) - anint( rc(:) )
               CALL cryst_to_cart( 1, rc, at, 1 )
               rc(:) = rc(:) * alat
               distsqc = rc(1)**2 + rc(2)**2 + rc(3)**2
               if(distsqc.le.rcsq) THEN
                distc=sqrt(distsqc)
!                ir1=int(distc/dr)+1
!                rhoat=rho_atc(ir1,1,i)
                CALL do_core(i,distc,rhoat)
                sum1=sum1+rhoat
                ir2=ir2+1
!                write(stdout,'(a,f,f,f)') "RC= ",rc(1),rc(2),rc(3)
!                WRITE(stdout,'(a,i,f,e)') "rho_atc: ",ia,distc,rhoat
               ENDIF
              ENDDO
             ENDDO
            ENDDO
            sum1=sum1*volc
!            WRITE(stdout,'(a,i,f,e)') "rho_atc: ",ia,sum1
            sum2=sum2+sum1
            !
            ! generate the atomic charge on point posi(:), which means
            ! sum over effective STOs for each atom.
            ! Interpolate the radial function at distance |posi(:)|
            !
            DO is = 1,nspin
             rho%of_r(ir,is)= rho%of_r(ir,is) + sum1
            ENDDO
         ENDDO rsp_point
         !
!   WRITE(stdout,*) 'Number of points: ',ir2,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
!   WRITE(stdout,*) 'volume: ',ir2*dvol*volc,fpi*sqrt(rcsq)**3/3.0_dp
!   WRITE(stdout,'(a,i3,a,f10.5)') 'Core: ',ia,' Number of electrons: ',sum2*dvol
   ENDDO atoms
   !
 ENDIF ifUS

END SUBROUTINE US_make_ae_charge

SUBROUTINE init_cores(ntyp,rcsq)
 
! Initializes cores of each atomic type.
! Distance units are bohr
 
 USE kinds,          ONLY : DP
 USE uspp_param,     ONLY : upf
 USE io_global,      ONLY : stdout, ionode
 USE ions_base,      ONLY : atom_label => atm
 USE constants,      ONLY : pi

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ntyp
 REAL (DP), INTENT(OUT) :: rcsq

 LOGICAL   :: any_uspp
 INTEGER   :: i,ia,nlm,ir,isc,ic
 REAL (DP) :: Z,Zc,Zv,Zsc,r,sum1,sum2,sum3,val,val2,xp,norm_1s

 rcsq = rc_all*rc_all
 norm_1s = 1./sqrt(4.d0*pi)

 any_uspp = any(upf(1:ntyp)%tvanp)
 IF(.not.any_uspp) THEN
  CALL errore('init_cores',' USPP not found.')
 ENDIF
 nlm=1
 DO i=1,ntyp
  iSTO(i)=nlm-1
  ia=atomic_number(atom_label(i))
  Z=REAL(ia)
  Zv=upf(i)%zp
  Zc=Z-Zv
! 2019-10-11
  ic=int(Zc)
!
  IF(ia.le.2) THEN
   ! no core
   nSTO(i)=0
   Zcmax(i)=0
  ELSE IF (ia.gt.2 .and. ia.le.10) THEN
   nSTO(i)=1
   Zcmax(i)=2
   !1s
   nstar(nlm)=0.
   zeta(nlm)=Z-0.30
   normc(nlm)=norm_1s*2.*zeta(nlm)*sqrt(zeta(nlm))
   occ(nlm)=2.
   nlm=nlm+1
  ELSE IF(ia.gt.10 .and. ia.le.18) THEN
   IF(ic.gt.2) THEN
    nSTO(i)=2
    Zcmax(i)=10
    ! 2s,2p 
    nstar(nlm)=1.
    zeta(nlm)=Z-2.*0.85-7.*0.35
    normc(nlm)=4.*norm_1s*zeta(nlm)*zeta(nlm)*sqrt(zeta(nlm)/12.)
    occ(nlm)=8.
    nlm=nlm+1
    ! 1s
    nstar(nlm)=0.
    zeta(nlm)=Z-0.30
    normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
    occ(nlm)=2.
    nlm=nlm+1
   ELSE
    ! 1s
    nSTO(i)=1
    Zcmax(i)=2
    nstar(nlm)=0.
    zeta(nlm)=Z-0.30
    normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
    occ(nlm)=2.
    nlm=nlm+1
   ENDIF
  ELSE IF(ia.gt.18 .and. ia.le.36) THEN
   IF(ic.gt.2 .and. ic.le.10) THEN
    nSTO(i)=3
    Zcmax(i)=18
    ! 3s,3p
    nstar(nlm)=2.
    zeta(nlm)=Z-2.*1.00-8.*0.85-7*0.35
    normc(nlm)=8.*norm_1s*((zeta(nlm))**3)*sqrt(zeta(nlm)/360.)
    occ(nlm)=8.
    nlm=nlm+1
    ! 2s,2p
    nstar(nlm)=1.
    zeta(nlm)=Z-2.*0.85-7.*0.35
    normc(nlm)=4.*norm_1s*zeta(nlm)*zeta(nlm)*sqrt(zeta(nlm)/12.)
    occ(nlm)=8.
    nlm=nlm+1
    ! 1s
    nstar(nlm)=0.
    zeta(nlm)=Z-0.30
    normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
    occ(nlm)=2.
    nlm=nlm+1
   ELSE IF(ic.gt.10) THEN
    nSTO(i)=2
    Zcmax(i)=10
    ! 2s,2p
    nstar(nlm)=1.
    zeta(nlm)=Z-2.*0.85-7.*0.35
    normc(nlm)=4.*norm_1s*zeta(nlm)*zeta(nlm)*sqrt(zeta(nlm)/12.)
    occ(nlm)=8.
    nlm=nlm+1
    ! 1s
    nstar(nlm)=0.
    zeta(nlm)=Z-0.30
    normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
    occ(nlm)=2.
    nlm=nlm+1
   ELSE
    nSTO(i)=1
    Zcmax(i)=2
    ! 1s
    nstar(nlm)=0.
    zeta(nlm)=Z-0.30
    normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
    occ(nlm)=2.
    nlm=nlm+1
   ENDIF
  ELSE IF(ia.gt.36 .and. ia.le.54) THEN
   nSTO(i)=5
   Zcmax(i)=36
   ! 4s,4p
!   nstar(nlm)=2.7
   nstar(nlm)=3.
   zeta(nlm)=Z-10.*1.00-18.*0.85-7.*0.35
   normc(nlm)=16.*norm_1s*((zeta(nlm))**4)*sqrt(zeta(nlm)/20160.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 3d
   nstar(nlm)=2.
   zeta(nlm)=Z-10*1.00-8.*1.00-9.*0.35
   normc(nlm)=8.*norm_1s*((zeta(nlm))**3)*sqrt(zeta(nlm)/360.)
   occ(nlm)=10.
   nlm=nlm+1
   ! 3s,3p
   nstar(nlm)=2.
   zeta(nlm)=Z-2.*1.00-8.*0.85-7*0.35
   normc(nlm)=8.*norm_1s*((zeta(nlm))**3)*sqrt(zeta(nlm)/360.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 2s,2p
   nstar(nlm)=1.
   zeta(nlm)=Z-2.*0.85-7.*0.35
   normc(nlm)=4.*norm_1s*zeta(nlm)*zeta(nlm)*sqrt(zeta(nlm)/12.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 1s
   nstar(nlm)=0
   zeta(nlm)=Z-0.30
   normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
   occ(nlm)=2.
   nlm=nlm+1
  ELSE IF(ia.gt.54 .and. ia.le.86) THEN
   nSTO(i)=7
   Zcmax(i)=68
   ! 5s,5p
   nstar(nlm)=3.
   zeta(nlm)=Z-28.*1.00-32.*0.85-7.*0.35
   normc(nlm)=16.*norm_1s*((zeta(nlm))**4)*sqrt(zeta(nlm)/20160.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 4d,4f
!   nstar(nlm)=2.7
   nstar(nlm)=3.
   zeta(nlm)=Z-28.*1.00-8.*1.00-13.*0.35
   normc(nlm)=16.*norm_1s*((zeta(nlm))**4)*sqrt(zeta(nlm)/20160.)
   occ(nlm)=24.
   nlm=nlm+1
   ! 4s,4p 
!   nstar(nlm)=2.7
   nstar(nlm)=3.
   zeta(nlm)=Z-10.*1.00-18.*0.85-7.*0.35
   normc(nlm)=16.*norm_1s*((zeta(nlm))**4)*sqrt(zeta(nlm)/20160.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 3d
   nstar(nlm)=2.
   zeta(nlm)=Z-10*1.00-8.*1.00-9.*0.35
   normc(nlm)=8.*norm_1s*((zeta(nlm))**3)*sqrt(zeta(nlm)/360.)
   occ(nlm)=10.
   nlm=nlm+1
   ! 3s,3p
   nstar(nlm)=2.
   zeta(nlm)=Z-2.*1.00-8.*0.85-7*0.35
   normc(nlm)=8.*norm_1s*((zeta(nlm))**3)*sqrt(zeta(nlm)/360.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 2s,2p
   nstar(nlm)=1.
   zeta(nlm)=Z-2.*0.85-7.*0.35
   normc(nlm)=4.*norm_1s*zeta(nlm)*zeta(nlm)*sqrt(zeta(nlm)/12.)
   occ(nlm)=8.
   nlm=nlm+1
   ! 1s
   nstar(nlm)=0.
   zeta(nlm)=Z-0.30
   normc(nlm)=2.*norm_1s*zeta(nlm)*sqrt(zeta(nlm))
   occ(nlm)=2.
   nlm=nlm+1
  ELSE
   CALL errore('init_cores','Z higher than Rn, not yet implemented.')
  ENDIF

  Zsc=Zcmax(i)-Zc
  isc=Zsc/2
! WRITE(stdout,'(a,i4,f4.0,f4.0,f4.0,i4,i4,i4,i4)') &
!  "INIT_CORES: core in pseudo ",&
!  ia,Z,Zv,Zc,ic,isc,iSTO(i),nSTO(i)
! WRITE(stdout,'(a,f5.2)') "INIT_CORES: core radius ",rc_all
! DO nlm=iSTO(i)+isc+1,iSTO(i)+nSTO(i)
!  WRITE(stdout,'(a,i4,i4,f10.3,e,f10.3,f10.3)') &
!   "INIT_CORES: Slater parameters (type,...) ",&
!   i,nlm,occ(nlm),normc(nlm),nstar(nlm),zeta(nlm)
! ENDDO

 ENDDO

 RETURN

 END SUBROUTINE init_cores

 SUBROUTINE do_core(i,r,rho_atc)
! PURE FUNCTION do_core(i,r) RESULT(sum3)

 USE kinds,          ONLY : DP
 USE uspp_param,     ONLY : upf
 USE io_global,      ONLY : stdout, ionode
 USE ions_base,      ONLY : atom_label => atm

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: i
 REAL (DP), INTENT(IN) :: r
 REAL (DP), INTENT(OUT) :: rho_atc

 INTEGER, EXTERNAL :: atomic_number

 INTEGER   :: ia,nlm,isc
 REAL (DP) :: Z,Zc,Zv,Zsc,sum1,sum2,sum3,val,val2,xp
 
! build rho_atc

 ia=atomic_number(atom_label(i))
 Z=REAL(ia)
 Zv=upf(i)%zp
 Zc=Z-Zv
 Zsc=Zcmax(i)-Zc
 isc=Zsc/2
!
! WRITE(stdout,'(a,i,f4.0,f4.0,f4.0,i4,i4,i4)') &
!  "DO_CORE: core in pseudo ",&
!  ia,Z,Zv,Zc,isc,iSTO(i),nSTO(i)
! DO nlm=iSTO(i)+isc+1,iSTO(i)+nSTO(i)
!  WRITE(stdout,'(a,i4,i4,f10.3,e,f10.3,f10.3)') &
!   "DO_CORE: Slater parameters (type,...) ",&
!   i,nlm,occ(nlm),normc(nlm),nstar(nlm),zeta(nlm)
! ENDDO
!
 sum3=0.d0
 DO nlm=iSTO(i)+isc+1,iSTO(i)+nSTO(i)
  sum1=normc(nlm)
  sum2=nstar(nlm)
  xp = -r*zeta(nlm)
  val = sum1*(r**sum2)*exp(xp)
  val2=val*val
  sum3=sum3+occ(nlm)*val2
 ENDDO
 rho_atc=sum3
! do_core=sum3
! write(stdout,'(i,f,f)') i,r,sum3
 
 RETURN

 END SUBROUTINE do_core
! END FUNCTION do_core

END MODULE adduscore

