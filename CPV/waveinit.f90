!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  Last modified: Sat Dec  4 07:03:35 MET 1999
!  ----------------------------------------------
!  BEGIN manual

      MODULE wave_init

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE pw_rand_init(nbeg,cm,c0,gv)
!  SUBROUTINE calphi(pslm,philm,gv,eigr,nchan)
!  SUBROUTINE calpslm(pslm,gv,nchan,ns,np,nd)
!  SUBROUTINE pw_atomic_init(nbeg,cm,c0,gv,eigr)
!  ----------------------------------------------
!  END manual

! ...   include modules
        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        INTEGER, PARAMETER :: maxchan = 4, maxgh = 10

        PUBLIC :: pw_atomic_init

! ...   end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE pw_rand_init(nbeg, cm, c0, wfill, ce, wempt, gv, kp)

!  this routine sets the initial wavefunctions at random
!  ----------------------------------------------

! ... declare modules
      USE wave_functions, ONLY: gram, dft_kinetic_energy
      USE wave_types, ONLY: wave_descriptor
      USE mp, ONLY: mp_sum
      USE mp_wave, ONLY: splitwf
      USE mp_global, ONLY: mpime, nproc, root
      USE reciprocal_vectors, ONLY: ig_l2g, ngw, ngwt
      USE brillouin, ONLY: kpoints
      USE cp_types, ONLY: recvecs
      USE io_base, ONLY: stdout
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

! ... declare module-scope variables

! ... declare subroutine arguments
      COMPLEX(dbl), INTENT(OUT) :: cm(:,:,:,:), c0(:,:,:,:), ce(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
      TYPE (recvecs), INTENT(IN) :: gv
      TYPE (kpoints), INTENT(IN) :: kp
      INTEGER, INTENT(IN) :: nbeg

      REAL(dbl) :: rranf
      EXTERNAL rranf

! ... declare other variables
      INTEGER :: ig_local
      INTEGER :: ntest, i, ig, j, ib, ik, ispin, nspin
      REAL(dbl) ::  rranf1, rranf2, ampre, rsum
      REAL(dbl) ::  rc0rc0, anorm
      COMPLEX(dbl) :: ctmp
      COMPLEX(dbl), ALLOCATABLE :: pwt( : )

! ... end of declarations
!  ----------------------------------------------

!
! ... Check array dimensions
!
      IF( ( SIZE( cm, 1 ) /= SIZE( c0, 1 ) ) .OR. &
          ( SIZE( cm, 2 ) /= SIZE( c0, 2 ) ) .OR. &
          ( SIZE( cm, 3 ) /= SIZE( c0, 3 ) ) .OR. &
          ( SIZE( cm, 4 ) /= SIZE( c0, 4 ) ) ) &
        CALL errore(' pw_rand_init ', ' wrong dimensions ', 1)

      IF( .NOT. force_pairing ) THEN
        IF( ( SIZE( ce, 1 ) /= SIZE( c0, 1 ) ) .OR. &
            ( SIZE( ce, 3 ) /= SIZE( c0, 3 ) ) .OR. &
            ( SIZE( ce, 4 ) /= SIZE( c0, 4 ) ) ) &
          CALL errore(' pw_rand_init ', ' wrong dimensions ', 3)
      END IF


! ... Reset them to zero
!
      cm = 0.0d0
      c0 = 0.0d0
      ce = 0.0d0

! ... initialize the wave functions in such a way that the values
! ... of the components are independent on the number of processors
!
      nspin = SIZE( cm, 4)
      IF( force_pairing ) nspin = 1

      IF( nbeg <  1) THEN

        ampre = 0.01d0
        ALLOCATE( pwt( ngwt ) )

        DO ispin = 1, nspin
          DO ik = 1, SIZE( cm, 3)
            ntest = ngwt / 4
            IF( ntest < wfill%nbt( ispin ) ) THEN 
              ntest = ngwt
            END IF
! ...       assign random values to wave functions
            DO ib = 1, SIZE( cm, 2 )
              pwt( : ) = 0.0d0
              DO ig = 3, ntest
                rranf1 = 0.5d0 - rranf()
                rranf2 = rranf()
                pwt( ig ) = ampre * DCMPLX(rranf1, rranf2)
              END DO
              CALL splitwf ( cm( :, ib, ik, ispin ), pwt, ngw, ig_l2g, mpime, nproc, 0 )
            END DO
            IF ( .NOT. kp%gamma_only ) THEN
! ..  .       set to zero all elements outside the cutoff sphere
              DO ib = 1, SIZE( cm, 2 )
                cm(:, ib, ik, ispin) = cm(:, ib, ik, ispin) * gv%kg_mask_l(:,ik)
              END DO
            END IF
          END DO
          DO ik = 1, SIZE( cm, 3 )
            IF ( gv%gzero ) THEN
              cm( 1, :, ik, ispin ) = (0.0d0, 0.0d0)
            END IF
          END DO
        END DO

        DEALLOCATE( pwt )
    
! DEBUG
!        ctmp = SUM( cm )
!        CALL mp_sum( ctmp ) 
!        WRITE( stdout,*) ' *** Wave init check ', ctmp
! DEBUG

!
! ...   orthonormalize wave functions
!
        CALL gram( cm, wfill )

        c0 = cm

! DEBUG
!        ctmp = SUM( cm )
!        CALL mp_sum( ctmp ) 
!        WRITE( stdout,*) ' *** Wave init check ', ctmp
! DEBUG

      END IF

      RETURN
      END SUBROUTINE pw_rand_init

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE calphi(pslm, philm, gv, eigr, nchan, ik)

!  this routine computes the atomic wavefunctions from pseudowf
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: nax, nsp, na, nat
      USE gvecw, ONLY: gkcut
      USE cp_types, ONLY: recvecs
      USE spherical_harmonics, ONLY: gkl

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: ik
      REAL(dbl)     pslm(:,:,:)
      COMPLEX(dbl) philm(:,:,:),eigr(:,:)
      TYPE (recvecs) gv
      INTEGER nchan(nsp)

! ... declare other variables
      INTEGER ig,ia,is,igh,i,isa
      INTEGER l,ll,m
      REAL(dbl)  arg, gg, gcutby4

! ... end of declarations
!  ----------------------------------------------

        igh=1
        isa = 0
        philm = 0.0d0
        gcutby4 = gkcut / 4.0d0

        DO is = 1, nsp
          ll = 1
          IF(nchan(is).GE.1) THEN
            DO ia=1,na(is)
              isa = isa + 1
              DO ig=1,gv%ngw_l
                IF(  gv%khg_l(ig,ik) <  gcutby4 ) THEN
                  philm(ig,isa,igh) = pslm(ig,ll,is) * eigr(ig,isa)
                END IF
              END DO
            END DO
          ELSE
            isa = isa + na(is)
          END IF
        END DO

        DO m = 1, 3
          igh=igh+1
          isa = 0
          DO is = 1, nsp
            ll=2
            IF(nchan(is).GE.2) THEN
              DO ia=1,na(is)
                isa = isa + 1
                DO ig=1,gv%ngw_l
                  gg = 0.0d0
                  IF(  gv%khg_l(ig,ik) < gcutby4 ) THEN
                    gg = gv%kgx_l(m,ig,ik) * eigr(ig,isa)
                  END IF
                  philm(ig,isa,igh) = CMPLX(0.d0,pslm(ig,ll,is)) * gg
                END DO
              END DO
            ELSE
              isa = isa + na(is)
            END IF
          END DO
        END DO

        DO m=1,5
          igh=igh+1
          isa = 0
          DO is=1,nsp
            ll=3
            IF(nchan(is).GE.3) THEN
              DO ia=1,na(is)
                isa = isa + 1
                DO ig=1,gv%ngw_l
                  gg = 0.0d0
                  IF(  gv%khg_l(ig,ik) .GT. 1.0d-12 .AND. gv%khg_l(ig,ik) < gcutby4 ) THEN
                    gg = gkl(gv%kgx_l(1,ig,ik), gv%kgx_l(2,ig,ik),gv%kgx_l(3,ig,ik), &
                      gv%khg_l(ig,ik),m)
                  END IF
                  philm(ig,isa,igh) = - pslm(ig,ll,is) * gg * eigr(ig,isa)
                END DO
              END DO
            ELSE
              isa = isa + na(is)
            END IF
          END DO
        END DO

      RETURN
      END SUBROUTINE calphi

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE calpslm(pslm,gv,nchan,ns,np,nd,ik)

!  this routine computes the temporary arrays pslm
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: nsp
      USE pseudopotential, ONLY: pseudo_wave_info
      USE cell_base, ONLY: tpiba
      USE bessel_functions, ONLY: bessel2
      USE cp_types, ONLY: recvecs
      USE parameters, ONLY: ndmx

      IMPLICIT  NONE

! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: ik
      REAL(dbl)  pslm(:,:,:)
      INTEGER nchan(:),ns(:),np(:),nd(:)
      TYPE (recvecs) gv

! ... declare other variables
      REAL(dbl)  xg, xrg, xrgm1, dum, gg
      INTEGER ig,is,l1,l2,mmax,igs,me,l,ll,ir
      INTEGER, ALLOCATABLE :: tindl(:)
      REAL(dbl), ALLOCATABLE :: oc(:,:)
      REAL(dbl),  ALLOCATABLE :: fint(:,:)
      REAL(dbl),  ALLOCATABLE :: func(:,:)
      REAL(dbl),  ALLOCATABLE :: dx(:),rw(:,:), rps(:,:,:)
      INTEGER, ALLOCATABLE :: mesh(:)

! ... end of declarations
!  ----------------------------------------------

      ALLOCATE(dx(nsp))
      ALLOCATE(mesh(nsp))
      ALLOCATE(oc(maxchan,nsp))
      ALLOCATE(rw(ndmx,nsp))
      ALLOCATE(rps(ndmx,nsp,maxchan))
      ALLOCATE(tindl(maxchan))

      DO l=1,maxchan
        tindl(l)=l
      END DO
      CALL pseudo_wave_info(oc,nchan,mesh,dx,rw,rps)

      DO is=1,nsp

        ns(is) = 0
        np(is) = 0
        nd(is) = 0
        IF(nchan(is).LE.0) CYCLE
        IF(nchan(is).GT.0) ns(is) = oc(1,is)
        IF(nchan(is).GT.1) np(is) = oc(2,is)
        IF(nchan(is).GT.2) nd(is) = oc(3,is)

        IF(nchan(is).GT.maxchan) THEN
          CALL errore(' module wave_init ', ' nchan out of range ',nchan(is))
        END IF

        mmax=mesh(is)
        ALLOCATE(fint(mmax,maxchan))
        ALLOCATE(func(mmax,maxchan))

        DO l=1,nchan(is)
          func(:,l) = rw(:,is)**2 * rps(:,is,l)
        END DO

        DO ig = 1, SIZE( pslm, 1 )

          gg = gv%khg_l(ig,ik)

! ...     G=0 (Only if L=1, since otherwise the radial Bessel function JL=0)
          IF(gg .LT. 1.0d-12) THEN
            DO l = 1, nchan(is)
              IF( l.EQ.1 ) THEN
                fint(:,l)=func(:,l)
                call simpson_fpmd(mmax,fint(1,l),dx(is),pslm(ig,l,is))
              ELSE
                pslm(ig,l,is)=0.d0
              END IF
            END DO
          ELSE
            xg=sqrt(gg)*tpiba
            CALL bessel2(xg,rw(:,is),fint,nchan(is),tindl,mmax)
            DO l=1,nchan(is)
              fint(1:mmax,l) = fint(1:mmax,l) * func(1:mmax,l)
              call simpson_fpmd(mmax,fint(1,l),dx(is),pslm(ig,l,is))
            END DO
          END IF
        END DO

        DEALLOCATE(fint, func)

      END DO

      DEALLOCATE(dx, mesh, rw, rps, tindl, oc)

      RETURN
      END SUBROUTINE calpslm

!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE pw_atomic_init(nbeg, cm, c0, wfill, ce, wempt, gv, kp, eigr)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: nat, nsp, na
      USE wave_functions, ONLY: gram, rande
      USE wave_types, ONLY: wave_descriptor
      USE control_flags, ONLY: tatomicwfc
      USE gvecw, ONLY: gkcut
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE brillouin, ONLY: kpoints
      USE cp_types, ONLY: recvecs
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: nbeg
      TYPE (recvecs), INTENT(IN) :: gv
      TYPE (kpoints), INTENT(IN) :: kp
      COMPLEX(dbl), INTENT(OUT) :: cm(:,:,:,:), c0(:,:,:,:), ce(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: wfill, wempt
      COMPLEX(dbl), INTENT(IN) :: eigr(:,:)
      INTEGER ns(nsp),np(nsp),nd(nsp)

      COMPLEX(dbl), ALLOCATABLE :: philm(:,:,:)
      REAL(dbl), ALLOCATABLE :: pslm(:,:,:)

      REAL(dbl) :: rranf
      EXTERNAL rranf

! ... declare other variables
      INTEGER ii, i, ik, ig, ia, is, igh, ngw, n, isa, nspin, ispin
      INTEGER ie, nga
      INTEGER at_occ(9,nsp)
      INTEGER nchan(nsp)
      REAL(dbl) :: rranf1, rranf2, gcutby4

! ... end of declarations
!  ----------------------------------------------

      CALL pw_rand_init(nbeg, cm, c0, wfill, ce, wempt, gv, kp)

      IF( nbeg >= 0 ) THEN
        RETURN
      END IF

      IF( .NOT. tatomicwfc ) THEN

        IF( ionode ) WRITE( stdout, fmt = '(/,3X, &
          & "Wave Initialization: random initial wave-functions" &
        & )' )

      ELSE

        nspin   = wfill%nspin
        IF( force_pairing ) nspin = 1

        gcutby4 = gkcut / 4.0d0


        IF( ionode ) WRITE( stdout, fmt = '(/,3X, &
          & "Wave Initialization: Initial wave-functions ", &
          & "from superposition of atomic orbitals" &
        & )' )

        SPIN: DO ispin = 1, nspin

          KAPPA: DO ik = 1, wfill%nkl
    
            ngw = wfill%ngwl
            n   = wfill%nbl( ispin )
            ALLOCATE(philm(ngw,nat,maxgh))
            ALLOCATE(pslm(ngw,maxchan,nsp))
    
            CALL calpslm(pslm,gv,nchan,ns,np,nd,ik)
            CALL calphi(pslm,philm,gv,eigr,nchan,ik)
    
            at_occ = 0
            DO is = 1, nsp
              at_occ(1,is) = ns(is)
              DO igh = 2, 4 
                IF( np(is) >= 2 ) THEN
                  at_occ(igh, is) = 2
                  np(is) = np(is) - 2
                ELSE IF( np(is) >= 1 ) THEN
                  at_occ(igh, is) = 1
                  np(is) = np(is) - 1
                ELSE
                  at_occ(igh,is) = 0
                END IF
              END DO
              DO igh = 5, 9
                IF( nd(is) >= 2 ) THEN
                  at_occ(igh, is) = 2
                  nd(is) = nd(is) - 2
                ELSE IF( nd(is) >= 1 ) THEN
                  at_occ(igh, is) = 1
                  nd(is) = nd(is) - 1
                ELSE
                  at_occ(igh,is) = 0
                END IF
              END DO
            END DO
    
            i  = 1
            ie = 0
            cm( :, :, ik, ispin) = 0.0d0
            DO igh = 1, 9
              isa = 0
              DO is = 1, nsp
                DO ia = 1, na(is)
                  isa = isa + 1
                  IF( nspin < 2 ) THEN
                    DO ii = 1, at_occ(igh,is)
                      cm( :, i, ik, ispin) = cm( :, i, ik, ispin) + 0.5d0 * philm(:,isa,igh)
                      ie = ie + 1
                      IF( MOD(ie, 2) == 0) i = i + 1
                    END DO
                  ELSE
                    DO ii = 1, at_occ(igh,is)
                      cm( :, i, ik, ispin) = philm(:,isa,igh)
                      ie = ie + 1
                      i  = i  + 1
                    END DO
                  END IF
                END DO
              END DO
            END DO
            DO ii = i, wfill%nbl( ispin )
              DO ig = 2, wfill%ngwl
                rranf1 = 0.5 - rranf()
                rranf2 = rranf()
                cm( ig, ii, ik, ispin ) =  0.001d0 * CMPLX(rranf1,rranf2)
              END DO
            END DO
            IF( .NOT. kp%gamma_only ) THEN
              DO i = 1, wfill%nbl( ispin )
                cm( :, i, ik, ispin) = cm( :, i, ik, ispin) * gv%kg_mask_l( :, ik )
              END DO
            END IF
            DEALLOCATE(philm, pslm)
          END DO KAPPA
        END DO SPIN
  
        IF( force_pairing ) THEN
          CALL rande( 1, cm(:,:,:,1), wfill, 0.01d0)
        ELSE
          CALL rande(cm, wfill, 0.01d0)
        END IF

        SPIN2: DO ispin = 1, nspin
          KAPPA2: DO ik = 1, wfill%nkl
            DO i = 1, wfill%nbl( ispin )
              DO ig = 1, wfill%ngwl
                IF( gv%khg_l(ig,ik) > gcutby4 ) THEN
                  cm( ig, i, ik, ispin ) = 0.0d0
                END IF
              END DO
            END DO
          END DO KAPPA2
        END DO SPIN2

        CALL gram( cm, wfill )

      END IF

      RETURN
   END SUBROUTINE pw_atomic_init

!  ----------------------------------------------
!  ----------------------------------------------

      END MODULE wave_init

