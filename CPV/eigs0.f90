!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
   subroutine eigs0( ei, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nudx, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv
      USE sic_module,        only : self_interaction
      USE cp_main_variables, only : nlax, nlam, la_proc
      USE descriptors,       ONLY : nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , descla_siz_
      USE mp,                only : mp_sum
      USE mp_global,         only : intra_image_comm

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin)
      integer, intent(in) :: desc( descla_siz_ , 2 )
      real(DP), intent(in) :: lambda( nlam, nlam, nspin ), f( nx )
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: lambdar(:), wr(:)
      real(DP) zr(1)
      integer :: iss, j, i, ierr, k, n, nspin_eig, npaired
      INTEGER :: ir, ic, nr, nc
      logical :: tsic
      CHARACTER(LEN=80) :: msg
!
      tsic = ( ABS( self_interaction) /= 0 )

      IF( tsic ) THEN
         nspin_eig = 1
         npaired   = nupdwn(2)
      ELSE
         nspin_eig = nspin
         npaired   = 0
      END IF

      do iss = 1, nspin_eig

         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF

         allocate( lambdar( n * ( n + 1 ) / 2 ), wr(n) )

         lambdar = 0.0d0

         k = 0

         !  lambda is distributed across processors

         IF( la_proc ) THEN
            ir = desc( ilar_ , iss )
            ic = desc( ilac_ , iss )
            nr = desc( nlar_ , iss )
            nc = desc( nlac_ , iss )
            do i = 1, n
               do j = i, n
                  k = k + 1
                  IF( i >= ic .AND. i - ic + 1 <= nc ) THEN
                     IF( j >= ir .AND. j - ir + 1 <= nr ) THEN
                        lambdar( k ) = lambda( j - ir + 1, i - ic + 1, iss )
                     END IF
                  END IF
               end do
            end do
         END IF

         CALL mp_sum( lambdar, intra_image_comm )

         CALL dspev_drv( 'N', 'L', n, lambdar, wr, zr, 1 )

         if( lf ) then
            do i = 1, n
               if ( f(iupdwn(iss)-1+i).gt.1.e-6) then
                  wr(i)=wr(i)/f(iupdwn(iss)-1+i)
               else
                  wr(i)=wr(i)/2.0d0 * nspin  ! fake occupation factor to print empty states
               end if
            end do
         end if
         !
         !     store eigenvalues
         !
         IF( nudx < nupdwn(iss) ) THEN 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
         END IF

         ei( 1:n, iss ) = wr( 1:n )

         IF( tsic ) THEN
            !
            !  store unpaired state
            !
            ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
            ei( nupdwn(1), 1 ) = lambda( nupdwn(1), nupdwn(1), 1 )
            !
         END IF

         ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
         ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug

         deallocate( lambdar, wr )

      end do
      !
      !
      do iss = 1, nspin

         IF( tsic .AND. iss == 2 ) THEN
            ei( 1:npaired, 2 ) = ei( 1:npaired, 1 )
         END IF

         IF( tprint ) THEN
            !
            !     print out eigenvalues
            !
            WRITE( stdout,12) 0.d0, 0.d0, 0.d0
            WRITE( stdout,14) ( ei( i, iss ) * autoev, i = 1, nupdwn(iss) )

         ENDIF

      end do

      IF( tprint ) WRITE( stdout,*)

   12 format(//' eigenvalues at k-point: ',3f6.3)
   14 format(10f8.2)
!
      return
   end subroutine eigs0


!-----------------------------------------------------------------------
   SUBROUTINE rpackgam_x( gam, f, aux )
!-----------------------------------------------------------------------
      USE kinds,            ONLY: DP
      USE mp_global, ONLY: me_image, nproc_image, intra_image_comm
      USE mp, ONLY: mp_sum
      IMPLICIT NONE
          REAL(DP), INTENT(INOUT)  :: gam(:,:)
          REAL(DP), INTENT(OUT), OPTIONAL :: aux(:)
          REAL(DP), INTENT(IN)  :: f(:)
          INTEGER n, nrl, i, j, k, jl
          nrl = SIZE(gam, 1)
          n   = SIZE(gam, 2)
          IF( PRESENT( aux ) ) THEN
            aux = 0.0d0
            IF( me_image < n ) THEN
              DO i = 1, n
                j = me_image + 1
                DO jl = 1, nrl
                  IF( j >= i ) THEN
                    !   maps (j,i) index to low-tri packed (k) index
                    k = (i-1)*n + j - i*(i-1)/2  
                    aux(k) = gam(jl,i) / f(j)
                  END IF
                  j = j + nproc_image
                END DO
              END DO
            END IF
            CALL mp_sum(aux, intra_image_comm)
          ELSE
            IF( me_image < n ) THEN
              DO i = 1, n
                j = me_image + 1
                DO jl = 1, nrl
                  gam(jl,i) = gam(jl,i) / f(j)
                  j = j + nproc_image
                END DO
              END DO
            END IF
          END IF
      RETURN
   END SUBROUTINE rpackgam_x



!-----------------------------------------------------------------------
   SUBROUTINE fermi_energy_x(eig, occ, wke, ef, qtot, temp, sume)
!-----------------------------------------------------------------------

!  this routine computes Fermi energy and weights of occupied states
!  using an improved Gaussian-smearing method
!  refs: C.L.Fu and K.M.Ho, Phys.Rev. B28, 5480 (1983)
!        M.Methfessel and A.T.Paxton Phys.Rev. B40 (15 aug. 89).
!
!  taken from APW code by J. Soler and A. Williams (jk+ss)
!  added computation of occupation numbers without k-point weight

      USE kinds,          ONLY: DP
      USE io_global,      ONLY: stdout
      USE electrons_base, ONLY: nspin, iupdwn

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP) :: occ(:)
      REAL(DP) ef, qtot, temp, sume
      REAL(DP) eig(:,:), wke(:,:)
      REAL(DP), PARAMETER  :: tol = 1.d-10
      INTEGER,   PARAMETER  :: nitmax = 100
      INTEGER ne, nk

! ... declare functions
      REAL(DP) stepf

! ... declare other variables
      REAL(DP) sumq,emin,emax,fac,t,drange
      INTEGER ik,ispin,ie,iter

!  end of declarations
!  ----------------------------------------------

      nk    = 1
      ne    = SIZE( occ, 1)
      sumq=0.d0
      sume=0.d0
      emin=eig(1,1)
      emax=eig(1,1)
      fac=2.d0
      IF (nspin.EQ.2) fac=1.d0

      DO ik=1,nk
        DO ispin=1,nspin
          DO ie=1,ne
            wke(ie,ispin) = fac
            occ(ie+iupdwn(ispin)-1) = fac
            sumq=sumq+wke(ie,ispin)
            sume=sume+wke(ie,ispin)*eig(ie,ispin)
            emin=MIN(emin,eig(ie,ispin))
            emax=MAX(emax,eig(ie,ispin))
          END DO
        END DO
      END DO
      ef=emax
      IF (abs(sumq-qtot).LT.tol) RETURN
      IF (sumq.LT.qtot) THEN
        WRITE( stdout,*) 'FERMIE: NOT ENOUGH STATES'
        WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
        STOP
      END IF
      t = MAX(temp,1.d-6)
      drange = t * SQRT( - LOG( tol*.01d0) )
      emin = emin - drange
      emax = emax + drange
      DO iter = 1, nitmax
        ef   = 0.5d0 * (emin+emax)
        sumq = 0.d0
        sume = 0.d0
        DO ik = 1, nk
          DO ispin = 1, nspin
            DO ie = 1, ne
              wke(ie,ispin) = fac / 2.d0 * stepf((eig(ie,ispin)-ef)/t)
              occ(ie+iupdwn(ispin)-1) = fac / 2.d0 * stepf((eig(ie,ispin)-ef)/t)
              sumq = sumq + wke(ie,ispin)
              sume = sume + wke(ie,ispin) * eig(ie,ispin)
            END DO
          END DO
        END DO
        IF (ABS(sumq-qtot).LT.tol) RETURN
        IF (sumq.LE.qtot) emin=ef
        IF (sumq.GE.qtot) emax=ef
      END DO

      WRITE( stdout,*) 'FERMIE: ITERATION HAS NOT CONVERGED.'
      WRITE( stdout,*) 'FERMIE: QTOT,SUMQ=',qtot,sumq
      STOP

   END SUBROUTINE fermi_energy_x

!
!
!

!-----------------------------------------------------------------------
   SUBROUTINE cp_eigs_x( nfi, lambdap, lambda )
!-----------------------------------------------------------------------

      USE kinds,             ONLY: DP
      use ensemble_dft,      only: tens
      use electrons_base,    only: nx => nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use io_global,         only: stdout
      USE cp_main_variables, only: descla

      IMPLICIT NONE

      INTEGER :: nfi
      REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )

      if( .not. tens ) then
         call eigs0( ei, .false. , nspin, nupdwn, iupdwn, .true. , f, nx, lambda, nudx, descla )
      else
         call eigs0( ei, .false. , nspin, nupdwn, iupdwn, .false. , f, nx, lambdap, nudx, descla )
      endif

      WRITE( stdout, * )
 
      RETURN
   END SUBROUTINE cp_eigs_x
