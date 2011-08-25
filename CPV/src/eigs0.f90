!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
   subroutine eigs0( ei, nudx, tprint, nspin, nupdwn, iupdwn, lf, f, nx, lambda, nlam, desc )
!-----------------------------------------------------------------------
!     computes eigenvalues (wr) of the real symmetric matrix lambda
!     Note that lambda as calculated is multiplied by occupation numbers
!     so empty states yield zero. Eigenvalues are printed out in eV
!
      use kinds,             only : DP
      use io_global,         only : stdout
      use constants,         only : autoev
      use dspev_module,      only : dspev_drv, pdspev_drv
      USE sic_module,        only : self_interaction
      USE descriptors,       ONLY : la_descriptor
      USE mp,                only : mp_sum, mp_bcast
      USE mp_global,         only : intra_bgrp_comm, root_bgrp, me_bgrp

      implicit none
! input
      logical, intent(in) :: tprint, lf
      integer, intent(in) :: nspin, nx, nudx, nupdwn(nspin), iupdwn(nspin), nlam
      type(la_descriptor), intent(in) :: desc( 2 )
      real(DP), intent(in) :: lambda( nlam, nlam, nspin ), f( nx )
      real(DP), intent(out) :: ei( nudx, nspin )
! local variables
      real(DP), allocatable :: ap(:), wr(:)
      real(DP) zr(1)
      integer :: iss, j, i, ierr, k, n, ndim, nspin_eig, npaired
      INTEGER :: ir, ic, nr, nc, nrl, nrlx, comm, np, me
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

         IF( nudx < nupdwn(iss) ) THEN 
            WRITE( msg, 100 ) nudx, SIZE( ei, 1 ), nupdwn(iss)
100         FORMAT( ' wrong dimension array ei = ', 3I10 )
            CALL errore( ' eigs0 ', msg, 1 )
         END IF

         IF( tsic ) THEN
            n = npaired
         ELSE
            n = nupdwn(iss)
         END IF

         allocate( wr( n ) )

         IF( desc( iss )%active_node > 0 ) THEN

            np = desc( iss )%npc * desc( iss )%npr

            IF( np > 1 ) THEN

               !  matrix is distributed

               CALL qe_pdsyevd( .false., n, desc(iss), lambda(1,1,iss), nlam, wr )

            ELSE

               !  matrix is not distributed

               allocate( ap( n * ( n + 1 ) / 2 ) )

               k = 0
               do i = 1, n
                  do j = i, n
                     k = k + 1
                     ap( k ) = lambda( j, i, iss )
                  end do
               end do

               CALL dspev_drv( 'N', 'L', n, ap, wr, zr, 1 )

               deallocate( ap )

            END IF

         END IF

         call mp_bcast( wr, root_bgrp, intra_bgrp_comm )

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
         ei( 1:n, iss ) = wr( 1:n )

         IF( tsic ) THEN
            !
            !  store unpaired state
            !
            ei( 1:n,       1 ) = ei( 1:n, 1 ) / 2.0d0
            ei( nupdwn(1), 1 ) = 0.0d0
            if( desc( iss )%active_node > 0 ) then
               IF( desc( iss )%myc == desc( iss )%myr ) THEN
                  ir = desc( iss )%ir
                  nr = desc( iss )%nr
                  IF( nupdwn(1) >= ir .AND. nupdwn(1) < ir + nr ) then
                     ei( nupdwn(1), 1 ) = lambda( nupdwn(1)-ir+1, nupdwn(1)-ir+1, 1 )
                  end if
               END IF
            endif
            call mp_sum( ei( nupdwn(1), 1 ), intra_bgrp_comm )
            !
         END IF

         ! WRITE( stdout,*)  '---- DEBUG ----' ! debug
         ! WRITE( stdout,14) ( wr( i ) * autoev / 2.0d0, i = 1, nupdwn(iss) ) ! debug

         deallocate( wr )

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
      USE mp_global, ONLY: me_bgrp, nproc_bgrp, intra_bgrp_comm
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
            IF( me_bgrp < n ) THEN
              DO i = 1, n
                j = me_bgrp + 1
                DO jl = 1, nrl
                  IF( j >= i ) THEN
                    !   maps (j,i) index to low-tri packed (k) index
                    k = (i-1)*n + j - i*(i-1)/2  
                    aux(k) = gam(jl,i) / f(j)
                  END IF
                  j = j + nproc_bgrp
                END DO
              END DO
            END IF
            CALL mp_sum(aux, intra_bgrp_comm)
          ELSE
            IF( me_bgrp < n ) THEN
              DO i = 1, n
                j = me_bgrp + 1
                DO jl = 1, nrl
                  gam(jl,i) = gam(jl,i) / f(j)
                  j = j + nproc_bgrp
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

   CONTAINS

      DOUBLE PRECISION FUNCTION stepf(x)
      USE kinds
      IMPLICIT NONE
      REAL(DP) :: x
      REAL(DP), PARAMETER :: c=0.5641895835D0
!      REAL(DP), EXTERNAL :: qe_erfc
!     stepf=qe_erfc(x)
      stepf=1.d0/(exp(min(x,100.d0))+1.d0)
      END FUNCTION stepf 


   END SUBROUTINE fermi_energy_x

!
!
!

!-----------------------------------------------------------------------
   SUBROUTINE cp_eigs_x( nfi, lambdap, lambda, descla )
!-----------------------------------------------------------------------

      USE kinds,             ONLY: DP
      use ensemble_dft,      only: tens
      use electrons_base,    only: nbspx, f, nspin
      use electrons_base,    only: iupdwn, nupdwn, nudx
      use electrons_module,  only: ei
      use io_global,         only: stdout
      USE descriptors,       ONLY: la_descriptor

      IMPLICIT NONE

      INTEGER :: nfi
      REAL(DP) :: lambda( :, :, : ), lambdap( :, :, : )
      TYPE(la_descriptor), INTENT(IN) :: descla( : )

      if( .not. tens ) then
         call eigs0( ei, nudx, .false. , nspin, nupdwn, iupdwn, .true. , f, nbspx, lambda, SIZE(lambda,1), descla )
      else
         call eigs0( ei, nudx, .false. , nspin, nupdwn, iupdwn, .false. , f, nbspx, lambdap, SIZE(lambdap,1), descla )
      endif

      RETURN
   END SUBROUTINE cp_eigs_x
