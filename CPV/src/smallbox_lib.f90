!
! Copyright (C) 2002-2011 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      SUBROUTINE initbox ( tau0, alat, at, ainv, taub, irb, iabox, nabox )
!-----------------------------------------------------------------------
      !! Sets the indexes \(\text{irb}\) and positions \(\text{taub}\) for
      !! the small boxes around atoms.
!
      USE kinds,                    ONLY: DP
      USE ions_base,                ONLY: nat, ityp
      USE uspp_param,               ONLY: upf
      USE control_flags,            ONLY: iverbosity
      USE io_global,                ONLY: stdout
      USE mp_global,                ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, my_bgrp_id, nbgrp
      USE fft_base,                 ONLY: dfftb, dfftp, fft_type_descriptor
      USE fft_smallbox_type,        ONLY: fft_box_set

      IMPLICIT NONE
! input
      REAL(DP), INTENT(in)  :: tau0(3,nat), at(3,3), ainv(3,3), alat
! output
      INTEGER,  INTENT(out) :: irb(3,nat)
      REAL(DP), INTENT(out) :: taub(3,nat)
      INTEGER,  INTENT(out) :: iabox(nat), nabox
! local
      REAL(DP) :: x(3), xmod
      INTEGER  :: nr(3), nrb(3), xint, ia, i, is
!
      IF ( dfftb%nr1 < 1) CALL errore ('initbox', 'incorrect value for box grid dimensions', 1)
      IF ( dfftb%nr2 < 1) CALL errore ('initbox', 'incorrect value for box grid dimensions', 2)
      IF ( dfftb%nr3 < 1) CALL errore ('initbox', 'incorrect value for box grid dimensions', 3)

      nr (1)=dfftp%nr1
      nr (2)=dfftp%nr2
      nr (3)=dfftp%nr3
      nrb(1)=dfftb%nr1
      nrb(2)=dfftb%nr2
      nrb(3)=dfftb%nr3
!
      DO ia=1,nat
!
            DO i=1,3
!
! bring atomic positions to crystal axis
!
               x(i) = ainv(i,1)*tau0(1,ia) +                         &
     &                ainv(i,2)*tau0(2,ia) +                         &
     &                ainv(i,3)*tau0(3,ia)
!
! bring x in the range between 0 and 1
!
               x(i) = MOD(x(i),1.d0)
               IF (x(i).LT.0.d0) x(i)=x(i)+1.d0
!
! case of nrb(i) even
!
               IF (MOD(nrb(i),2).EQ.0) THEN
!
! find irb = index of the grid point at the corner of the small box
!           (the indices of the small box run from irb to irb+nrb-1)
!
                  xint=INT(x(i)*nr(i))
                  irb (i,ia)=xint+1-nrb(i)/2+1
                  IF(irb(i,ia).LT.1) irb(i,ia)=irb(i,ia)+nr(i)
!
! x(i) are the atomic positions in crystal coordinates, where the
! "crystal lattice" is the small box lattice and the origin is at
! the corner of the small box. Used to calculate phases exp(iG*taub)
!
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+nrb(i)/2-1)/nr(i)
               ELSE
!
! case of nrb(i) odd - see above for comments
!
                  xint=NINT(x(i)*nr(i))
                  irb (i,ia)=xint+1-(nrb(i)-1)/2
                  IF(irb(i,ia).LT.1) irb(i,ia)=irb(i,ia)+nr(i)
                  xmod=x(i)*nr(i)-xint
                  x(i)=(xmod+(nrb(i)-1)/2)/nr(i)
               END IF
            END DO
!
! bring back taub in cartesian coordinates
!
            DO i=1,3
               taub(i,ia)=(x(1)*at(i,1) + x(2)*at(i,2) + x(3)*at(i,3))*alat
            END DO
      END DO

      ! initialize FFT descriptor

      CALL fft_box_set( dfftb, nat, irb, dfftp )

      ! build the local list of atom 
      nabox = 0
      DO ia = 1, nat
         IF( .NOT. upf(ityp(ia))%tvanp )  CYCLE
#if defined(__MPI)
         IF(  ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
            CYCLE
         END IF
#endif
         nabox = nabox + 1
         iabox( nabox ) = ia
      END DO

      IF( iverbosity > 1 ) THEN
           DO ia=1,nat
              WRITE( stdout,2000) ia, (irb(i,ia),i=1,3)
           END DO
2000       FORMAT(2x, 'atom= ', i3, ' irb1= ', i3, ' irb2= ', i3, ' irb3= ', i3)
      ENDIF

#if defined(__MPI)
      ! 
      ! for processor that do not call fft on the box
      ! artificially start the clock
      ! 
      CALL start_clock( 'fftb' )
      CALL stop_clock( 'fftb' )
      !
#endif
!
      RETURN
   END SUBROUTINE initbox
!
!-----------------------------------------------------------------------
      SUBROUTINE phbox( taub, iverbosity, eigrb )
!-----------------------------------------------------------------------
      !! Calculates the phase factors for the g's of the little box:
      !! \[ \text{eigrt}=\exp(-i\ g\ \text{tau})\ . \]
      !! Uses the same logic for fast calculation as in \(\texttt{phfac}\).
!
      USE kinds,         only: DP
      use io_global,     only: stdout
      use ions_base,     only: nsp, na, nat
      use cell_base,     only: r_to_s
      use cp_interfaces, only: phfacs
      use small_box,     only: bgb, alatb
      use smallbox_gvec, only: ngb, mill_b
      use fft_base,      only: dfftb
!                 
      IMPLICIT NONE    
      REAL(DP), INTENT(IN)    :: taub(3,nat)
      COMPLEX(DP), INTENT(OUT) :: eigrb(ngb,nat)
      INTEGER, INTENT(IN) :: iverbosity
! local           
      REAL(DP)    :: ainvb(3,3)
      integer :: i, ia, ig
      complex(dp), allocatable:: ei1b(:,:), ei2b(:,:), ei3b(:,:)
      real(dp), allocatable :: taus(:,:)
!
      allocate(ei1b(-dfftb%nr1:dfftb%nr1,nat))
      allocate(ei2b(-dfftb%nr2:dfftb%nr2,nat))
      allocate(ei3b(-dfftb%nr3:dfftb%nr3,nat))
      allocate( taus( 3, nat ) )
!
      if(iverbosity > 2) then
         WRITE( stdout,*) ' phbox: taub '
         WRITE( stdout,*) ( (taub(i,ia), i=1, 3 ), ia=1, nat )
      endif

      ainvb(1,:) = bgb(:,1)/alatb
      ainvb(2,:) = bgb(:,2)/alatb
      ainvb(3,:) = bgb(:,3)/alatb

      CALL r_to_s( taub, taus, nat, ainvb )
      CALL phfacs( ei1b, ei2b, ei3b, eigrb, mill_b, taus, dfftb%nr1,dfftb%nr2,dfftb%nr3, nat )
!
      if(iverbosity > 2) then
         WRITE( stdout,*)
         do ia=1,nat
            WRITE( stdout,'(33x,a,i4)') ' ei1b, ei2b, ei3b (ia)',ia
            do ig=1,4
               WRITE( stdout,'(6f9.4)') ei1b(ig,ia),ei2b(ig,ia),ei3b(ig,ia)
            end do
            WRITE( stdout,*)
         end do
      endif
!
      deallocate(ei3b)
      deallocate(ei2b)
      deallocate(ei1b)
      deallocate( taus )
!
      RETURN
      END SUBROUTINE phbox

