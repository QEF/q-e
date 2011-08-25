!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE berry_phase

        USE io_global,  ONLY : stdout
      
        IMPLICIT NONE

        PRIVATE
        SAVE

        PUBLIC :: berry_setup, berry_closeup

        INTEGER, allocatable, target :: indi_l(:,:)     !  list of G-vec index to be exchanged
        INTEGER, allocatable, target :: sour_indi(:,:)  !  the list of source processors
        INTEGER, allocatable, target :: dest_indi(:,:)  !  the list of destination processors
        INTEGER :: n_indi_rcv(8) = 0   !   number of G-vectors to be received
        INTEGER :: n_indi_snd(8) = 0   !   number of G-vectors to be sent
        INTEGER :: icntix(8)     = 0   !   total number of G-vec to be exchanged
        LOGICAL :: lsetup = .FALSE.

        PUBLIC :: indi_l, sour_indi, dest_indi, n_indi_rcv, n_indi_snd, icntix

      CONTAINS


        SUBROUTINE ln_setup( mill, ngwt )

          !  setups the "C" functions that will manage the 
          !  mapping between the miller index and the g-vector
          !  index 

          INTEGER :: mill( :, : )
          INTEGER :: ngwt
          INTEGER :: ig

          CALL LN_ALLOC( ngwt )
          DO IG = 1, ngwt
            call LN_SET( mill(1,ig), mill(2,ig), mill(3,ig), ig )
          ENDDO
          CALL LN_ACTIVATE

          lsetup = .TRUE.

          RETURN
        END SUBROUTINE ln_setup


        SUBROUTINE ln_closeup
          IF( lsetup ) CALL LN_DEALLOC
          RETURN
        END SUBROUTINE ln_closeup


    SUBROUTINE berry_setup( ngw, mill )

      USE io_global, only: ionode, stdout
      USE mp_global, ONLY: nproc_image, me_image, intra_image_comm
      USE mp, ONLY: mp_max, mp_sum
      USE stick_base, ONLY : sticks_owner
      USE gvect, ONLY: ig_l2g, sortedig_l2g

      IMPLICIT NONE
      integer LN_IND
      integer ig_local
      external ln_ind, ig_local

      integer :: mill(:,:), ngw

      integer in(8)
      integer, allocatable :: icnt_snd(:,:) ! icnt_snd(nproc_image,8)
      integer, allocatable :: icnt_rcv(:,:) ! icnt_rcv(nproc_image,8)
      integer :: i, j, ig, itmp, in_l, ngwt


      IF( ionode ) THEN
         WRITE( stdout, fmt="(3X,'Polarizability using berry phase')" )
      END IF

      allocate( icnt_snd( nproc_image, 8 ) )
      allocate( icnt_rcv( nproc_image, 8 ) )

      !  compute global number of G vectors
      !
      ngwt = ngw
      CALL mp_sum( ngwt )

      CALL ln_setup( mill, ngwt )

      allocate( indi_l( ngw, 8 ) )
      allocate( sour_indi( ngw, 8 ) )
      allocate( dest_indi( ngw, 8 ) )
      n_indi_rcv = 0
      n_indi_snd = 0

      DO IG = 1, ngwt

        !  compute the indexes "in" of the G + 1 vectors

        call indi_of_ig( mill(:,ig), in )

        do i = 1, 8

          if( in(i) > 0 ) then

            !   find out local index in_l corresponding to the global index in(i)

            in_l = ig_local( in(i), ig_l2g, sortedig_l2g, SIZE( ig_l2g ) )

            if( in_l > 0 ) then

              n_indi_snd(i) = n_indi_snd(i) + 1

              !   find out the processor that own the G vector in(i)
              !   and fill in the array of destination procs

              dest_indi( n_indi_snd(i), i ) = sticks_owner( mill(1,ig), mill(2,ig) )

              !   array of index to of G-vecs to be sent to the processor
              !   whose index is stored in dest_indi

              indi_l( n_indi_snd(i), i ) = in_l

            end if

          end if
          if( sticks_owner( mill(1,ig),  mill(2,ig) ) == ( me_image+1 ) ) then
            n_indi_rcv(i) = n_indi_rcv(i) + 1
            if( in(i) > 0 ) then
              sour_indi( n_indi_rcv(i), i ) = sticks_owner( mill( 1 , in(i) ), mill( 2 , in(i) ) )
            else
              sour_indi( n_indi_rcv(i), i ) = -1
            end if
          end if

        end do

      end do
!     calculate dimension for the variable to be allocated
      icnt_snd = 0
      do i = 1,8
        do ig = 1,n_indi_snd(i)
          itmp = dest_indi(ig,i)
          if(itmp.ne.(me_image+1)) then
            icnt_snd(itmp,i) = icnt_snd(itmp,i) + 1
          end if
        end do
      end do
      do i = 1,8
        icntix(i) = 0
        do j=1,nproc_image
          if(icnt_snd(j,i).gt.icntix(i)) then
            icntix(i) = icnt_snd(j,i)
          end if
        end do
      end do


      call mp_max( icntix(1:8), intra_image_comm )
      WRITE( stdout, fmt="(3X,'Dipole init ')" )
      DO i = 1, 8
        WRITE( stdout, fmt="(3X,'icntix ',I3,' = ',I5)" ) i, icntix(i)
      END DO

      CALL ln_closeup( )

      DEALLOCATE(icnt_snd)
      DEALLOCATE(icnt_rcv)
      ! workaround: sortedig_l2g no longer needed after this routine
      DEALLOCATE(sortedig_l2g)

      RETURN
      END SUBROUTINE berry_setup


      SUBROUTINE berry_closeup( )
        IF( allocated( indi_l    ) )  deallocate(INDI_L   )
        IF( allocated( sour_indi ) )  deallocate(SOUR_INDI)
        IF( allocated( dest_indi ) )  deallocate(DEST_INDI)
        RETURN
      END SUBROUTINE berry_closeup
     

        SUBROUTINE indi_of_ig( mill, indi )

! compute the array "indi" containing the position of
! translated G vectors, given the array of miller ( mill ) indexes of the
! G vectors.
! mill( 1 : 3 )  miller index of a G vectors
! indi( 1 ) = index of G whose miller index are:  mill(1) + 1, mill(2), mill(3)

          IMPLICIT NONE
          INTEGER :: LN_IND
          EXTERNAL LN_IND
!
          INTEGER, INTENT(IN) :: mill(:)
          INTEGER, INTENT(OUT) :: indi(:)
!
          INTEGER :: iri1, iri2, iri3, iricheck
!
          iri1 = mill(1)
          iri2 = mill(2)
          iri3 = mill(3)
          iricheck = iri1**2 + iri2**2 + iri3**2

          if( iricheck == 0 ) then

            !  only positive directions for Gamma point when Gamma symmetry is used

            INDI(1) = LN_IND(1,0,0)
            INDI(2) = 0
            INDI(3) = 0
            INDI(4) = LN_IND(0,1,0)
            INDI(5) = 0
            INDI(6) = 0
            INDI(7) = LN_IND(0,0,1)
            INDI(8) = 0

          ELSE

            !  for gamma symmetry  iri1 >= 0

            INDI(1) = LN_IND( IRI1 + 1, IRI2, IRI3 )

            IF( IRI1 > 0 ) THEN
              INDI(2) = LN_IND( IRI1 - 1, IRI2, IRI3 )
            ELSE
              INDI(2) = -1   !  LN_IND( IRI1 + 1, IRI2, IRI3 )
            ENDIF

            iricheck = iri2**2 + iri3**2
            IF( ( IRI1 < 2 ) .and. ( iricheck /= 0 ) ) THEN
              INDI(3) = LN_IND( 1 - IRI1, -IRI2, -IRI3 )
            ELSE
              INDI(3) = -1
            ENDIF

            INDI(4) = LN_IND(IRI1,IRI2+1,IRI3)
            INDI(5) = LN_IND(IRI1,IRI2-1,IRI3)
            IF( ( IRI1 == 0 ) .AND. ( IRI2 < 2 ) .and. ( iri3 /= 0 ) ) THEN
              INDI(6) = LN_IND( 0, 1-IRI2, -IRI3 )
            ELSE
              INDI(6) = -1
            ENDIF

            INDI(7)=LN_IND(IRI1,IRI2,IRI3+1)
            INDI(8)=LN_IND(IRI1,IRI2,IRI3-1)

          END IF

          RETURN
        END SUBROUTINE indi_of_ig



      END MODULE berry_phase
