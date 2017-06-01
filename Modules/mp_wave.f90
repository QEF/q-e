!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
    MODULE mp_wave

      IMPLICIT NONE
      SAVE

    CONTAINS

      SUBROUTINE mergewf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... This subroutine merges the pieces of a wave functions (pw) splitted across 
! ... processors into a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used)

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      COMPLEX(DP), intent(in) :: PW(:)
      COMPLEX(DP), intent(out) :: PWT(:)
      INTEGER, INTENT(IN) :: mpime     ! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc     ! number of processors
      INTEGER, INTENT(IN) :: root      ! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( igwx, itmp, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF ( mpime == root .AND. igwx > SIZE( pwt ) ) &
        CALL errore(' mergewf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid, IERR )
            CALL MPI_SEND( pw(1), ngwl, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_RECV( pw_ip, ngw_lmax, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ngw_ip, ierr ) 
            DO I = 1, ngw_ip
              PWT(ig_ip(i)) = pw_ip(i)
            END DO
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
            DO I = 1, ngwl
              PWT(ig_l2g(i)) = pw(i)
            END DO
          END IF

        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
        PWT( ig_l2g(i) ) = pw(i)
      END DO

#else

      CALL errore(' MERGEWF ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE mergewf

!=----------------------------------------------------------------------------=!
      
      SUBROUTINE mergekg ( mill, millt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... Same logic as for mergewf, for Miller indices:
!...  mill = distributed input, millt = collected output

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      INTEGER, intent(in) :: mill(:,:)
      INTEGER, intent(out):: millt(:,:)
      INTEGER, INTENT(IN) :: mpime     ! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc     ! number of processors
      INTEGER, INTENT(IN) :: root      ! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      INTEGER, ALLOCATABLE :: mill_ip(:,:)

      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( igwx, itmp, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif
      IF ( mpime == root .AND. igwx > SIZE( millt, 2 ) ) &
        CALL errore(' mergekg',' wrong size for millt ',SIZE(millt,2) )

#if defined __MPI

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid, IERR )
            CALL MPI_SEND( mill,3*ngwl, MPI_INTEGER, ROOT, IP+NPROC, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(mill_ip(3,ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_INTEGER, ngw_ip, ierr ) 
            CALL MPI_RECV( mill_ip,3*ngw_lmax, MPI_INTEGER, (ip-1), IP+NPROC, gid, istatus, IERR )
            DO I = 1,ngw_ip
              millt(:,ig_ip(i)) = mill_ip(:,i)
            END DO
            DEALLOCATE(ig_ip)
            DEALLOCATE(mill_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
            DO I = 1, ngwl
              millt(:,ig_l2g(i)) = mill(:,i)
            END DO
          END IF

        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
         millt(:,ig_l2g(i) ) = mill(:,i)
      END DO

#else

      CALL errore(' mergekg ',' no communication protocol ',0)

#endif

      RETURN
    END SUBROUTINE mergekg

!=----------------------------------------------------------------------------=!

      SUBROUTINE splitwf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... This subroutine splits a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used), across 
! ... processors (pw).

      USE kinds
      USE parallel_include
      IMPLICIT NONE

      COMPLEX(DP), INTENT(OUT) :: PW(:)
      COMPLEX(DP), INTENT(IN) :: PWT(:)
      INTEGER, INTENT(IN) :: mpime, nproc, root
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE(igwx, itmp    , 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF ( mpime == root .AND. igwx > SIZE( pwt ) ) &
        CALL errore(' splitwf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc
! ...   In turn each processor send to root the the indexes of its wavefunction conponents
! ...   Root receive the indexes and send the componens of the wavefunction read from the disk (pwt)
        IF ( (ip-1) /= root ) THEN
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid,IERR)
            CALL MPI_RECV( pw(1), ngwl, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, istatus, IERR )
          END IF
          IF ( mpime == root ) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
            DO i = 1, ngw_ip
              pw_ip(i) = PWT(ig_ip(i))
            END DO
            CALL MPI_SEND( pw_ip, ngw_ip, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, IERR )
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF
        ELSE
          IF ( mpime == root ) THEN
            DO i = 1, ngwl
              pw(i) = PWT(ig_l2g(i)) 
            END DO
          END IF
        END IF
        CALL MPI_BARRIER(gid, IERR)
      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
        pw(i) = pwt( ig_l2g(i) )
      END DO

#else

      CALL errore(' SPLITWF ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE splitwf

      !=----------------------------------------------------------------------------=!

      SUBROUTINE splitkg ( mill, millt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... Same logic as for splitwf, for Miller indices:
!...  mill = distributed output, millt = collected input

      USE kinds
      USE parallel_include
      IMPLICIT NONE

      INTEGER, INTENT(OUT):: mill(:,:)
      INTEGER, INTENT(IN) :: millt(:,:)
      INTEGER, INTENT(IN) :: mpime, nproc, root
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      INTEGER, ALLOCATABLE :: mill_ip(:,:)

      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE(igwx, itmp    , 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF ( mpime == root .AND. igwx > SIZE( millt, 2 ) ) &
        CALL errore(' splitkg ',' wrong size for millt ',SIZE(millt,2) )

#if defined __MPI

      DO ip = 1, nproc
! ...   In turn each processor send to root the the indexes of its wavefunction conponents
! ...   Root receive the indexes and send the componens of the wavefunction read from the disk (pwt)
        IF ( (ip-1) /= root ) THEN
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid,IERR)
            CALL MPI_RECV( mill(1,1),3*ngwl, MPI_INTEGER, ROOT, IP+NPROC, gid, istatus, IERR )
          END IF
          IF ( mpime == root ) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(mill_ip(3,ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
            DO i = 1, ngw_ip
              mill_ip(:,i) = millt(:,ig_ip(i))
            END DO
            CALL MPI_SEND( mill_ip, 3*ngw_ip, MPI_INTEGER, (ip-1), IP+NPROC, gid, IERR )
            DEALLOCATE(ig_ip)
            DEALLOCATE(mill_ip)
          END IF
        ELSE
          IF ( mpime == root ) THEN
            DO i = 1, ngwl
              mill(:,i) = millt(:,ig_l2g(i)) 
            END DO
          END IF
        END IF
        CALL MPI_BARRIER(gid, IERR)
      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
         mill(:,i) = millt(:,ig_l2g(i)) 
      END DO

#else

      CALL errore(' SPLITWF ',' no communication protocol ',0)

#endif

      RETURN
    END SUBROUTINE splitkg

      SUBROUTINE mergeig(igl, igtot, ngl, mpime, nproc, root, comm)

! ... This subroutine merges the pieces of a vector splitted across 
! ... processors into a total vector (igtot) containing al the components
! ... in a pre-defined order (the same as if only one processor is used)

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      INTEGER, intent(in)  :: igl(:)
      INTEGER, intent(out) :: igtot(:)
      INTEGER, INTENT(IN) :: mpime     ! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc     ! number of processors
      INTEGER, INTENT(IN) :: root      ! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ngl

      INTEGER, ALLOCATABLE :: ig_ip(:)

      INTEGER :: ierr, i, ip, ng_ip, ng_lmax, ng_g, gid, igs

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngl, ng_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( ngl, ng_g   , 1, MPI_INTEGER, MPI_SUM, gid, IERR )
      IF( ng_g > SIZE( igtot ) ) THEN
        CALL errore(' mergeig ',' wrong size for igtot ',SIZE(igtot) )
      END IF

      igs = 1

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( igl(1), ngl, MPI_INTEGER, ROOT, IP, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE( ig_ip(ng_lmax) )
            CALL MPI_RECV( ig_ip, ng_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_INTEGER, ng_ip, ierr ) 
            DO i = 1, ng_ip
              igtot( igs + i - 1 ) = ig_ip( i )
            END DO
            DEALLOCATE(ig_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
            ng_ip = ngl
            DO i = 1, ngl
              igtot( igs + i - 1 ) = igl( i )
            END DO
          END IF

        END IF

        IF(mpime == root) THEN
          igs = igs + ng_ip
        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#elif ! defined __MPI

      igtot( 1:ngl ) = igl( 1:ngl )

#else

      CALL errore(' mergeig ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE mergeig

!=----------------------------------------------------------------------------=!

      SUBROUTINE splitig(igl, igtot, ngl, mpime, nproc, root, comm)

! ... This subroutine splits a replicated vector (igtot) stored on the root proc
! ... across processors (igl).

      USE kinds
      USE parallel_include
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: igl(:)
      INTEGER, INTENT(IN)  :: igtot(:)
      INTEGER, INTENT(IN)  :: mpime, nproc, root
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ngl

      INTEGER ierr, i, ng_ip, ip, ng_lmax, ng_g, gid, igs

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

      INTEGER, ALLOCATABLE :: ig_ip(:)

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE(ngl, ng_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE(ngl, ng_g   , 1, MPI_INTEGER, MPI_SUM, gid, IERR )
      IF( ng_g > SIZE( igtot ) ) THEN
        CALL errore(' splitig ',' wrong size for igtot ', SIZE(igtot) )
      END IF

      igs = 1

      DO ip = 1, nproc

! ...   In turn each processor sends to root the indices of its wavefunction components
! ...   Root receives the indices and sends the components of the wavefunction read from the disk (pwt)

        IF ( (ip-1) /= root ) THEN

          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ngl, 1  , MPI_INTEGER, ROOT, IP, gid,IERR)
            CALL MPI_RECV( igl, ngl, MPI_INTEGER, ROOT, IP+NPROC, gid, istatus, IERR )
          END IF

          IF ( mpime == root ) THEN
            ALLOCATE(ig_ip(ng_lmax))
            CALL MPI_RECV( ng_ip, 1, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            DO i = 1, ng_ip
              ig_ip(i) = igtot( igs + i - 1)
            END DO
            CALL MPI_SEND( ig_ip, ng_ip, MPI_INTEGER, (ip-1), IP+NPROC, gid, IERR )
            DEALLOCATE(ig_ip)
          END IF

        ELSE

          IF ( mpime == root ) THEN
            ng_ip = ngl
            DO i = 1, ng_ip
              igl(i) = igtot( igs + i - 1)
            END DO
          END IF

        END IF

        IF( mpime == root ) igs = igs + ng_ip

        CALL MPI_BARRIER(gid, IERR)

      END DO

#elif ! defined __MPI

      igl( 1:ngl ) = igtot( 1:ngl )

#else

      CALL errore(' splitig ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE splitig

!=----------------------------------------------------------------------------=!

   SUBROUTINE pwscatter( c, ctmp, ngw, indi_l, sour_indi, dest_indi, &
      n_indi_rcv, n_indi_snd, icntix, mpime, nproc, group )

      USE kinds
      USE parallel_include

      implicit none

      integer :: indi_l(:)     !  list of G-vec index to be exchanged
      integer :: sour_indi(:)  !  the list of source processors
      integer :: dest_indi(:)  !  the list of destination processors
      integer :: n_indi_rcv    !   number of G-vectors to be received
      integer :: n_indi_snd    !   number of G-vectors to be sent
      integer :: icntix        !   total number of G-vec to be exchanged
      INTEGER, INTENT(IN) :: nproc, mpime, group

      COMPLEX(DP) :: c(:)
      COMPLEX(DP) :: ctmp(:)
      integer  ::  ngw

      integer :: ig, icsize
      INTEGER :: me, idest, isour, ierr

      COMPLEX(DP), ALLOCATABLE :: my_buffer( : )
      COMPLEX(DP), ALLOCATABLE :: mp_snd_buffer( : )
      COMPLEX(DP), ALLOCATABLE :: mp_rcv_buffer( : )
      INTEGER, ALLOCATABLE :: ibuf(:)

      !
      ! ... SUBROUTINE BODY
      !

      me = mpime + 1

      if( icntix .lt. 1 ) then
        icsize = 1
      else
        icsize = icntix
      endif

      ALLOCATE( mp_snd_buffer( icsize * nproc ) )
      ALLOCATE( mp_rcv_buffer( icsize * nproc ) )
      ALLOCATE( my_buffer( ngw ) )
      ALLOCATE( ibuf( nproc ) )
      ctmp = ( 0.0_DP, 0.0_DP )

      ! WRITE( stdout,*) 'D: ', nproc, mpime, group

      ibuf = 0
      DO IG = 1, n_indi_snd
        idest = dest_indi(ig)
        ibuf(idest) = ibuf(idest) + 1;
        if(idest .ne. me) then
          mp_snd_buffer( ibuf(idest) + (idest-1)*icsize ) = C( indi_l( ig ) )
        else
          my_buffer(ibuf(idest)) = C(indi_l(ig))
        end if
      end do

#if defined __MPI
      call MPI_ALLTOALL( mp_snd_buffer(1), icsize, MPI_DOUBLE_COMPLEX, &
                         mp_rcv_buffer(1), icsize, MPI_DOUBLE_COMPLEX, &
                         group, ierr)
#else

      CALL errore(' pwscatter ',' no communication protocol ',0)

#endif

      ibuf = 0
      DO IG = 1, n_indi_rcv
        isour = sour_indi(ig)
        if(isour.gt.0 .and. isour.ne.me) then
          ibuf(isour) = ibuf(isour) + 1
          CTMP(ig) = mp_rcv_buffer(ibuf(isour) + (isour-1)*icsize)
        else if(isour.gt.0) then
          ibuf(isour) = ibuf(isour) + 1
          CTMP(ig) = my_buffer(ibuf(isour))
        else
          CTMP(ig) = (0.0_DP,0.0_DP)
        end if
      end do

      DEALLOCATE( mp_snd_buffer )
      DEALLOCATE( mp_rcv_buffer )
      DEALLOCATE( my_buffer )
      DEALLOCATE( ibuf )

      RETURN
    END SUBROUTINE pwscatter




!=----------------------------------------------------------------------------=!

SUBROUTINE redistwf( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !  Redistribute wave function.
   !  c_dist_pw are the wave functions with plane waves distributed over processors 
   !  c_dist_st are the wave functions with electronic states distributed over processors 
   !
   USE kinds
   USE parallel_include

   implicit none

   COMPLEX(DP) :: c_dist_pw(:,:)
   COMPLEX(DP) :: c_dist_st(:,:)
   INTEGER, INTENT(IN) :: npw_p(:)  !  the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)  !  the number of states on each processor
   INTEGER, INTENT(IN) :: comm      !  group communicator
   INTEGER, INTENT(IN) :: idir      !  direction of the redistribution 
                                    !  idir > 0  c_dist_pw --> c_dist_st
                                    !  idir < 0  c_dist_pw <-- c_dist_st

   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww, ii
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   COMPLEX(DP), ALLOCATABLE :: ctmp( : )

#if defined(__MPI)
   CALL mpi_comm_rank( comm, mpime, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_rank ', ierr )
   CALL mpi_comm_size( comm, nproc, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_size ', ierr )

   ALLOCATE( rdispls( nproc ), recvcount( nproc ), sendcount( nproc ), sdispls( nproc ) )

   npw_t = 0
   nst_t = 0
   DO proc=1,nproc
      sendcount(proc) = npw_p(mpime+1) * nst_p(proc)
      recvcount(proc) = npw_p(proc) * nst_p(mpime+1)
      npw_t = npw_t + npw_p(proc)
      nst_t = nst_t + nst_p(proc)
   END DO
   sdispls(1)=0
   rdispls(1)=0
   DO proc=2,nproc
      sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
      rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
   END DO

   ALLOCATE( ctmp( npw_t * nst_p( mpime + 1 ) ) )

   IF( idir > 0 ) THEN
      !
      ! ... Step 1. Communicate to all Procs so that each proc has all
      ! ... G-vectors and some states instead of all states and some
      ! ... G-vectors. This information is stored in the 1-d array ctmp.
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )
      !
      CALL MPI_ALLTOALLV( c_dist_pw, sendcount, sdispls, MPI_DOUBLE_COMPLEX,             &
           &             ctmp, recvcount, rdispls, MPI_DOUBLE_COMPLEX, comm, ierr)
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )
      !
      !   Step 2. Convert the 1-d array ctmp into a 2-d array consistent with the
      !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
      !   and nstat states instead of all nbsp states
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1)
            ii = (i-1) * npw_p(proc) 
            DO j = 1, npw_p(proc)
               c_dist_st( j + ngpww, i ) = ctmp( rdispls(proc) + j + ii )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO

   ELSE
      !
      !   Step 4. Convert the 2-d array c_dist_st into 1-d array
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1) 
            ii = (i-1) * npw_p(proc)
            DO j = 1, npw_p(proc)
               ctmp( rdispls(proc) + j + ii ) = c_dist_st( j + ngpww, i )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO
      !        
      !   Step 5. Redistribute among processors. The result is stored in 2-d
      !   array c_dist_pw consistent with the notation c(ngw,nbsp)
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )

      CALL MPI_ALLTOALLV( ctmp, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
          &               c_dist_pw, sendcount , sdispls, MPI_DOUBLE_COMPLEX, comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )


   END IF

   DEALLOCATE( ctmp )
   DEALLOCATE( rdispls, recvcount, sendcount, sdispls )
#endif
   RETURN
END SUBROUTINE redistwf

!=----------------------------------------------------------------------------=!

SUBROUTINE redistwfr( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !  Redistribute wave function.
   !  c_dist_pw are the wave functions with plane waves distributed over processors 
   !  c_dist_st are the wave functions with electronic states distributed over processors 
   !
   USE kinds
   USE parallel_include

   implicit none

   REAL(DP) :: c_dist_pw(:,:)
   REAL(DP) :: c_dist_st(:,:)
   INTEGER, INTENT(IN) :: npw_p(:)  !  the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)  !  the number of states on each processor
   INTEGER, INTENT(IN) :: comm      !  group communicator
   INTEGER, INTENT(IN) :: idir      !  direction of the redistribution 
                                    !  idir > 0  c_dist_pw --> c_dist_st
                                    !  idir < 0  c_dist_pw <-- c_dist_st

   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   REAL(DP), ALLOCATABLE :: ctmp( : )

#if defined(__MPI)
   CALL mpi_comm_rank( comm, mpime, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_rank ', ierr )
   CALL mpi_comm_size( comm, nproc, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_size ', ierr )

   ALLOCATE( rdispls( nproc ), recvcount( nproc ), sendcount( nproc ), sdispls( nproc ) )

   npw_t = 0
   nst_t = 0
   DO proc=1,nproc
      sendcount(proc) = npw_p(mpime+1) * nst_p(proc)
      recvcount(proc) = npw_p(proc) * nst_p(mpime+1)
      npw_t = npw_t + npw_p(proc)
      nst_t = nst_t + nst_p(proc)
   END DO
   sdispls(1)=0
   rdispls(1)=0
   DO proc=2,nproc
      sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
      rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
   END DO

   ALLOCATE( ctmp( npw_t * nst_p( mpime + 1 ) ) )

   IF( idir > 0 ) THEN
      !
      ! ... Step 1. Communicate to all Procs so that each proc has all
      ! ... G-vectors and some states instead of all states and some
      ! ... G-vectors. This information is stored in the 1-d array ctmp.
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )
      !
      CALL MPI_ALLTOALLV( c_dist_pw, sendcount, sdispls, MPI_DOUBLE_PRECISION,             &
           &             ctmp, recvcount, rdispls, MPI_DOUBLE_PRECISION, comm, ierr)
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )
      !
      !   Step 2. Convert the 1-d array ctmp into a 2-d array consistent with the
      !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
      !   and nstat states instead of all nbsp states
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1)
            DO j = 1, npw_p(proc)
               c_dist_st( j + ngpww, i ) = ctmp( rdispls(proc) + j + (i-1) * npw_p(proc) )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO

   ELSE
      !
      !   Step 4. Convert the 2-d array c_dist_st into 1-d array
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1) 
            DO j = 1, npw_p(proc)
               ctmp( rdispls(proc) + j + (i-1) * npw_p(proc) ) = c_dist_st( j + ngpww, i )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO
      !        
      !   Step 5. Redistribute among processors. The result is stored in 2-d
      !   array c_dist_pw consistent with the notation c(ngw,nbsp)
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )

      CALL MPI_ALLTOALLV( ctmp, recvcount, rdispls, MPI_DOUBLE_PRECISION,          &
          &               c_dist_pw, sendcount , sdispls, MPI_DOUBLE_PRECISION, comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )


   END IF

   DEALLOCATE( ctmp )
   DEALLOCATE( rdispls, recvcount, sendcount, sdispls )
#endif
   RETURN
END SUBROUTINE redistwfr

!=----------------------------------------------------------------------------=!

    END MODULE mp_wave

