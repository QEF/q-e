!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 MODULE mp_wave_parallel

      IMPLICIT NONE
      SAVE

    CONTAINS

      SUBROUTINE mergewfp ( npw,pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

! ... This subroutine merges the pieces of a wave functions (pw) splitted across 
! ... processors into a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used)

      USE kinds
      USE parallel_include
      USE io_global, ONLY :stdout
     

      IMPLICIT NONE

      INTEGER, INTENT(in) :: npw,ngwl
      COMPLEX(DP), intent(in) :: PW(npw,nproc)
      COMPLEX(DP), intent(out) :: PWT(:)
      INTEGER, INTENT(IN) :: mpime     ! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc     ! number of processors
      INTEGER, INTENT(IN) :: root      ! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      
     

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)


      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid, req

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

      INTEGER :: iorig, idest
     

     

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

      IF( igwx > SIZE( pwt ) ) &
        CALL errore(' mergewf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI
      ALLOCATE(ig_ip(ngw_lmax))
      ALLOCATE(pw_ip(ngw_lmax))
       
      do ip = 0, nproc-1

         if( ip/=0) then
            
! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
            idest=mpime+ip
            if(idest>nproc-1)idest=idest-nproc
         
            iorig=mpime-ip
            if(iorig<0)iorig=iorig+nproc


            CALL MPI_ISEND( ig_l2g, ngwl, MPI_INTEGER, idest, IP, gid, req,IERR )


            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, iorig, IP, gid, istatus, IERR )
    
            CALL MPI_WAIT(req,istatus,ierr)
       
            CALL MPI_ISEND( pw(1,idest+1), ngwl, MPI_DOUBLE_COMPLEX, idest, IP, gid, req,IERR )
        

            CALL MPI_RECV( pw_ip, ngw_lmax, MPI_DOUBLE_COMPLEX, iorig, IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ngw_ip, ierr ) 

            CALL MPI_WAIT(req,istatus,ierr)

            
            DO I = 1, ngw_ip
               PWT(ig_ip(i)) = pw_ip(i)
            END DO
         

            
        ELSE

         
           DO I = 1, ngwl
              PWT(ig_l2g(i)) = pw(i,mpime+1)
           END DO
        END IF

    

        CALL MPI_BARRIER( gid, IERR )
            
     END DO

     DEALLOCATE(ig_ip)
     DEALLOCATE(pw_ip)
    

        
#elif ! defined __MPI

     DO I = 1, ngwl
        PWT( ig_l2g(i) ) = pw(i,1)
     END DO

#else

      CALL errore(' MERGEWF ',' no communication protocol ',0)

#endif

      RETURN
    END SUBROUTINE mergewfp

    SUBROUTINE splitwfp (npw, pw, pwt, ngwl, ig_l2g, mpime, nproc,root, comm )

! ... This subroutine splits a total wave function (pwt) containing al the components
! ... in a pre-defined order (the same as if only one processor is used), across 
! ... processors (pw).

      USE kinds
      USE parallel_include
      USE io_global, ONLY : stdout 
      IMPLICIT NONE

      INTEGER, INTENT(in) :: npw,nproc
      COMPLEX(DP), INTENT(OUT) :: PW(npw,nproc)
      COMPLEX(DP), INTENT(IN) :: PWT(:)
      INTEGER, INTENT(IN) :: mpime,  root
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp,len, req

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

      INTEGER :: iorig, idest
  
     
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

      IF( igwx > SIZE( pwt ) ) &
        CALL errore(' splitwf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI
      ALLOCATE(ig_ip(ngw_lmax))
      ALLOCATE(pw_ip(ngw_lmax))
     


      DO ip = 0, nproc-1

      

         idest=mpime+ip
         if(idest>nproc-1)idest=idest-nproc
         
         iorig=mpime-ip
         if(iorig<0)iorig=iorig+nproc


         if(ip/=0) then
          
            
            CALL MPI_ISEND( ig_l2g, ngwl, MPI_INTEGER, iorig, IP, gid,req,IERR)
       
           
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, idest, IP, gid, istatus, IERR )

            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
            DO i = 1, ngw_ip
               pw_ip(i) = PWT(ig_ip(i))
            END DO
           
            CALL MPI_WAIT(req,istatus,ierr)

            CALL MPI_ISEND( pw_ip, ngw_ip, MPI_DOUBLE_COMPLEX, idest, IP, gid,req, IERR )
            CALL MPI_RECV( pw(1,iorig+1), ngwl, MPI_DOUBLE_COMPLEX, iorig, IP, gid, istatus, IERR )
            !CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
                     
            CALL MPI_WAIT(req,istatus,ierr)
         ELSE
          

            DO i = 1, ngwl
              pw(i,mpime+1) = PWT(ig_l2g(i)) 
            END DO
            
         END IF
         CALL MPI_BARRIER(gid, IERR)
        

      END DO
      DEALLOCATE(ig_ip)
      DEALLOCATE(pw_ip)
    
       

#elif ! defined __MPI

      DO I = 1, ngwl
        pw(i,1) = pwt( ig_l2g(i) )
      END DO

#else

      CALL errore(' SPLITWF ',' no communication protocol ',0)

#endif

      RETURN
    END SUBROUTINE splitwfp


    


  END MODULE mp_wave_parallel

  SUBROUTINE reorderwfp (nbands,npw1, npw2,pw1,pw2, ngwl1,ngwl2, ig_l2g1,ig_l2g2,n_g,mpime, nproc,root, comm )
    
      USE kinds
      USE parallel_include
      USE io_global, ONLY : stdout
      USE mp_wave_parallel

      IMPLICIT NONE

      INTEGER, INTENT(in) :: npw1,npw2,nbands
      COMPLEX(DP), INTENT(OUT) :: pw1(npw1,nbands),pw2(npw2,nbands)
      INTEGER, INTENT(IN) :: mpime,  root, nproc
      INTEGER, INTENT(IN) :: comm    ! communicator
      INTEGER, INTENT(IN) :: ig_l2g1(ngwl1),ig_l2g2(ngwl2)
      INTEGER, INTENT(IN) :: ngwl1,ngwl2
      INTEGER, INTENT(in) :: n_g!global maximum number of G vectors for both grids
      

      COMPLEX(kind=DP), ALLOCATABLE :: cbuf1(:,:),cbuf2(:,:), pwt(:)
      INTEGER :: ii, ilast


      allocate(cbuf1(npw1,nproc),cbuf2(npw2,nproc))
      allocate(pwt(n_g))
      cbuf1(:,:)=(0.d0,0.d0)
      cbuf2(:,:)=(0.d0,0.d0)
      !loop on block of states
      
      do ii=1,nbands,nproc
         ilast=min(nbands,ii+nproc-1)
         cbuf1(1:npw1,1:(ilast-ii+1))=pw1(1:npw1,ii:ilast)
         call mergewfp ( npw1,cbuf1, pwt, ngwl1, ig_l2g1, mpime, nproc, root, comm )
         call splitwfp (npw2, cbuf2, pwt, ngwl2, ig_l2g2, mpime, nproc,root, comm )
         pw2(1:npw2,ii:ilast)=cbuf2(1:npw2,1:(ilast-ii+1))
      enddo
      deallocate(cbuf1,cbuf2)
      deallocate(pwt)
      return
    END SUBROUTINE reorderwfp

    SUBROUTINE reorderwfp_col (nbands,npw1, npw2,pw1,pw2, ngwl1,ngwl2, ig_l2g1,ig_l2g2,n_g,mpime, nproc, comm )
!routine using collective mpi calls
      USE kinds
      USE parallel_include
      USE io_global, ONLY : stdout
      USE mp_wave_parallel

      IMPLICIT NONE

      INTEGER, INTENT(in) :: npw1,npw2,nbands
      COMPLEX(kind=DP) :: pw1(npw1,nbands),pw2(npw2,nbands)
      INTEGER, INTENT(IN) :: mpime, nproc
      INTEGER, INTENT(IN) :: comm    ! communicator                           
      INTEGER, INTENT(IN) :: ig_l2g1(ngwl1),ig_l2g2(ngwl2)
      INTEGER, INTENT(IN) :: ngwl1,ngwl2
      INTEGER, INTENT(in) :: n_g!global maximum number of G vectors for both grids

      INTEGER :: ngwl1_max,ngwl2_max,npw1_max,npw2_max
      INTEGER :: gid,ierr
      INTEGER, ALLOCATABLE :: npw1_loc(:),npw2_loc(:)
      INTEGER, ALLOCATABLE :: ig_l2g1_tot(:,:),ig_l2g2_tot(:,:), itmp(:)

      INTEGER :: ii,ip,ilast,iband
      COMPLEX(kind=DP), ALLOCATABLE :: pw1_tot(:,:),pw2_tot(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: pw1_tmp(:),pw2_tmp(:), pw_global(:)

      gid=comm

#if defined __MPI
      allocate(npw1_loc(nproc),npw2_loc(nproc))
!all procs gather correspondance arrays
      CALL MPI_ALLREDUCE( ngwl1, ngwl1_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( ngwl2, ngwl2_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( npw1, npw1_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( npw2, npw2_max, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLGATHER (npw1,1,MPI_INTEGER,npw1_loc,1,MPI_INTEGER,gid,IERR)      
      CALL MPI_ALLGATHER (npw2,1,MPI_INTEGER,npw2_loc,1,MPI_INTEGER,gid,IERR)

      allocate(ig_l2g1_tot(ngwl1_max,nproc),ig_l2g2_tot(ngwl2_max,nproc))
      allocate(itmp(ngwl1_max))
      itmp(1:ngwl1)=ig_l2g1(1:ngwl1)
      CALL MPI_ALLGATHER (itmp,ngwl1_max,MPI_INTEGER,ig_l2g1_tot,ngwl1_max,MPI_INTEGER,gid,IERR)
      deallocate(itmp)
      allocate(itmp(ngwl2_max))
      itmp(1:ngwl2)=ig_l2g2(1:ngwl2)
      CALL MPI_ALLGATHER (itmp,ngwl2_max,MPI_INTEGER,ig_l2g2_tot,ngwl2_max,MPI_INTEGER,gid,IERR)
      deallocate(itmp)

      allocate(pw1_tot(npw1_max,nproc),pw2_tot(npw2_max,nproc))
      allocate(pw1_tmp(npw1_max),pw2_tmp(npw2_max))
      allocate(pw_global(n_g))

      do ii=1,nbands,nproc
         ilast=min(nbands,ii+nproc-1)
         do iband=ii,ilast
            ip=mod(iband,nproc)!ip starts from 1 to nproc-1
            pw1_tmp(1:npw1)=pw1(1:npw1,iband)
            CALL MPI_GATHER (pw1_tmp,npw1_max,MPI_DOUBLE_COMPLEX,pw1_tot,npw1_max,MPI_DOUBLE_COMPLEX,ip,gid,ierr)
         enddo
         pw_global=0.d0
         do ip=1,nproc
            pw_global(ig_l2g1_tot(1:npw1_loc(ip),ip))=pw1_tot(1:npw1_loc(ip),ip)
         enddo
         do ip=1,nproc
            pw2_tot(1:npw2_loc(ip),ip)=pw_global(ig_l2g2_tot(1:npw2_loc(ip),ip))
         enddo
         do iband=ii,ilast
            ip=mod(iband,nproc)
            CALL MPI_SCATTER (pw2_tot,npw2_max,MPI_DOUBLE_COMPLEX,pw2_tmp,npw2_max ,MPI_DOUBLE_COMPLEX,ip,gid,ierr)
            pw2(1:npw2,iband)=pw2_tmp(1:npw2)
         enddo
      enddo

      deallocate(npw1_loc,npw2_loc)
      deallocate(ig_l2g1_tot,ig_l2g2_tot)
      deallocate(pw1_tot,pw2_tot)
      deallocate(pw1_tmp,pw2_tmp)
      deallocate(pw_global)
#endif
      return
    END SUBROUTINE reorderwfp_col
