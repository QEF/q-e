!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


  MODULE  fft_gw
!this modules contains the structures and subroutine
!which permits ffts: the structures polaw are read from
!disk in a series of rows a then ffts are performed
!SEE INTERNAL NOTES
   USE kinds, only : DP

   TYPE fft_data
!this structures described a series of row from polaw
!at all times/frequencies in imaginary time
!to construct from polaw the convention for the labels is:
!label=i , i=-n,+n
     INTEGER :: label!label to read/write from/to disk
     LOGICAL :: ontime!if .true. data are  on imaginary time
     INTEGER :: numpw!number of coloumns same as in polaw
     INTEGER :: numrows!number of rows
     INTEGER :: firstrow!first row included
     INTEGER :: lastrow!last row included
     REAL(kind=DP) :: period!max time tau data are from O to tau
     INTEGER :: n !number of campions on time/frequency T
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: fd!data in format (numpw,numrows,2*n+1)
     COMPLEX :: factor!used for real matrices
!we take advantage of the symmetry t ==> -t
! the format will be for GL grid (numpw,numrows,n+1)
  END TYPE fft_data


 CONTAINS

 SUBROUTINE free_memory_fft_data(fftd)
!this subroutine  deallocates the fft descriptor
   implicit none
   TYPE(fft_data) :: fftd
   if(associated(fftd%fd)) deallocate(fftd%fd)
   nullify(fftd%fd)
   return
 END SUBROUTINE
  


 SUBROUTINE read_fft_data(label,fftd,debug)
!this subroutine reads the fft descriptor from file
!we take care of the t ==> -t symmetry
      USE io_files,  ONLY : prefix, tmp_dir
   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   TYPE(fft_data) :: fftd
   INTEGER :: label !label for the corresponding file
   LOGICAL :: debug !if true formatted files
   
   INTEGER :: iunf, iw,it, jw
   CHARACTER(5) :: nfile

   CALL free_memory_fft_data(fftd)

!open file
    if(label >= 0 ) then
      write(nfile,'(5i1)') &
        & label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
      iunf = find_free_unit()
      if(.not.debug) then
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.'// nfile, status='old',form='unformatted')
      else
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.'// nfile, status='old',form='formatted')
      endif
    else
      write(nfile,'(5i1)') &
        & -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
      iunf = find_free_unit()
      if(.not.debug) then
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.-'// nfile, status='old',form='unformatted')
      else
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.-'// nfile, status='old',form='formatted')
      endif

    endif
    if(.not.debug) then
      read(iunf)  fftd%label
      read(iunf)  fftd%ontime
      read(iunf)  fftd%numpw
      read(iunf)  fftd%numrows
      read(iunf)  fftd%firstrow
      read(iunf)  fftd%lastrow
      read(iunf)  fftd%period
      read(iunf)  fftd%n   
    else
      read(iunf,*)  fftd%label
      read(iunf,*)  fftd%ontime
      read(iunf,*)  fftd%numpw
      read(iunf,*)  fftd%numrows
      read(iunf,*)  fftd%firstrow
      read(iunf,*)  fftd%lastrow
      read(iunf,*)  fftd%period
      read(iunf,*)  fftd%n
    endif

!    allocate(fftd%fd(fftd%numpw,fftd%numrows,2*fftd%n+2))
    allocate(fftd%fd(fftd%numpw,fftd%numrows,fftd%n+1))
    if(.not.debug) then
!      do it=1,2*fftd%n+2
       do it=1,fftd%n+1
        do iw=1,fftd%numrows
            read(iunf) fftd%fd(1:fftd%numpw,iw,it)
        enddo
      enddo
    else
!      do it=1,2*fftd%n+2
       do it=1,fftd%n+1
       do iw=1,fftd%numrows
          do jw=1,fftd%numpw
            read(iunf,*) fftd%fd(jw,iw,it)
          enddo
        enddo
      enddo
    endif

    close(iunf)
   
    return
  END SUBROUTINE
   
 SUBROUTINE write_fft_data(fftd,debug)
!this subroutine writes the fft descriptor on file
!we take care of the t ==> -t symmetry
   USE io_files,  ONLY : prefix, tmp_dir
   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   TYPE(fft_data) :: fftd
   LOGICAL :: debug!if true formatted output

   INTEGER :: iunf, iw,it,jw
   CHARACTER(5) :: nfile


!open file
    if(fftd%label >= 0 ) then
      write(nfile,'(5i1)') &
        & fftd%label/10000,mod(fftd%label,10000)/1000,mod(fftd%label,1000)/100,mod(fftd%label,100)/10,mod(fftd%label,10)
      iunf = find_free_unit()
      if(.not.debug) then
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.'// nfile, status='unknown',form='formatted')
      endif
    else
      write(nfile,'(5i1)') &
        & -fftd%label/10000,mod(-fftd%label,10000)/1000,mod(-fftd%label,1000)/100,mod(-fftd%label,100)/10,mod(-fftd%label,10)
      iunf = find_free_unit()
      if(.not.debug) then
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.-'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iunf, file=trim(tmp_dir)//trim(prefix)//'-'//'fftdata.-'// nfile, status='unknown',form='formatted')
      endif
    endif

    if(.not.debug) then
      write(iunf)  fftd%label
      write(iunf)  fftd%ontime
      write(iunf)  fftd%numpw
      write(iunf)  fftd%numrows
      write(iunf)  fftd%firstrow
      write(iunf)  fftd%lastrow
      write(iunf)  fftd%period
      write(iunf)  fftd%n

!      do it=1,2*fftd%n+2
      do it=1,fftd%n+1
        do iw=1,fftd%numrows
            write(iunf) fftd%fd(1:fftd%numpw,iw,it)
        enddo
      enddo
    else
      write(iunf,*)  fftd%label
      write(iunf,*)  fftd%ontime
      write(iunf,*)  fftd%numpw
      write(iunf,*)  fftd%numrows
      write(iunf,*)  fftd%firstrow
      write(iunf,*)  fftd%lastrow
      write(iunf,*)  fftd%period
      write(iunf,*)  fftd%n

!      do it=1,2*fftd%n+2
      do it=1,fftd%n+1
        do iw=1,fftd%numrows
          do jw=1,fftd%numpw
            write(iunf,*) fftd%fd(jw,iw,it)
          enddo
        enddo
      enddo
    endif

    close(iunf)

    return
  END SUBROUTINE

  SUBROUTINE create_fft_data(tf,firstr,lastr,period,n,label,fftd,debug)
! this subroutine creates the descriptor for the fftw reading
!data from polaw on disk
!data is put on appropiate order for FFT
!total period=2*T+T/n
    USE  polarization,        ONLY : polaw,read_polaw, free_memory_polaw,&
                                 & read_polaw_range  
    USE  io_global,           ONLY : stdout
    USE  constants,           ONLY : eps8, pi
    USE  times_gw,            ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
    INTEGER :: firstr !first row to be read (included)
    INTEGER :: lastr  !last row to be read  (included)
    REAL(kind=DP) :: period!period tau (data from -T to T)
    INTEGER :: n ! number of steps on T
    TYPE(fft_data) :: fftd!structure to be initialized
    LOGICAL :: debug!if true formatted files

    TYPE(polaw) :: pw
    INTEGER :: label, il, iw, ipos
    REAL(kind=DP) :: tfcheck, totalperiod

    LOGICAL, PARAMETER ::  direct_access = .true.

!first dealloacate and set


    write(stdout,*) 'VALUE TF', tf%l_fft_timefreq!ATTENZIONE

    totalperiod=2.d0*period+2.d0*period/real(n)

    CALL free_memory_fft_data(fftd)
    fftd%label=label
    fftd%period=period
    fftd%firstrow=firstr
    fftd%lastrow=lastr
    fftd%numrows=abs(lastr-firstr)+1
    fftd%n=n

!read the -n polaw
    if(.not.direct_access) then
!       CALL read_polaw(-n,pw,debug)
       CALL read_polaw(n,pw,debug,.false.)
    else
!       CALL read_polaw_range(-n,pw,debug,firstr,firstr+fftd%numrows-1, .true.)
       CALL read_polaw_range(n,pw,debug,firstr,firstr+fftd%numrows-1, .true. )
    endif
    fftd%ontime=pw%ontime
    fftd%numpw=pw%numpw
!we take care of the t ==> -t symmetry
!    allocate(fftd%fd(fftd%numpw,fftd%numrows,2*fftd%n+2))
    allocate(fftd%fd(fftd%numpw,fftd%numrows,fftd%n+1))
    fftd%fd(:,:,:)=(0.d0,0.d0)
!test time frequency
!does not stop any more to permit prallel execution
    if(tf%l_fft_timefreq) then
       if(fftd%ontime) then!check imaginary time
          tfcheck=-period
          if(abs(pw%time-tfcheck) >= eps8) then
             write(stdout,*) 'routine create_fft_data: times do not correspond  ',n
             !stop
          endif
       else !check imaginary frequency
          tfcheck=(2.d0*pi/totalperiod)*real(-n)
          if(abs(pw%time-tfcheck) >= eps8) then
             write(stdout,*) 'routine create_fft_data: frequencies do not correspond  ',n
             !stop
          endif
       endif
    endif
!put in data at the right position
!the position at 1 is zero for definition
!we take advantage of the t ==> -t symmetry
!    do iw=1,fftd%numrows 
!        fftd%fd(1:fftd%numpw,iw,2) = pw%pw(1:fftd%numpw,firstr+iw-1)   
!    enddo

!read in the other times/frequencies

!we take advantage of the t ==> -t symmetry
!    do il=-n+1,n!loop on time/frequency
    do il=0,n
       CALL read_polaw(il,pw,debug,.false.)

!consistency test

!in case of parallel the check on ontime must not be done
     ! if(pw%ontime .NEQV. fftd%ontime .OR. pw%numpw /= fftd%numpw ) then
       if( pw%numpw /= fftd%numpw ) then
         write(stdout,*) 'routine create_fft_data: consistency failed'
         write(stdout,*) 'il', il
         write(stdout,*) 'ontime',pw%ontime,fftd%ontime, pw%numpw
         stop
      endif

!test time frequency
      if(tf%l_fft_timefreq) then
         if(fftd%ontime) then!check imaginary time
            tfcheck=period/real(n)*real(il)
            if(abs(pw%time-tfcheck) >= eps8) then
               write(stdout,*) 'routine create_fft_data: times do not correspond  ',n
               !stop
            endif
         else !check imaginary frequency
            tfcheck=(2.d0*pi/totalperiod)*real(il)
            if(abs(pw%time-tfcheck) >= eps8) then
               write(stdout,*) 'routine create_fft_data: frequencies do not correspond  ',n
               !stop
            endif
         endif
      endif
!put in data at the right position
!we take care of the t ==> -t symmetry
!       ipos=il+n+2
      ipos=il+1
       do iw=1,fftd%numrows
         fftd%fd(1:fftd%numpw,iw,ipos) = pw%pw(1:fftd%numpw,firstr+iw-1)
       enddo

    enddo
   

   CALL free_memory_polaw(pw)

  END SUBROUTINE



  SUBROUTINE create_fft_data2( tf, firstr, lastr, period, n, fftd, debug )
! this subroutine creates the descriptor for the fftw reading
!data from polaw on disk
!data is put on appropiate order for FFT
!total period=2*T+T/n
    USE polarization,        ONLY : polaw,read_polaw, free_memory_polaw,&
                                 & read_polaw_range  
    USE io_global,           ONLY : stdout
    USE constants,           ONLY : eps8, pi
    USE times_gw,            ONLY : times_freqs
    USE mp_world,            ONLY : nproc, mpime,world_comm! group
    USE parallel_include

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
    INTEGER :: firstr !first row to be read (included)
    INTEGER :: lastr  !last row to be read  (included)
    REAL(kind=DP) :: period!period tau (data from -T to T)
    INTEGER :: n ! number of steps on T
    TYPE(fft_data) :: fftd!structure to be initialized
    LOGICAL :: debug!if true formatted files

    TYPE(polaw) :: pw
    INTEGER :: il, iw, ipos, numrows_read, nblk_siz, k, ierr
    INTEGER :: nbegin, nend, nbegin_ip, nend_ip, ip
    REAL(kind=DP) :: tfcheck, totalperiod
    COMPLEX(kind=DP), ALLOCATABLE :: rcvbuf( : ), sndbuf( : )

    LOGICAL, PARAMETER ::  direct_access = .true.

!first dealloacate and set


    write(stdout,*) 'VALUE TF', tf%l_fft_timefreq!ATTENZIONE

    totalperiod=2.d0*period+2.d0*period/real(n)

    CALL free_memory_fft_data(fftd)
    fftd%label   =0
    fftd%period  =period
    fftd%firstrow=firstr
    fftd%lastrow =lastr

    numrows_read = lastr - firstr + 1 

    fftd%numrows = numrows_read / nproc   
    if( MOD( numrows_read, nproc ) /= 0 ) fftd%numrows = fftd%numrows + 1

    fftd%n = n

    allocate( fftd%fd( fftd%numpw, fftd%numrows, fftd%n + 1 ) )

    fftd%fd(:,:,:) = (0.d0,0.d0)

    nblk_siz = (n + 1) / nproc
    if( MOD( (n + 1), nproc ) /= 0 ) nblk_siz = nblk_siz + 1

    nbegin = 0 + mpime * nblk_siz 
    nend   = nbegin + nblk_siz - 1

    ALLOCATE( sndbuf( fftd%numpw * fftd%numrows * nproc ) )
    ALLOCATE( rcvbuf( fftd%numpw * fftd%numrows * nproc ) )

    do il = nbegin, nend

      if( il <= n ) then

         CALL read_polaw_range( il, pw, debug, firstr, lastr, .false. )

         do iw = 1, numrows_read
            do k = 1, fftd%numpw
               sndbuf( k + fftd%numpw * ( iw - 1 ) ) = pw%pw( k, iw ) ! pw%pw( k, firstr + iw - 1 )
            end do
         end do
         do iw = numrows_read + 1, fftd%numrows * nproc
            do k = 1, fftd%numpw
               sndbuf( k + fftd%numpw * ( iw - 1 ) ) = 0.0d0
            end do
         end do

      else

         sndbuf = 0.0d0

      end if

#if defined(__MPI)
      CALL MPI_ALLTOALL( sndbuf, fftd%numrows * fftd%numpw, MPI_DOUBLE_COMPLEX,  &
                         rcvbuf, fftd%numrows * fftd%numpw, MPI_DOUBLE_COMPLEX, world_comm, ierr )
#else
      rcvbuf(:)=sndbuf(:)
#endif

      do ip = 0, nproc - 1

         nbegin_ip = 0 + ip * nblk_siz 

         ipos = il - nbegin + nbegin_ip 

         if( ipos <= n ) then
            do iw = 1, fftd%numrows
               do k = 1, fftd%numpw
                  fftd%fd( k, iw, ipos + 1 ) = rcvbuf( k + fftd%numpw * (iw-1) + fftd%numrows * fftd%numpw * ip )
               end do
            end do
         end if

      enddo

    enddo
   
    DEALLOCATE( rcvbuf )
    DEALLOCATE( sndbuf )

    CALL free_memory_polaw( pw )

  END SUBROUTINE




  SUBROUTINE save_fft_data(tf, fftd,debug)
! this subroutine writes the descriptor for the fftw 
! on the  polaw on disk
!data is put on appropiate order for FFT
!total period=2*T+T/n
!if we are ON TIME: THE ORDER IS REVERSED
    USE  polarization,        ONLY : polaw,read_polaw,write_polaw, free_memory_polaw, &
                                 &read_polaw_range,write_polaw_range
    USE  io_global,           ONLY : stdout
    USE  constants,           ONLY : eps8, pi
    USE  mp,                  ONLY : mp_barrier
    USE  mp_world,            ONLY : world_comm, mpime, nproc
    USE  times_gw,            ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
    TYPE(fft_data) :: fftd!structure to be written
    LOGICAL :: debug!if true formatted files

    TYPE(polaw) :: pw
    INTEGER :: label, il, iw, ipos, iil
    REAL(kind=DP) :: tfcheck, totalperiod

    LOGICAL, PARAMETER :: direct_access = .true.

     write(stdout,*) 'VALUE TF', tf%l_fft_timefreq!ATTENZIONE
!first dealloacate and set

    totalperiod=2.d0*fftd%period+2.d0*fftd%period/real(fftd%n)


!read in  times/frequencies

!    do iil=-fftd%n,fftd%n!loop on time/frequency the order is the physical one
! we take advantage of the symmetry t ==> -t 
     do iil=0,fftd%n

!the following is in order to avoid same processor working with the same polaw
       il=mpime+iil
!       if(il>fftd%n) il=il-2*fftd%n-1
! we take advantage of the symmetry t ==> -t
       if(il>fftd%n) il=il-fftd%n-1
!we take care of the symmetry t ==> -t
       if(.not.direct_access) then
          CALL read_polaw(il,pw,debug,.false.)
       else
          CALL read_polaw_range(il,pw,debug,fftd%firstrow,fftd%firstrow+fftd%numrows-1,.true.)
       endif

!consistency test

      if( pw%numpw /= fftd%numpw ) then
         write(stdout,*) 'routine save_fft_data: consistency failed'
         stop
      endif

!check if ontime does not correspond
      if(pw%ontime .NEQV. fftd%ontime) then!update
         pw%ontime = fftd%ontime
         if(tf%l_fft_timefreq) then
            if(pw%ontime) then
               pw%time=fftd%period/real(fftd%n)*real(il)
            else
               pw%time=(2.d0*pi/totalperiod)*real(il)
            endif
         else
            if(pw%ontime) then
               pw%time=tf%times(il)
            else
               pw%time=tf%freqs(il)
            endif
         endif
      endif

!put in data at the right position
      if(tf%l_fft_timefreq) then
         if(.not.fftd%ontime) then !freqeuncy case 
            if(il>=0) then
               ipos=il+1
            else
               ipos=2*fftd%n+2+il+1
            endif
         else !time case 
            if(il>0) then
               ipos=2*fftd%n+2-il+1
            else
               ipos=-il+1
            endif
         endif
      else
!we take care of the t ==> -t symmetry
!         ipos=il+tf%n+2
         ipos = il +1
      endif
      do iw=1,fftd%numrows
         pw%pw(1:fftd%numpw,fftd%firstrow+iw-1) = fftd%fd(1:fftd%numpw,iw,ipos)
      enddo

!write on disk
      if(.not.direct_access) then
         CALL write_polaw(pw,debug)
      else
         CALL write_polaw_range(pw,debug,fftd%firstrow,fftd%firstrow+fftd%numrows-1,.true.)
      endif
      call mp_barrier( world_comm )
   enddo


   CALL free_memory_polaw(pw)

  END SUBROUTINE



  SUBROUTINE save_fft_data2(tf, fftd,debug)
! this subroutine writes the descriptor for the fftw 
! on the  polaw on disk
!data is put on appropiate order for FFT
!total period=2*T+T/n
!if we are ON TIME: THE ORDER IS REVERSED
    USE  polarization,        ONLY : polaw,read_polaw,write_polaw, free_memory_polaw, &
                                 &read_polaw_range,write_polaw_range
    USE  io_global,           ONLY : stdout
    USE  constants,           ONLY : eps8, pi
    USE  mp,                  ONLY : mp_barrier
    USE  mp_world,            ONLY : mpime, nproc, world_comm!group
    USE  times_gw,            ONLY : times_freqs
    USE  parallel_include


    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
    TYPE(fft_data) :: fftd!structure to be written
    LOGICAL :: debug!if true formatted files

    TYPE(polaw) :: pw
    INTEGER :: il, iw, ipos, numrows_read, nblk_siz, k, ierr
    INTEGER :: nbegin, nend, nbegin_ip, nend_ip, ip
    REAL(kind=DP) :: tfcheck, totalperiod
    COMPLEX(kind=DP), ALLOCATABLE :: rcvbuf( : ), sndbuf( : )

    LOGICAL, PARAMETER :: direct_access = .true.

     write(stdout,*) 'VALUE TF', tf%l_fft_timefreq!ATTENZIONE
!first dealloacate and set

    totalperiod=2.d0*fftd%period+2.d0*fftd%period/real(fftd%n)

    numrows_read = fftd%lastrow - fftd%firstrow + 1

    nblk_siz = ( fftd%n + 1) / nproc
    if( MOD( ( fftd%n + 1), nproc ) /= 0 ) nblk_siz = nblk_siz + 1

    nbegin = 0 + mpime * nblk_siz
    nend   = nbegin + nblk_siz - 1

    ALLOCATE( sndbuf( fftd%numpw * fftd%numrows * nproc ) )
    ALLOCATE( rcvbuf( fftd%numpw * fftd%numrows * nproc ) )

    ALLOCATE( pw%pw( fftd%numpw, numrows_read ) )
    pw%numpw = fftd%numpw

    DO il =  nbegin, nend

      pw%ontime = fftd%ontime
      if(tf%l_fft_timefreq) then
            if(pw%ontime) then
               pw%time=fftd%period/real(fftd%n)*real(il)
            else
               pw%time=(2.d0*pi/totalperiod)*real(il)
            endif
      else
            if(pw%ontime) then
               pw%time=tf%times(il)
            else
               pw%time=tf%freqs(il)
            endif
      endif

      do ip = 0, nproc - 1

         nbegin_ip = 0 + ip * nblk_siz

         ipos = il - nbegin + nbegin_ip

         if( ipos <= fftd%n ) then
            do iw = 1, fftd%numrows
               do k = 1, fftd%numpw
                  sndbuf( k + fftd%numpw * (iw-1) + fftd%numrows * fftd%numpw * ip ) = fftd%fd( k, iw, ipos + 1 )
               end do
            end do
         else
            do iw = 1, fftd%numrows
               do k = 1, fftd%numpw
                  sndbuf( k + fftd%numpw * (iw-1) + fftd%numrows * fftd%numpw * ip ) = 0.0d0 
               end do
            end do
         end if

      enddo

#if defined(__MPI)
      CALL MPI_ALLTOALL( sndbuf, fftd%numrows * fftd%numpw, MPI_DOUBLE_COMPLEX,  &
                         rcvbuf, fftd%numrows * fftd%numpw, MPI_DOUBLE_COMPLEX, world_comm, ierr )
#else
      rcvbuf(:)=sndbuf(:)
#endif


      if( il <= fftd%n ) then

         do iw = 1, numrows_read
            do k = 1, fftd%numpw
               pw%pw( k, iw ) = dble(rcvbuf( k + fftd%numpw * ( iw - 1 ) ) )
            end do
         end do

         

         pw%label = il

         pw%factor=fftd%factor
         CALL write_polaw_range( pw, debug, fftd%firstrow, fftd%lastrow, .false. )

      end if

   enddo

   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )

   CALL free_memory_polaw( pw )

  END SUBROUTINE



  SUBROUTINE transform_fft_data(fftd)
!this subroutine performs a FFT transform from imaginary time 
!to imaginary frequency and viceversa
!uses FFTW machinery
!does not reorder data but puts appropriate factors
   USE  constants,           ONLY :  pi
   USE  fft_scalar,          ONLY : cft_1z
   USE  fft_support,         ONLY : good_fft_order


  implicit none
  TYPE(fft_data) :: fftd!structure to be transformed



  INTEGER :: iw,jw, il,kw
  COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: in,out!temporary arrays
  INTEGER*8  :: plan
  REAL(kind=DP) :: omega,time,totalperiod,totalfrequency
  COMPLEX(kind=DP) :: fact
  INTEGER :: good_dim


  good_dim = good_fft_order(2*fftd%n+2)


  allocate(in(good_dim),out(good_dim))
  
  totalperiod=2.d0*fftd%period+2.d0*fftd%period/real(fftd%n)
  totalfrequency=(2.d0*pi/totalperiod)*real(2*fftd%n+2)
  
  if(fftd%ontime) then!time to frequency transform
    fftd%ontime=.false.

    do iw=1,fftd%numrows
      do jw=1,fftd%numpw
        in(1:2*fftd%n+2)=fftd%fd(jw,iw,1:2*fftd%n+2)


        call cft_1z(in,1,2*fftd%n+2,good_dim, -1,out)

        fftd%fd(jw,iw,1:2*fftd%n+2)=out(1:2*fftd%n+2)*dble(2*fftd%n+2) !ATTENZIONE
      enddo
    enddo

!set appropriate factors
   do il=0,2*fftd%n+2-1
      if(il <= (2*fftd%n+1)) then
        omega=(2.d0*pi/totalperiod)*real(il)
      else
        omega=(2.d0*pi/totalperiod)*real(il-2*fftd%n-2)
      endif
      !fact=exp((0.d0,-1.d0)*omega*totalperiod/2.d0)*(0.d0,-1.d0)*(fftd%period/real(fftd%n))
      fact=exp((0.d0,-1.d0)*omega*totalperiod/2.d0)*(fftd%period/real(fftd%n))
      fftd%fd(:,:,il+1)=fftd%fd(:,:,il+1)*fact
   enddo
   fftd%factor=fftd%factor*(0.d0,-1.d0)
  else!frequency to time transform
!alternative approach
    fftd%ontime=.true.

    do iw=1,fftd%numrows
      do jw=1,fftd%numpw
        in(1:2*fftd%n+2)=fftd%fd(jw,iw,1:2*fftd%n+2)

        call cft_1z(in,1,2*fftd%n+2,good_dim, 1,out)

        fftd%fd(jw,iw,1:2*fftd%n+2)=out(1:2*fftd%n+2)
      enddo
    enddo

!set appropriate factors
   do il=0,2*fftd%n+2-1
      if(il<= (2*fftd%n+1)) then
       time=(fftd%period/real(fftd%n))*real(il)
      else
        time=(fftd%period/real(fftd%n))*real(il-2*fftd%n-2)
      endif
!      fact=exp((0.d0,+1.d0)*time*totalfrequency/2.d0)*(0.d0,+1.d0)/totalperiod
      fact=exp((0.d0,+1.d0)*time*totalfrequency/2.d0)/totalperiod
      fftd%fd(:,:,il+1)=fftd%fd(:,:,il+1)*fact
   enddo
   fftd%factor=fftd%factor*(0.d0,+1.d0)

  endif

  deallocate(in,out)
  return
  END SUBROUTINE



  SUBROUTINE transform_fft_data_grid(tf, fftd)
!this subroutine performs a Fourier transform from imaginary time
!to imaginary frequency and viceversa
!uses user defined grids
   USE  constants,           ONLY :  pi
   USE  times_gw,            ONLY : times_freqs
   USE  io_global,           ONLY : stdout

   implicit none

   TYPE(times_freqs) :: tf! time frequency grids and factors
   TYPE(fft_data) :: fftd!structure to be transformed

   COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: fd_new
   COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: factors
   COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: tmpc
   INTEGER :: ii,jj, iw, jw

   INTEGER, PARAMETER :: nmesh=50
   INTEGER, PARAMETER :: nmesh_g=50
   INTEGER, PARAMETER :: nn=2

   REAL(kind=DP) :: b_p,b_m,r_p,r_m
   COMPLEX(kind=DP) :: a_p,a_m, cor_1,cor_2
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
   COMPLEX(kind=DP), DIMENSION(nmesh_g) :: tmpg
   COMPLEX(kind=DP), ALLOCATABLE ::  fij(:,:), fp(:),fm(:)
   
!we take advantage of the t ==> -t symmetry
  
!   allocate(fd_new(fftd%numpw,fftd%numrows,2*fftd%n+1))
   allocate(fd_new(fftd%numpw,fftd%numrows,fftd%n+1))
   allocate(factors(-tf%n:tf%n), tmpc(fftd%numpw,-tf%n:tf%n))
 
!check for consistency
   if(fftd%n /= tf%n) then
      write(stdout,*) 'Routine transform_fft_data_grid: consistency failed'
      stop
   endif

   if(fftd%ontime) then!time to frequency transform
      fftd%factor=fftd%factor*(0.d0,-1.d0)
   else
      fftd%factor=fftd%factor*(0.d0,+1.d0)
   endif


!we take care of the t ==> -t symmetry
!   do ii=-tf%n,tf%n
   do ii=0,tf%n
      if(fftd%ontime) then!time to frequency transform
         do jj=-tf%n,tf%n
            factors(jj)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs(ii)*tf%times(jj))
         enddo
         !factors(:)=factors(:)*(0.d0,-1.d0)
      else!frequency to time transform
         do jj=-tf%n,tf%n
            factors(jj)=tf%weights_freq(jj)*exp((0.d0,1.d0)*tf%times(ii)*tf%freqs(jj))
         enddo
         !factors(:)=factors(:)*(0.d0,1.d0)/(2.d0*pi)
         factors(:)=factors(:)/(2.d0*pi)
      endif
      do jw=1,fftd%numrows
!we take care of the t ==> -t symmetry
!         do jj=-tf%n,tf%n
         do jj=0,tf%n
!            tmpc(:,jj)=fftd%fd(:,jw,jj+tf%n+2)*factors(jj)
            tmpc(:,jj)=fftd%fd(:,jw,jj+1)*factors(jj)
         enddo
         do jj=-tf%n,-1
!            tmpc(:,jj)=fftd%fd(:,jw,abs(jj)+tf%n+2)*factors(jj)
            tmpc(:,jj)=fftd%fd(:,jw,abs(jj)+1)*factors(jj)
         enddo
         do iw=1,fftd%numpw
!            fd_new(iw,jw,ii+tf%n+1)=sum(tmpc(iw,-tf%n:tf%n))
            fd_new(iw,jw,ii+1)=sum(tmpc(iw,-tf%n:tf%n))
         enddo
      enddo
   enddo
   if(fftd%ontime .and. tf%l_fourier_fit_time) then

      allocate(fij(-tf%n:tf%n,nmesh))
      allocate(fp(nmesh),fm(nmesh))
      allocate(x(nmesh),w(nmesh))
      x(:)=0.d0
      w(:)=0.d0
      call legzo(nmesh,x,w)
                    
      x(:)=x(:)*tf%tau/2.d0
      x(:)=x(:)+tf%tau/2.d0
      w(:)=w(:)*tf%tau/2.d0
      
      !x(:)=x(:)*(tf%times(tf%n)-tf%tau)/2.d0
      !x(:)=x(:)+(tf%times(tf%n)-tf%tau)/2.d0+tf%tau
      !w(:)=w(:)*(tf%times(tf%n)-tf%tau)/2.d0
         


      do jj=1,nmesh
         write(stdout,*)'MESH', jj, x(jj),w(jj)
      enddo
      
      do ii=-tf%n,tf%n
         do jj=1,nmesh
            fij(ii,jj)=exp((0.d0,-1.d0)*tf%freqs(ii)*x(jj))
         enddo
      enddo
      
      do iw=1,fftd%numpw
         do jw=1,fftd%numrows
            r_p=dble(fftd%fd(iw,jw,2*tf%n+1)/fftd%fd(iw,jw,2*tf%n+2))
            if(r_p <= 1.d0) r_p = tf%g_tau
            b_p=log(r_p)/(tf%times(tf%n)-tf%times(tf%n-1))
            a_p=fftd%fd(iw,jw,2*tf%n+1)/(exp(-b_p*tf%times(tf%n-1)))
            if(r_p == tf%g_tau) a_p=0.d0
            
            r_m=dble(fftd%fd(iw,jw,3)/fftd%fd(iw,jw,2))
            if(r_m <= 1.d0) r_m = tf%g_tau
            b_m=log(r_m)/(tf%times(-tf%n+1)-tf%times(-tf%n))
            a_m=fftd%fd(iw,jw,3)/(exp(b_m*tf%times(-tf%n+1)))
            if(r_m == tf%g_tau) a_m=0.d0
            
            do jj=1,nmesh
               fp(jj)=a_p*exp(-b_p*x(jj))*w(jj)
               fm(jj)=a_m*exp(-b_m*x(jj))*w(jj)
            enddo

            do ii=-tf%n,tf%n
               !cor_1=(0.d0,-1.d0)*(a_p/(b_p+(0.d0,1.d0)*tf%freqs(ii)))+&
               !     &(0.d0,-1.d0)*(a_m/(b_m-(0.d0,1.d0)*tf%freqs(ii)))
               cor_1=(0.d0,-1.d0)*2.d0*a_p*b_p/(b_p**2.d0+tf%freqs(ii)**2.d0)
               cor_2=0.d0
               do jj=1,nmesh
                  cor_2=cor_2-fij(ii,jj)*fp(jj)
                  cor_2=cor_2-conjg(fij(ii,jj))*fp(jj)
               enddo
               cor_2=cor_2*(0.d0,-1.d0)
               if(ii==0) write(stdout,*) 'COR', cor_1,cor_2
               fd_new(iw,jw,ii+tf%n+1)=fd_new(iw,jw,ii+tf%n+1)+cor_1+cor_2
            enddo
         enddo
      enddo
      deallocate(fij,fp,fm)
      deallocate(x,w)
   else if(.not.fftd%ontime .and. tf%l_fourier_fit_freq) then
      
      allocate(fij(-tf%n:tf%n,nmesh_g))
      allocate(fp(nmesh_g),fm(nmesh_g))
      allocate(x(nmesh_g),w(nmesh_g))
      x(:)=0.d0
      w(:)=1.d0
      do jj=1,nmesh
         write(stdout,*)'MESH', jj, x(jj),w(jj)
      enddo

      x(:)=0.d0
      w(:)=0.d0
      call legzo(nmesh_g,x,w)
      do jj=1,nmesh
         write(stdout,*)'MESH', jj, x(jj),w(jj)
      enddo

      x(:)=x(:)*(tf%omega/2.d0)*dble(nn)
      x(:)=x(:)+(tf%omega/2.d0)*dble(nn)+tf%omega
      w(:)=w(:)*(tf%omega/2.d0)*dble(nn)

      do ii=-tf%n,tf%n
         do jj=1,nmesh_g
            fij(ii,jj)=exp((0.d0,1.d0)*tf%times(ii)*x(jj))*w(jj)
         enddo
      enddo
      
      do jj=1,nmesh
         write(stdout,*)'MESH', jj, x(jj),w(jj)
      enddo


      do iw=1,fftd%numpw
         do jw=1,fftd%numrows
            
            r_p=dble(fftd%fd(iw,jw,2*tf%n+1)/fftd%fd(iw,jw,2*tf%n+2))
            b_p=(tf%freqs(tf%n)**2.d0-r_p*tf%freqs(tf%n-1)**2.d0)/(r_p-1.d0)
            if(b_p < -tf%freqs(tf%n-1)**2.d0) b_p=-tf%g_omega*tf%freqs(tf%n-1)**2.d0
            a_p=fftd%fd(iw,jw,2*tf%n+1)*(b_p+tf%freqs(tf%n-1)**2.d0)
            
            r_m=dble(fftd%fd(iw,jw,3)/fftd%fd(iw,jw,2))
            b_m=(tf%freqs(-tf%n)**2.d0-r_m*tf%freqs(-tf%n+1)**2.d0)/(r_m-1.d0)
            if(b_m < -tf%freqs(-tf%n+1)**2.d0) b_m=-tf%g_omega*tf%freqs(-tf%n+1)**2.d0
            a_m=fftd%fd(iw,jw,3)*(b_m+tf%freqs(-tf%n+1)**2.d0)
            
            do jj=1,nmesh_g
               fp(jj)=a_p/(b_p+x(jj)**2.d0)
               fm(jj)=a_m/(b_m+x(jj)**2.d0)
            enddo

            do ii=-tf%n,tf%n
               do jj=1,nmesh_g
                  tmpg(jj)=fij(ii,jj)*fp(jj)
               enddo
               cor_1=sum(tmpg(:))
               do jj=1,nmesh_g
                  tmpg(jj)=conjg(fij(ii,jj))*fm(jj)
               enddo
               cor_2=sum(tmpg(:))
               cor_1=cor_1*(0.d0,+1.d0)/(2.d0*pi)
               cor_2=cor_2*(0.d0,+1.d0)/(2.d0*pi)
               fd_new(iw,jw,ii+tf%n+1)=fd_new(iw,jw,ii+tf%n+1)+cor_1+cor_2
            enddo
         enddo
      enddo
      deallocate(x,w)
      deallocate(fij,fp,fm)

   endif
      
   if(fftd%ontime) then
      fftd%ontime=.false.
   else
      fftd%ontime=.true.
   endif

!we take care of the t ==> -t symmetry
!   do jj=-tf%n,tf%n
   do jj=0,tf%n
!      fftd%fd(:,:,jj+tf%n+2)=fd_new(:,:,jj+tf%n+1)
      fftd%fd(:,:,jj+1)=fd_new(:,:,jj+1)
   enddo
   



   deallocate(fd_new)
   deallocate(factors)
   deallocate(tmpc)

   return

 END SUBROUTINE transform_fft_data_grid



 END MODULE fft_gw
     
   
