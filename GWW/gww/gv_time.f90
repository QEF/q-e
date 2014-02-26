!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


MODULE w_divergence
!for the treatment of the G=0,G=0 divergence of the W operator
  USE kinds, ONLY : DP

  TYPE gv_time
!this structure contains the data for the treatment of the w_divergence
    INTEGER :: n!number of time/frequency steps
    REAL(kind=DP) :: omega!frequency range
    REAL(kind=DP) :: tau!time range
    INTEGER :: max_i!number of states
    COMPLEX(kind=DP), DIMENSION(:,:),POINTER :: ex!terms <Psi_i|iG(r,r')v(r,r')|Psi_i>
    LOGICAL :: ontime!if .true. is on imaginary time, otherwise frequency
    COMPLEX(kind=DP), DIMENSION(:), POINTER :: inv_epsi!head (G=0,G=0) of the inverse dielectric matrix
 END TYPE gv_time

 CONTAINS

   SUBROUTINE initialize_gv_time(gt)
!initialize
     implicit none

     TYPE(gv_time) :: gt
     
     nullify(gt%ex)
     nullify(gt%inv_epsi)
     return
   END SUBROUTINE initialize_gv_time

   SUBROUTINE free_memory_gv_time(gt)
     implicit none                                            
                                                            
     TYPE(gv_time) :: gt                                      
         
     if(associated(gt%ex)) deallocate(gt%ex)
     if(associated(gt%inv_epsi)) deallocate(gt%inv_epsi)
     nullify(gt%ex)  
     return
   END SUBROUTINE free_memory_gv_time

!read data from PW from file
   SUBROUTINE read_data_pw_gv_time(gt, prefix)

     USE io_global,            ONLY : stdout, ionode, ionode_id
     USE mp,                   ONLY : mp_bcast
     USE mp_world,             ONLY : world_comm
     USE io_files,             ONLY : tmp_dir

     implicit none

     INTEGER, EXTERNAL :: find_free_unit
     TYPE(gv_time) :: gt!to be read from PW file
     CHARACTER(LEN=256) ::  prefix!to designate the PW files
     
     REAL(kind=DP), ALLOCATABLE :: buf(:)
     INTEGER :: iun,i


     if(ionode) then
        iun = find_free_unit()
        open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'.gv_time', status='old',form='unformatted')
        read(iun) gt%max_i
        read(iun) gt%n
        read(iun) gt%tau
     endif

     call mp_bcast(gt%max_i, ionode_id,world_comm)
     call mp_bcast(gt%n, ionode_id,world_comm)
     call mp_bcast(gt%tau, ionode_id,world_comm)
     
     allocate(gt%ex(gt%max_i,2*gt%n+2))
     allocate(gt%inv_epsi(2*gt%n+1))
     gt%inv_epsi(:)=0.d0

     if(ionode) then
        allocate(buf(gt%max_i))
        do i=1,2*gt%n+2
           read(iun) buf(1:gt%max_i)
           gt%ex(:,i)=cmplx(buf(:),0.d0)
        enddo
        close(iun)
        deallocate(buf)
     endif

     call mp_bcast(gt%ex(:,:),ionode_id,world_comm)

     return
   END SUBROUTINE read_data_pw_gv_time



   SUBROUTINE write_gv_time(gt)
!save on file
     USE io_global,            ONLY : ionode
     USE io_files,             ONLY : prefix,tmp_dir
         
     implicit none
     INTEGER, EXTERNAL :: find_free_unit
     TYPE(gv_time) :: gt

     INTEGER :: iun,i
     
     if(ionode) then
         iun = find_free_unit()
         open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'gv_time', status='unknown',form='unformatted')

         write(iun) gt%n
         write(iun) gt%omega
         write(iun) gt%tau
         write(iun) gt%max_i
         write(iun) gt%ontime
         do i=1,2*gt%n+2
            write(iun) gt%ex(1:gt%max_i,i)
         enddo
         write(iun) gt%inv_epsi(1:2*gt%n+1)
         close(iun)
      end if
      return
   END SUBROUTINE write_gv_time

   SUBROUTINE read_gv_time(gt)
!read from file
     USE io_global,            ONLY : ionode, ionode_id
     USE mp,                   ONLY : mp_bcast
     USE mp_world,             ONLY : world_comm
     USE io_files,             ONLY : prefix,tmp_dir
    

     implicit none
     INTEGER, EXTERNAL :: find_free_unit
     TYPE(gv_time) :: gt

     INTEGER :: iun,i

     if(ionode) then
         iun = find_free_unit()
         open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'gv_time', status='old',form='unformatted')

         read(iun) gt%n
         read(iun) gt%omega
         read(iun) gt%tau
         read(iun) gt%max_i
         read(iun) gt%ontime
      endif

      call mp_bcast(gt%n, ionode_id,world_comm)
      call mp_bcast(gt%omega, ionode_id,world_comm)
      call mp_bcast(gt%tau, ionode_id,world_comm)
      call mp_bcast(gt%max_i, ionode_id,world_comm)
      call mp_bcast(gt%ontime, ionode_id,world_comm)

      allocate(gt%ex(gt%max_i,2*gt%n+2))
      allocate(gt%inv_epsi(2*gt%n+1))


      if(ionode) then
         do i=1,2*gt%n+2
            read(iun) gt%ex(1:gt%max_i,i)
         enddo
         read(iun) gt%inv_epsi(1:2*gt%n+1)
         close(iun)
      endif
      
      call mp_bcast(gt%ex(:,:), ionode_id,world_comm)
      call mp_bcast(gt%inv_epsi(:), ionode_id,world_comm)

      return
    END SUBROUTINE read_gv_time


    SUBROUTINE fft_gv_time(gt,tf)
!performs fft transform of gv inv_epsi data

      USE  constants,           ONLY :  pi
      USE  times_gw,            ONLY : times_freqs
      USE io_global,            ONLY : stdout

      implicit none
      
      TYPE(gv_time) :: gt
      TYPE(times_freqs) :: tf! time frequency grids and factors


      INTEGER :: ii,jj,iw,jw
      COMPLEX(kind=DP), ALLOCATABLE :: tmpc(:), factors(:)
      COMPLEX(kind=DP), ALLOCATABLE :: inv_epsi_new(:)
      

      allocate(factors(-tf%n:tf%n), tmpc(-tf%n:tf%n))
      allocate(inv_epsi_new(2*tf%n+1))

!check for consistency

      if(tf%n /= gt%n) then
         write(stdout,*) 'FFT_GV: not consistent n'
         stop
      endif
      
      if(tf%omega /= gt%omega) then
          write(stdout,*) 'FFT_GV: not consistent omega'
         stop
      endif

      if(tf%tau /= gt%tau) then
         write(stdout,*) 'FFT_GV: not consistent tau'
         stop
      endif




     do ii=-tf%n,tf%n
        write(*,*) 'ATTENZIONE',ii,gt%max_i!ATTENZIONE

      if(gt%ontime) then!time to frequency transform
         do jj=-tf%n,tf%n
            factors(jj)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs(ii)*tf%times(jj))
         enddo
      else!frequency to time transform
         do jj=-tf%n,tf%n
            factors(jj)=tf%weights_freq(jj)*exp((0.d0,1.d0)*tf%times(ii)*tf%freqs(jj))
         enddo
         factors(:)=factors(:)/(2.d0*pi)
      endif
      write(*,*) 'ATTENZIONE2',ii!ATTENZIONE


      do jj=-tf%n,tf%n
         tmpc(jj)=gt%inv_epsi(jj+tf%n+1)*factors(jj)
      enddo
      
      
      inv_epsi_new(ii+tf%n+1)=sum(tmpc(-tf%n:tf%n))
      


   enddo

   write(stdout,*) 'ATTENZIONE3'

!add i factor
   if(gt%ontime) then
      gt%ontime=.false.
      gt%inv_epsi(1:2*gt%n+1)=(0.d0,-1.d0)*inv_epsi_new(1:2*gt%n+1)
   else
      gt%ontime=.true.
      gt%inv_epsi(1:2*gt%n+1)=(0.d0,1.d0)*inv_epsi_new(1:2*gt%n+1)
   endif

   write(stdout,*) 'ATTENZIONE4'

   deallocate(factors,tmpc,inv_epsi_new)

   write(stdout,*) 'ATTENZIONE5'

   return
 end SUBROUTINE fft_gv_time

 SUBROUTINE setup_gv_time(gt)
!this subroutine set up the gv_time structure 
!with the head of the inverse dielectric matrices
!to be done in imaginary time

   USE io_global,  ONLY : stdout
   
   implicit none

   TYPE(gv_time) :: gt!the structure to be set up
   INTEGER :: ii,it

   
   if(.not.gt%ontime) then
      write(stdout,*) 'Routine setup_gv_time imaginary time required'
      stop
   endif

   do it=1,gt%n
      gt%ex(:,it)=gt%ex(:,it)*gt%inv_epsi(it)
   enddo

   do it=gt%n+2,2*gt%n+1
      gt%ex(:,it)=gt%ex(:,it)*gt%inv_epsi(it)
   enddo

!now the t=0 term

   gt%ex(:,gt%n+1)=0.5d0*gt%inv_epsi(gt%n+1)*(gt%ex(:,gt%n+1)+gt%ex(:,2*gt%n+2))

   
   return
 END SUBROUTINE setup_gv_time

 END MODULE w_divergence
