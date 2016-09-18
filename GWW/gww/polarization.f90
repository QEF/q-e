!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

  MODULE polarization
!this module describes the structure for the polarization P
! and dressed iteraction W, imaginary time/frequency
    USE kinds, ONLY : DP

      TYPE polaw
!this structure describe a generic P or W function
!usually in the space of orthonormalized products of wanniers
        INTEGER :: label!label to read/write to disk
        LOGICAL :: ontime!if .true. is on imaginary time, otherwise frequency
        REAL(kind=DP) :: time!imaginary time or frequency
        INTEGER :: numpw!number of states (products of wanniers)
        REAL(kind=DP), DIMENSION(:,:), POINTER :: pw!the P or W 
        COMPLEX(kind=DP) :: factor!complex factor to be multiplied to the real matrix pw
      END TYPE polaw



  CONTAINS

    SUBROUTINE initialize_polaw(pw)
!this subroutine initializes polaw
      implicit none
      TYPE(polaw) :: pw
      nullify(pw%pw)
      return
    END SUBROUTINE initialize_polaw

    SUBROUTINE free_memory_polaw(pw)
!this subroutine deallocates the green descriptor
      implicit none
      TYPE(polaw) pw 
      if(associated(pw%pw)) deallocate(pw%pw)
      nullify(pw%pw)
      return
    END SUBROUTINE

    SUBROUTINE conjugate_polaw(pw)
! this subroutine calculculates the conjugate of the polaw matrix
      implicit none

      TYPE(polaw) pw
      pw%label=-pw%label
      return
    END SUBROUTINE conjugate_polaw


   SUBROUTINE write_polaw(pw,debug)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(polaw) :: pw!the green function to be written
    LOGICAL :: debug! if .true. produces formatted output

    LOGICAL :: direct_file = .true.!if true uses direct- access file to write matrix on disk

    INTEGER :: iw, jw, iung
    CHARACTER(5) :: nfile

    if(pw%label >= 0 ) then
      write(nfile,'(5i1)') &
        & pw%label/10000,mod(pw%label,10000)/1000,mod(pw%label,1000)/100,mod(pw%label,100)/10,mod(pw%label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='unknown',form='unformatted')
      else 
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='unknown',form='formatted')
      endif
    else
      write(nfile,'(5i1)') &         
        & -pw%label/10000,mod(-pw%label,10000)/1000,mod(-pw%label,1000)/100,mod(-pw%label,100)/10,mod(-pw%label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='unknown',form='formatted')
      endif
    endif
    if(.not.debug) then
      write(iung) pw%label
      write(iung) pw%ontime
      write(iung) pw%time
      write(iung) pw%numpw
      write(iung) pw%factor
      if(.not. direct_file) then
         do iw=1,pw%numpw
            write(iung)  pw%pw(1:pw%numpw,iw)
         enddo
      endif
    else
      write(iung,*) pw%label
      write(iung,*) pw%ontime
      write(iung,*) pw%time
      write(iung,*) pw%numpw
      write(iung,*) pw%factor
      if(.not. direct_file) then
         do iw=1,pw%numpw
            do jw=1,pw%numpw
               write(iung,*)  pw%pw(jw,iw)
            enddo
         enddo
      endif
    endif

    close(iung)
    
    if(direct_file) then
       iung = find_free_unit()
       if(pw%label >= 0 ) then
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       else
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.-'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       endif
       do iw=1,pw%numpw
          write(unit=iung, rec=iw) pw%pw(:,iw)
       enddo
       close(iung)
    endif
       

    return
   END SUBROUTINE

   SUBROUTINE  read_polaw(label, pw,debug,l_verbose)
!this subroutine reads the green function from disk
!the file name is taken from the label

    USE io_files,             ONLY : prefix,tmp_dir
    USE io_global,            ONLY : stdout
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(polaw) :: pw!the green function to be read
    INTEGER :: label! the label identifing the required green function
    LOGICAL :: debug!if true formatted files
    LOGICAL, INTENT(in) :: l_verbose

    LOGICAL :: direct_file = .true.!if true uses direct- access file to read matrix from disk

    INTEGER :: iw, jw, iung
    CHARACTER(5) :: nfile

    if(l_verbose) write(stdout,*) 'Read polaw'!ATTENZIONE
!first deallocate
    call free_memory_polaw(pw)
     if(l_verbose) write(stdout,*) 'Read polaw2'!ATTENZIONE
    if(label >= 0 ) then
      write(nfile,'(5i1)') label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='old',form='formatted')
      endif
    else
      write(nfile,'(5i1)') -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='old',form='formatted')
      endif
    endif
    if(.not.debug) then
      read(iung) pw%label
      read(iung) pw%ontime
      read(iung) pw%time
      read(iung) pw%numpw
      read(iung) pw%factor
    else
      read(iung,*) pw%label
      read(iung,*) pw%ontime
      read(iung,*) pw%time
      read(iung,*) pw%numpw
      read(iung,*) pw%factor
    endif
   write(stdout,*) 'Read polaw',pw%numpw!ATTENZIONE
!now allocate
    allocate(pw%pw(pw%numpw,pw%numpw))
    if(.not. direct_file) then
       if(.not.debug) then
          do iw=1,pw%numpw
             read(iung)  pw%pw(1:pw%numpw,iw)
          enddo
       else
          do iw=1,pw%numpw
             do jw=1,pw%numpw
                read(iung,*)  pw%pw(jw,iw)
             enddo
          enddo
       endif
    endif
    close(iung)
     if(l_verbose) write(stdout,*) 'Read polaw4'!ATTENZIONE
    if(direct_file) then
       iung = find_free_unit()
       if(label >= 0 ) then
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       else
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.-'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       endif
       if(l_verbose) write(stdout,*) 'Read polaw5'!ATTENZIONE
       do iw=1,pw%numpw
          read(unit=iung, rec=iw) pw%pw(:,iw)
       enddo
       close(iung)
    endif
    if(l_verbose) write(stdout,*) 'Read polaw6'!ATTENZIONE

    return
   END SUBROUTINE


   SUBROUTINE  read_polaw_global(label, pw)
!this subroutine reads the green function from disk 
!the file name is taken from the label
!the ionode_id distribute to all the processors

    USE io_files,             ONLY : prefix,tmp_dir
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE mp,                   ONLY : mp_barrier, mp_bcast
    USE mp_world,             ONLY : world_comm

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(polaw) :: pw!the green function to be read
    INTEGER :: label! the label identifing the required green function 
   

    LOGICAL :: direct_file = .true.!if true uses direct- access file to read matrix from disk

    INTEGER :: iw, jw, iung
    CHARACTER(5) :: nfile


!first deallocate
    call free_memory_polaw(pw)
    if(ionode) then
       if(label >= 0 ) then
          write(nfile,'(5i1)') label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
          iung = find_free_unit()
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='old',form='unformatted')
       else
          write(nfile,'(5i1)') -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
          iung = find_free_unit()
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='old',form='unformatted')
       endif
       
   
       read(iung) pw%label
       read(iung) pw%ontime
       read(iung) pw%time
       read(iung) pw%numpw
       read(iung) pw%factor
    endif
    call mp_bcast(pw%label,ionode_id,world_comm)
    call mp_bcast(pw%ontime,ionode_id,world_comm)
    call mp_bcast(pw%time,ionode_id,world_comm)
    call mp_bcast(pw%numpw,ionode_id,world_comm)
    call mp_bcast(pw%factor,ionode_id,world_comm)
 
    allocate(pw%pw(pw%numpw,pw%numpw))
    if(ionode) then
       if(.not. direct_file) then
          do iw=1,pw%numpw
             read(iung)  pw%pw(1:pw%numpw,iw)
          enddo
       endif
       close(iung)

       if(direct_file) then
          iung = find_free_unit()
          if(label >= 0 ) then
             open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.'// nfile, &
                  &status='unknown',recl=pw%numpw*DP,access='direct')
          else
             open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.-'// nfile, &
                  &status='unknown',recl=pw%numpw*DP,access='direct')
          endif
          
          do iw=1,pw%numpw
             read(unit=iung, rec=iw) pw%pw(:,iw)
          enddo
          close(iung)
       endif
    endif
    do iw=1,pw%numpw
       call mp_barrier( world_comm )
       call mp_bcast(pw%pw(:,iw),ionode_id,world_comm)
    enddo

    return
  END SUBROUTINE read_polaw_global




   SUBROUTINE write_polaw_range( pw, debug, range_min, range_max, full_range )
!this subroutine writes the green function on disk
!the file name is taken from the label
!writes column from range_min to range_max

    USE io_files,             ONLY : prefix,tmp_dir
    USE io_global,            ONLY : stdout


    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(polaw) :: pw!the green function to be written
    LOGICAL :: debug! if .true. produces formatted output
    INTEGER, INTENT(in) :: range_min, range_max!range of column
    LOGICAL, INTENT(IN) :: full_range


    LOGICAL :: direct_file = .true.!if true uses direct- access file to write matrix on disk

    INTEGER :: iw, jw, iung, iww
    CHARACTER(5) :: nfile

!check range

    if(range_min<1 .or. range_max> pw%numpw) then
       write(stdout,*) 'write_polaw_range: out of range = ', range_min, range_max
       stop
    endif



    if(pw%label >= 0 ) then
      write(nfile,'(5i1)') &
        & pw%label/10000,mod(pw%label,10000)/1000,mod(pw%label,1000)/100,mod(pw%label,100)/10,mod(pw%label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='unknown',form='formatted')
      endif
    else
      write(nfile,'(5i1)') &
        & -pw%label/10000,mod(-pw%label,10000)/1000,mod(-pw%label,1000)/100,mod(-pw%label,100)/10,mod(-pw%label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='unknown',form='formatted')
      endif
    endif
    if(.not.debug) then
      write(iung) pw%label
      write(iung) pw%ontime
      write(iung) pw%time
      write(iung) pw%numpw
      write(iung) pw%factor
      if(.not. direct_file) then
         do iw=1,pw%numpw
            write(iung)  pw%pw(1:pw%numpw,iw)
         enddo
      endif
    else
      write(iung,*) pw%label
      write(iung,*) pw%ontime
      write(iung,*) pw%time
      write(iung,*) pw%numpw
      write(iung,*) pw%factor
      if(.not. direct_file) then
         do iw=1,pw%numpw
            do jw=1,pw%numpw
               write(iung,*)  pw%pw(jw,iw)
            enddo
         enddo
      endif
    endif

    close(iung)

    if(direct_file) then
       iung = find_free_unit()
       if(pw%label >= 0 ) then
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       else
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.-'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       endif
       do iw=range_min,range_max
          iww = iw
          if( .not. full_range ) iww = iw - range_min + 1
          write(unit=iung, rec=iw) pw%pw(:,iww)
       enddo
       close(iung)
    endif


    return
  END SUBROUTINE write_polaw_range

   SUBROUTINE  read_polaw_range(label, pw,debug,range_min,range_max, full_range )
!this subroutine reads the green function from disk
!the file name is taken from the label
!reads columns from range_min to range_max

    USE io_files,             ONLY : prefix,tmp_dir
    USE io_global,            ONLY : stdout

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(polaw) :: pw!the green function to be read
    INTEGER :: label! the label identifing the required green function
    LOGICAL :: debug!if true formatted files
    INTEGER, INTENT(in) :: range_min, range_max!defining range
    LOGICAL, INTENT(IN) :: full_range


    LOGICAL :: direct_file = .true.!if true uses direct- access file to read matrix from disk

    INTEGER :: iw, jw, iung, iww
    CHARACTER(5) :: nfile




!check
    if(range_min<1 ) then
       write(stdout,*) 'read_polaw_range: out of range ', range_min, range_max
       stop
    endif



!first deallocate
    call free_memory_polaw(pw)


    if(label >= 0 ) then
      write(nfile,'(5i1)') label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.'// nfile, status='old',form='formatted')
      endif
    else
      write(nfile,'(5i1)') -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polaw.-'// nfile, status='old',form='formatted')
      endif
    endif
    if(.not.debug) then
      read(iung) pw%label
      read(iung) pw%ontime
      read(iung) pw%time
      read(iung) pw%numpw
      read(iung) pw%factor
    else
      read(iung,*) pw%label
      read(iung,*) pw%ontime
      read(iung,*) pw%time
      read(iung,*) pw%numpw
      read(iung,*) pw%factor
    endif
!now allocate


    if( full_range ) then
       allocate(pw%pw(pw%numpw,pw%numpw))
    else
       allocate( pw%pw( pw%numpw, range_max - range_min + 1 ) )
    endif

    if(.not. direct_file) then
       if(.not.debug) then
          do iw=1,pw%numpw
             read(iung)  pw%pw(1:pw%numpw,iw)
          enddo
       else
          do iw=1,pw%numpw
             do jw=1,pw%numpw
                read(iung,*)  pw%pw(jw,iw)
             enddo
          enddo
       endif
    endif
    close(iung)

    if(direct_file) then
!check
    if(range_max > pw%numpw ) then
       write(stdout,*) 'read_polaw_range: out of range = ', range_min, range_max
       stop
    endif

       iung = find_free_unit()
       if(label >= 0 ) then
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       else
          open( unit=iung, file=trim(tmp_dir)//trim(prefix)//'-'//'polawd.-'// nfile, &
               &status='unknown',recl=pw%numpw*DP,access='direct')
       endif
       do iw=range_min, range_max
          iww = iw
          if( .not. full_range ) iww = iw - range_min + 1
          read(unit=iung, rec=iw) pw%pw(:,iww)
       enddo
       close(iung)
    endif


    return
  END SUBROUTINE read_polaw_range





   SUBROUTINE create_polarization(time,pr,gf_p,gf_m,qm,debug)
!this subroutine set the polarization in imaginary time 
! as P(r,r',it)=G(r,r',it)*G(r,r',-it)
! for our basis 
! P_{i,j}=\sum_{l,n,m,o} Q_{i,lm}Q_{j,mo}G_{lm}(it)G_{no}(-it)
!THERE IS ALSO A SPIN FACTOR 2


   USE basic_structures,     ONLY : q_mat,wannier_P
   USE green_function,       ONLY : green
   USE io_global,            ONLY : stdout
   USE constants,            ONLY : eps8

   implicit none

   REAL(kind=DP) :: time!imaginary time t, just a check
   TYPE(polaw)   :: pr!polarization P(it) to be created
   TYPE(green)   ::  gf_p!green function G(it)
   TYPE(green)   ::  gf_m!green function G(-it)
   TYPE(q_mat)   ::  qm!overlap matrix Q
   LOGICAL       :: debug!if true check for hermeticity

   INTEGER :: iw,jw, ip,jp
   INTEGER :: l,n,m,o

!first annihilation
   call free_memory_polaw(pr)

!check time
   if((abs(gf_p%time - time)>=eps8) .OR. (abs(gf_m%time + time) >= eps8)) then!times are wrong
      write(stdout,*) 'Subroutine polarization: times are wrong',gf_p%time,gf_m%time
      stop
   endif
    
!set pr
   pr%ontime=.true.
   pr%time=time
   pr%numpw=qm%numpw

!allocate
   allocate(pr%pw( pr%numpw,pr%numpw))
   pr%pw(:,:) =(0.d0,0.d0)

   do iw=1,pr%numpw
      do jw=iw,pr%numpw
         do ip=1,qm%wp(iw)%numij
            do jp=1,qm%wp(jw)%numij

               l=qm%wp(iw)%ij(1,ip)
               n=qm%wp(iw)%ij(2,ip)
               m=qm%wp(jw)%ij(1,jp)
               o=qm%wp(jw)%ij(2,jp)

               pr%pw(iw,jw)=pr%pw(iw,jw)+qm%wp(iw)%o(ip)*qm%wp(jw)%o(jp)* &
                  & gf_p%gf(l,m,1)*gf_m%gf(n,o,1)

!couples are NOT ordered
              if(l/=n) then
                pr%pw(iw,jw)=pr%pw(iw,jw)+qm%wp(iw)%o(ip)*qm%wp(jw)%o(jp)* &
                  & gf_p%gf(n,m,1)*gf_m%gf(l,o,1)
              endif

              if(m/=o) then
                pr%pw(iw,jw)=pr%pw(iw,jw)+qm%wp(iw)%o(ip)*qm%wp(jw)%o(jp)* &
                  & gf_p%gf(l,o,1)*gf_m%gf(n,m,1)
              endif

              if(l/=n .AND. m/=o) then
                pr%pw(iw,jw)=pr%pw(iw,jw)+qm%wp(iw)%o(ip)*qm%wp(jw)%o(jp)* &
                  & gf_p%gf(n,o,1)*gf_m%gf(l,m,1)
              endif
           enddo
           pr%pw(jw,iw)=pr%pw(iw,jw)
         enddo
      enddo
    enddo


    pr%factor=(0.d0,-1.d0)
!now spin factor
    pr%pw(:,:)=2.d0*pr%pw(:,:)
    
    return

  END subroutine

  SUBROUTINE calculate_w(vp,pp,ww,xc_together,l_symm_epsilon,l_head_epsilon,agz,head,l_divergence,inv_epsi, &
                    &l_wing_epsilon, awing, l_verbose)
!this subroutine calculates W=(1+vp)^{-1}v
!this is meaningful only on frequency domain
!lapack routines are used
    USE io_global,            ONLY : stdout
    USE basic_structures,     ONLY : v_pot, head_epsilon

    implicit none 
    TYPE(v_pot)  :: vp!coulomb potential
    TYPE(polaw)  :: pp!polarization on imaginary frequency, destroyed on exit
    TYPE(polaw)  :: ww!dressed interaction to be calculated
    LOGICAL :: xc_together!if true the entire W is taken, otherwise W-v
    LOGICAL :: l_symm_epsilon! if true uses the symmetrized form of the dielectric matrix
                             !for calculating W
    LOGICAL :: l_head_epsilon!if true add to the symmetrized form  of the dielectric matrix
                             !the head terms
    REAL(kind=DP), DIMENSION(:) :: agz!terms A_ij<\tilde{w^P_j}|G=0>
    REAL(kind=DP) :: head!term (G=0,G=0) of the symmetric dielectric matrix
    LOGICAL, INTENT(in) :: l_divergence!if true calculate the head of the inverse dielectric matrix
    REAL(kind=DP), INTENT(out) :: inv_epsi!head of the inverse dielectric matrix
    LOGICAL, INTENT(in) :: l_wing_epsilon!if true calculate the wings of the symmetrized dielectric matrix
    REAL(kind=DP), DIMENSION(:) :: awing!the terms A_ij wing_j
    LOGICAL, INTENT(in) :: l_verbose

    INTEGER iw,jw,kw
    REAL(kind=DP), ALLOCATABLE, DIMENSION(:,:) :: dtmp!temporary array
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
    INTEGER :: info
    REAL(kind=DP),ALLOCATABLE, DIMENSION(:) :: work
    INTEGER :: lwork
    REAL(kind=DP) sca
    REAL(kind=DP) :: workd

!deallocate if the case
    call free_memory_polaw(ww)


!check and set
   if(pp%ontime) then
      write(stdout,*) 'Routine calculate_w: frequencies required'
      stop
   endif
   if(pp%numpw /= vp%numpw) then
      write(stdout,*) 'Routine calculate_w: basis set does not correspond',pp%numpw,vp%numpw
      stop
   endif

   ww%ontime=.false.
   ww%time=pp%time
   ww%label=pp%label
   ww%numpw=pp%numpw
   allocate(ww%pw(ww%numpw,ww%numpw))

   allocate(dtmp(ww%numpw,ww%numpw))
   allocate(ipiv(ww%numpw))


   if(.not.l_symm_epsilon) then
!not symmetric case calculates -vP
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,-1.d0*dble(pp%factor),&
           & vp%vmat,ww%numpw,pp%pw,ww%numpw,0.d0,dtmp,ww%numpw)
   else
!symmetric case calculates -v^1/2 P v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,-1.d0*dble(pp%factor),&
           & vp%vmat,ww%numpw,pp%pw,ww%numpw,0.d0,dtmp,ww%numpw)
      pp%pw(:,:)=dtmp(:,:)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & pp%pw,ww%numpw,vp%vmat,ww%numpw,0.d0,dtmp,ww%numpw)
   endif

!if required add the head
   if(l_symm_epsilon .and.l_head_epsilon) then
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            dtmp(iw,jw)=dtmp(iw,jw)+agz(iw)*agz(jw)*head
         enddo
      enddo
   endif


!if required add the wings
   if(l_symm_epsilon .and.l_wing_epsilon) then
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            dtmp(iw,jw)=dtmp(iw,jw)+agz(iw)*awing(jw)+agz(jw)*awing(iw)
         enddo
      enddo
   endif

   do iw=1,ww%numpw
      dtmp(iw,iw)=dtmp(iw,iw)+1.d0
   enddo


!inverse zmat

   call dgetrf(ww%numpw,ww%numpw,dtmp,ww%numpw,ipiv,info)
   if(info /= 0) then
     write(stdout,*) 'Routine calculate_w: problem with dgetrf :', info
     stop
   endif

   call dgetri(ww%numpw,dtmp,ww%numpw,ipiv,workd,-1,info)
   if(l_verbose) write(stdout,*) 'Dimension', workd!ATTENZIONE
   allocate(work(int(workd)))
   call dgetri(ww%numpw,dtmp,ww%numpw,ipiv,work,int(workd),info)
   
   if(info /= 0) then
     write(stdout,*) 'Routine calculate_w: problem with zgetri :', info
     stop
   endif

   
   if(.not.xc_together) then
      do iw=1,ww%numpw
         dtmp(iw,iw)=dtmp(iw,iw)-1.d0
      enddo
   endif

!if required calculates the head (G=0,G=0) of \epsilon^-1

   if(l_divergence) then
      inv_epsi=0.d0
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            inv_epsi = inv_epsi+dtmp(iw,jw)*agz(iw)*agz(jw)
         enddo
      enddo
      do iw=1,ww%numpw
         dtmp(iw,iw)=dtmp(iw,iw)-inv_epsi
      enddo
   endif

   if(l_verbose) write(stdout,*) 'INV EPSI G=0,G=0', inv_epsi

   if(.not. l_symm_epsilon) then
!calculates (e-1 -1)v
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & dtmp,ww%numpw,vp%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)
   else
         !calculates v^1/2 (e-1-1)v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & vp%vmat,ww%numpw,dtmp,ww%numpw,0.d0,pp%pw,ww%numpw)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
              & pp%pw,ww%numpw,vp%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)
   endif


   ww%factor=(1.d0,0.d0)

!   if(.not.xc_together) then
!     do iw=1,ww%numpw
!       do jw=1,ww%numpw
!         ww%pw(iw,jw)=ww%pw(iw,jw)-vp%vmat(iw,jw)
!       enddo
!     enddo
!   endif


   deallocate(dtmp,ipiv,work)
   return
  END SUBROUTINE

  SUBROUTINE create_polarization_contraction(time,pr,cp,uu,l_hf_energies, ene_hf)
!this subroutine set the polarization in imaginary time
! as P(r,r',it)=G(r,r',it)*G(r,r',-it)
! for our basis
!uses contractions
!THERE IS ALSO A SPIN FACTOR 2
!if required uses HF energies

   USE io_global,            ONLY : stdout
   USE constants,            ONLY : eps8
   USE compact_product,     ONLY : contraction_pola
   USE basic_structures,     ONLY : wannier_u

   implicit none

   REAL(kind=DP) :: time!imaginary time t, just a check
   TYPE(polaw)   :: pr!polarization P(it) to be created
   TYPE(contraction_pola) :: cp!the contracted products descriptor
   TYPE(wannier_u) :: uu!for the KS energies
   LOGICAL, INTENT(in) :: l_hf_energies!if true uses HF energies
   REAL(kind=DP), INTENT(in) :: ene_hf(:)!HF energies

   INTEGER :: iw,jw, vv, cc
   INTEGER :: l,n,m,o
   REAL(kind=DP) :: offset
   REAL(kind=DP),ALLOCATABLE  :: expene(:)!to calculate the exponentials just once

!first annihilation
   call free_memory_polaw(pr)


!set pr
   pr%ontime=.true.
   pr%time=time
   pr%numpw=cp%numpw

!calculates energy offset

   if(.not.l_hf_energies) then
      if(cp%nums > cp%nums_occ) then
         offset=-(uu%ene(cp%nums_occ+1,1)+uu%ene(cp%nums_occ,1))/2.d0
      else
         offset=-uu%ene(cp%nums_occ,1)
      endif
   else
       if(cp%nums > cp%nums_occ) then
         offset=-(ene_hf(cp%nums_occ+1)+ene_hf(cp%nums_occ))/2.d0
      else
         offset=-ene_hf(cp%nums_occ)
      endif
   endif
!calcualte exponentials of ks energies

  allocate(expene(cp%nums))
  if(.not.l_hf_energies) then
     do vv=1,cp%nums_occ
        expene(vv)=exp((uu%ene(vv,1)+offset)*time)
     enddo
     do cc=cp%nums_occ+1,cp%nums
        expene(cc)=exp(-(uu%ene(cc,1)+offset)*time)
     enddo
  else
     do vv=1,cp%nums_occ
        expene(vv)=exp((ene_hf(vv)+offset)*time)
     enddo
     do cc=cp%nums_occ+1,cp%nums
        expene(cc)=exp(-(ene_hf(cc)+offset)*time)
     enddo
  endif

!allocate
   allocate(pr%pw( pr%numpw,pr%numpw))
   pr%pw(:,:)=(0.d0,0.d0)
   do iw=1,pr%numpw
     do jw=iw,pr%numpw
        do vv=1,cp%nums_occ
           do cc=1,cp%nums-cp%nums_occ
              pr%pw(iw,jw)=pr%pw(iw,jw) + cp%ou(iw,vv,cc)*conjg(cp%ou(jw,vv,cc))*&
                & expene(vv)*expene(cc+cp%nums_occ)
           enddo
        enddo
        pr%pw(jw,iw)=pr%pw(iw,jw)
      enddo
   enddo
   pr%factor=(0.d0,-1.d0)
!now spin factor
   pr%pw(:,:)=2.d0*pr%pw(:,:)

   deallocate(expene)
   return
  END SUBROUTINE

  SUBROUTINE invert_ortho_polaw(op,opi)
!this subroutine inverts the orthonormalization matrix
!acting of products of wanniers
    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : ortho_polaw, free_memory

    implicit none

    TYPE(ortho_polaw), INTENT(in) :: op !the descriptor of the orthonormalization matrix to be inverted
    TYPE(ortho_polaw), INTENT(out) :: opi !the descriptor of the orthonormalization matrix to be inverted

    INTEGER :: info,lwork
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
    REAL(kind=DP), ALLOCATABLE, DIMENSION(:) :: work

    lwork=op%numpw
    allocate(ipiv(op%numpw))
    allocate(work(lwork))

    call free_memory(opi)
    opi%numpw=op%numpw
    allocate(opi%on_mat( opi%numpw, opi%numpw))
    opi%on_mat(:,:)=op%on_mat(:,:)   

   call dgetrf(opi%numpw,opi%numpw,opi%on_mat,opi%numpw,ipiv,info)
   if(info /= 0) then
     write(stdout,*) 'Routine invert_ortho_polaw: problem with dgetrf :', info
     stop
   endif

   call dgetri(opi%numpw,opi%on_mat,opi%numpw,ipiv,work,lwork,info)
   if(info /= 0) then
     write(stdout,*) 'Routine invert_ortho_polaw: problem with dgetri :', info
     stop
   endif

   if(op%inverse) then
     opi%inverse=.false.
   else
      opi%inverse=.true.
   endif

   deallocate(ipiv,work)
   return
  END SUBROUTINE

  SUBROUTINE distribute_ortho_polaw(op,opd)
!this subroutine distributes the orthonormalization matrix
! among processors

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : ortho_polaw, free_memory
    USE mp_world,         ONLY : nproc,mpime
    
    implicit none

    TYPE(ortho_polaw), INTENT(in) :: op !the descriptor of the orthonormalization matrix to be distributed
    TYPE(ortho_polaw), INTENT(out) :: opd!distributed orthonormalization matrix

    INTEGER :: l_blk,nbegin,nend,ii

    call free_memory(opd)

    opd%numpw = op%numpw
    opd%inverse = op%inverse
    
     l_blk= op%numpw/nproc
     if(l_blk*nproc < op%numpw) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > op%numpw) nend = op%numpw
     
     allocate(opd%on_mat(op%numpw,l_blk))
     
     do ii=nbegin,nend
        opd%on_mat(:,ii-nbegin+1)=op%on_mat(:,ii)
     enddo

    return

  END SUBROUTINE distribute_ortho_polaw

    SUBROUTINE collect_ortho_polaw(op,opd)
!this subroutine collects the orthonormalization matrix
! among processors

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : ortho_polaw, free_memory
    USE mp_world,         ONLY : nproc,mpime,world_comm!group
    USE parallel_include

    implicit none

    TYPE(ortho_polaw), INTENT(out) :: op !the descriptor of the orthonormalization matrix to be distributed
    TYPE(ortho_polaw), INTENT(in) :: opd!distributed orthonormalization matrix

    INTEGER :: l_blk,nbegin,nend,ierr

    call free_memory(op)

    op%numpw = opd%numpw
    op%inverse = opd%inverse

     l_blk= op%numpw/nproc
     if(l_blk*nproc < op%numpw) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > op%numpw) nend = op%numpw

     allocate(op%on_mat(op%numpw,l_blk*nproc))

#if defined(__MPI)
      call MPI_ALLGATHER(opd%on_mat,l_blk*op%numpw,MPI_DOUBLE_PRECISION, op%on_mat, &
           &    l_blk*op%numpw, MPI_DOUBLE_PRECISION,world_comm, ierr)
#else
      op%on_mat(:,:)=opd%on_mat
#endif
     

     return
   END SUBROUTINE collect_ortho_polaw
     
  



  SUBROUTINE orthonormalize(op,pw)
!this subroutine rotates the pw data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i,i'}P_{i',j'}B_{j,j'}

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : ortho_polaw

    implicit none

    TYPE(polaw), INTENT(inout) :: pw!data 
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)

    

    if(op%numpw /= pw%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))

    call dgemm('N','N',op%numpw,op%numpw,op%numpw,1.d0,&
                 & op%on_mat,op%numpw,pw%pw,op%numpw,0.d0,mat,op%numpw)



    call dgemm('N','T',op%numpw,op%numpw,op%numpw,1.d0,&
         & mat,op%numpw,op%on_mat,op%numpw,0.d0,pw%pw,op%numpw)

     
    deallocate(mat)

    return
   
  END SUBROUTINE

  SUBROUTINE orthonormalize_inverse(op,pw)
!this subroutine rotates the pw data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i',i}P_{i',j'}B_{j',j}

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : ortho_polaw

    implicit none

    TYPE(polaw), INTENT(inout) :: pw!data
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)



    if(op%numpw /= pw%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))
 
    call dgemm('T','N',op%numpw,op%numpw,op%numpw,1.d0,&
                 & op%on_mat,op%numpw,pw%pw,op%numpw,0.d0,mat,op%numpw)



    call dgemm('N','N',op%numpw,op%numpw,op%numpw,1.d0,&
         & mat,op%numpw,op%on_mat,op%numpw,0.d0,pw%pw,op%numpw)



    deallocate(mat)
    return

  END SUBROUTINE

  SUBROUTINE orthonormalize_vpot_inverse(op,vp)
!this subroutine rotates the v_pot data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i',i}P_{i',j'}B_{j',j}

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, ortho_polaw


    implicit none

    TYPE(v_pot), INTENT(inout) :: vp!data
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)


    if(op%numpw /= vp%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))
!    mat(:,:)=0.d0
!    do iw=1,op%numpw
!       do jw=1,op%numpw
!           do kw=1,op%numpw
!              mat(iw,jw)=mat(iw,jw)+op%on_mat(kw,iw)*vp%vmat(kw,jw)
!           enddo
!       enddo
!    enddo

    call dgemm('T','N',op%numpw,op%numpw,op%numpw,1.d0,op%on_mat,op%numpw,vp%vmat,op%numpw,0.d0,mat,op%numpw)


!    vp%vmat(:,:)=0.d0
!    do iw=1,op%numpw
!       do kw=1,op%numpw
!          do jw=1,op%numpw
!              vp%vmat(iw,jw)=vp%vmat(iw,jw)+op%on_mat(kw,jw)*mat(iw,kw)
!           enddo
!       enddo
!    enddo

    call dgemm('N','N',op%numpw,op%numpw,op%numpw,1.d0,mat,op%numpw,op%on_mat,op%numpw,0.d0,vp%vmat,op%numpw)

    


    deallocate(mat)
    return

  END SUBROUTINE

 SUBROUTINE orthonormalize_vpot(op,vp)
!this subroutine rotates the v_pot data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i,i'}P_{i',j'}B_{j,j'}

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, ortho_polaw


    implicit none

    TYPE(v_pot), INTENT(inout) :: vp!data
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)


    if(op%numpw /= vp%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))

!    mat(:,:)=0.d0
!    do iw=1,op%numpw
!       do jw=1,op%numpw
!           do kw=1,op%numpw
!              mat(iw,jw)=mat(iw,jw)+op%on_mat(iw,kw)*vp%vmat(kw,jw)
!           enddo
!       enddo
!    enddo

    call dgemm('N','N',op%numpw,op%numpw,op%numpw,1.d0,op%on_mat,op%numpw,vp%vmat,op%numpw,0.d0,mat,op%numpw)

!    vp%vmat(:,:)=0.d0
!    do iw=1,op%numpw
!       do jw=1,op%numpw
!          do kw=1,op%numpw
!              vp%vmat(iw,jw)=vp%vmat(iw,jw)+op%on_mat(jw,kw)*mat(iw,kw)
!           enddo
!       enddo
!    enddo

    call dgemm('N','T',op%numpw,op%numpw,op%numpw,1.d0,mat,op%numpw,op%on_mat,op%numpw,0.d0,vp%vmat,op%numpw)


    deallocate(mat)
    return

  END SUBROUTINE orthonormalize_vpot

 SUBROUTINE orthonormalize_vpot_para(op,vp)
!this subroutine rotates the v_pot data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i,i'}P_{i',j'}B_{j,j'}
!parallel version

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, ortho_polaw
    USE mp_world,  ONLY : mpime, nproc, world_comm
    USE mp,         ONLY : mp_sum


    implicit none

    TYPE(v_pot), INTENT(inout) :: vp!data
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)


    if(op%numpw /= vp%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))
    mat(:,:)=0.d0
    do iw=1,op%numpw
       do jw=1,op%numpw
          if(mod(jw,nproc)==mpime) then
             do kw=1,op%numpw
                mat(iw,jw)=mat(iw,jw)+op%on_mat(iw,kw)*vp%vmat(kw,jw)
             enddo
          endif
       enddo
       call mp_sum(mat(iw,:),world_comm)
    enddo

    vp%vmat(:,:)=0.d0
    do iw=1,op%numpw
       do jw=1,op%numpw
          if(mod(jw,nproc)==mpime) then
             do kw=1,op%numpw
                vp%vmat(iw,jw)=vp%vmat(iw,jw)+op%on_mat(jw,kw)*mat(iw,kw)
             enddo
          endif
       enddo
       call mp_sum(vp%vmat(iw,:),world_comm)
    enddo




    deallocate(mat)
    return

  END SUBROUTINE orthonormalize_vpot_para




  SUBROUTINE invert_v_pot(vp,vpi)
!this subroutine inverts the coulombian matrix
!acting on space of orthonormal products of wanniers
    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, free_memory

    implicit none

    TYPE(v_pot), INTENT(in) :: vp !the descriptor of the coulombian matrix to be inverted
    TYPE(v_pot), INTENT(inout) :: vpi !the descriptor of the inverted coulombian matrix

    INTEGER :: info,lwork,i,j
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
    REAL(kind=DP), ALLOCATABLE, DIMENSION(:) :: work

    call free_memory(vpi)

    lwork=vp%numpw
    allocate(ipiv(vp%numpw))
    allocate(work(lwork))

    vpi%numpw=vp%numpw
    allocate(vpi%vmat( vpi%numpw, vpi%numpw))
!write(stdout,*) size(vpi%vmat,1), size(vpi%vmat,2), size(vp%vmat,1), size(vp%vmat,2)
!    vpi%vmat(:,:)=vp%vmat(:,:)  ! bug
    do j = 1, size(vpi%vmat,2)
      do i = 1, size(vpi%vmat,1)
        vpi%vmat(i,j)=vp%vmat(i,j)
      end do
    end do

   call dgetrf(vpi%numpw,vpi%numpw,vpi%vmat,vpi%numpw,ipiv,info)
   if(info /= 0) then
     write(stdout,*) 'Invert V: problem with dgetrf :', info
     stop
   endif

   call dgetri(vpi%numpw,vpi%vmat,vpi%numpw,ipiv,work,lwork,info)
   if(info /= 0) then
     write(stdout,*) 'Invert V: problem with dgetri :', info
     stop
   endif

   deallocate(ipiv,work)
   return
 END SUBROUTINE invert_v_pot

 SUBROUTINE fake_polarization_io(n)
!this subroutine just call mp_barrier 2*n+1 times
   USE mp,       ONLY : mp_barrier
   USE mp_world, ONLY : world_comm

   implicit none

   INTEGER :: n
   INTEGER :: i
!we take advantage of the t ==> -t symmetry
!   do i=-n,n
   do i=0,n
     call mp_barrier( world_comm )
   enddo
   return
 END SUBROUTINE fake_polarization_io
   
  SUBROUTINE orthonormalize_vpot_inverse_para(op,vp)
!this subroutine rotates the v_pot data on the basis of the trasform op
!perform the trasform \sum_{i',j'} B_{i',i}P_{i',j'}B_{j',j}
!parallel version

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, ortho_polaw
    USE mp_world,  ONLY : world_comm, mpime, nproc
    USE mp,         ONLY : mp_sum

    implicit none

    TYPE(v_pot), INTENT(inout) :: vp!data
    TYPE(ortho_polaw), INTENT(in) :: op!trasform

    INTEGER :: iw,jw,kw
    REAL(kind=DP), ALLOCATABLE :: mat(:,:)


    if(op%numpw /= vp%numpw) then
      write(stdout,*) 'ROUTINE ORTHONORMALIZE: BASIS INCONSISTENT'
      stop
    endif

    allocate(mat(op%numpw,op%numpw))
    mat(:,:)=0.d0
    do iw=1,op%numpw
       do jw=1,op%numpw
          if(mod(jw,nproc)==mpime) then
             do kw=1,op%numpw
                mat(iw,jw)=mat(iw,jw)+op%on_mat(kw,iw)*vp%vmat(kw,jw)
             enddo
          endif
       enddo
       call mp_sum(mat(iw,:),world_comm)
    enddo

    vp%vmat(:,:)=0.d0
    do iw=1,op%numpw
       do jw=1,op%numpw
          if(mod(jw,nproc)==mpime) then
             do kw=1,op%numpw
                vp%vmat(iw,jw)=vp%vmat(iw,jw)+op%on_mat(kw,jw)*mat(iw,kw)
             enddo
          endif
       enddo
       call mp_sum(vp%vmat(iw,:),world_comm)
    enddo




    deallocate(mat)
    return

  END SUBROUTINE orthonormalize_vpot_inverse_para

  SUBROUTINE create_polarization_contraction_state(time,pr,uu,l_hf_energies, ene_hf,options)
!this subroutine set the polarization in imaginary time
! as P(r,r',it)=G(r,r',it)*G(r,r',-it)
! for our basis
!uses contractions
!THERE IS ALSO A SPIN FACTOR 2
!if required uses HF energies
!use single occupied state contractions

   USE io_global,            ONLY : stdout
   USE constants,            ONLY : eps8
   USE compact_product,      ONLY : contraction_pola_state, free_memory_contraction_pola_state, &
                                               &read_contraction_pola_state
   USE basic_structures,     ONLY : wannier_u
   USE input_gw,             ONLY : input_options

   implicit none

   REAL(kind=DP) :: time!imaginary time t, just a check
   TYPE(polaw)   :: pr!polarization P(it) to be created
   TYPE(wannier_u) :: uu!for the KS energies
   LOGICAL, INTENT(in) :: l_hf_energies!if true uses HF energies
   REAL(kind=DP), INTENT(in) :: ene_hf(:)!HF energies
   TYPE(input_options) :: options!for i/o purpoose

   INTEGER :: iw,jw, vv, cc
   INTEGER :: l,n,m,o
   REAL(kind=DP) :: offset
   REAL(kind=DP),ALLOCATABLE  :: expene(:)!to calculate the exponentials just once
   REAL(kind=DP),ALLOCATABLE  :: exptmp(:)!temporary arrays for speed-up
   REAL(kind=DP),ALLOCATABLE :: outmp(:,:)
   REAL(kind=DP),ALLOCATABLE :: tmpc(:)
   TYPE(contraction_pola_state) :: cps

!first annihilation
   call free_memory_polaw(pr)


!set pr
   pr%ontime=.true.
   pr%time=time

!the following for accessing dimension of polarization matrix
   cps%state=1
   call read_contraction_pola_state(cps,options)
   pr%numpw=cps%numpw

!calculates energy offset

   if(.not.l_hf_energies) then
      if(uu%nums > uu%nums_occ(1)) then
         offset=-(uu%ene(uu%nums_occ(1)+1,1)+uu%ene(uu%nums_occ(1),1))/2.d0
      else
         offset=-uu%ene(uu%nums_occ(1),1)
      endif
   else
       if(uu%nums > uu%nums_occ(1)) then
         offset=-(ene_hf(uu%nums_occ(1)+1)+ene_hf(uu%nums_occ(1)))/2.d0
      else
         offset=-ene_hf(uu%nums_occ(1))
      endif
   endif
!calcualte exponentials of ks energies

  allocate(expene(uu%nums))

  if(.not.l_hf_energies) then
     do vv=1,uu%nums_occ(1)
        expene(vv)=exp((uu%ene(vv,1)+offset)*time)
     enddo
     do cc=uu%nums_occ(1)+1,uu%nums
        expene(cc)=exp(-(uu%ene(cc,1)+offset)*time)
     enddo
  else
     do vv=1,uu%nums_occ(1)
        expene(vv)=exp((ene_hf(vv)+offset)*time)
     enddo
     do cc=uu%nums_occ(1)+1,uu%nums
        expene(cc)=exp(-(ene_hf(cc)+offset)*time)
     enddo
  endif

!allocate
   allocate(pr%pw( pr%numpw,pr%numpw))
   pr%pw(:,:)=0.d0

   allocate(exptmp(cps%nums-cps%nums_occ))
   allocate(outmp(cps%nums-cps%nums_occ,pr%numpw))
   allocate(tmpc(cps%nums-cps%nums_occ))
  do vv=1,cps%nums_occ
      cps%state=vv
      write(stdout,*) 'read_contraction_pola_state'
      if(vv /= 1) call read_contraction_pola_state(cps,options)
      write(stdout,*) 'read_contraction_pola_state', vv
      do cc=1,cps%nums-cps%nums_occ
         exptmp(cc)=expene(vv)*expene(cc+cps%nums_occ)
      enddo
      do iw=1,pr%numpw
         do cc=1,cps%nums-cps%nums_occ
            outmp(cc,iw)=cps%ou(cc,iw)*exptmp(cc)
         enddo
      enddo
      write(stdout,*) 'calculus1'
!      do jw=1,pr%numpw
!         do iw=jw,pr%numpw
!            do cc=1,cps%nums-cps%nums_occ
!               tmpc(cc)=cps%ou(cc,iw)*outmp(cc,jw)
!            enddo
!            pr%pw(iw,jw)=pr%pw(iw,jw)+sum(tmpc(:))
!         enddo
!      enddo

      call dgemm('T','N',pr%numpw,pr%numpw,cps%nums-cps%nums_occ,1.d0,cps%ou,cps%nums-cps%nums_occ, &
                    outmp,cps%nums-cps%nums_occ,1.d0,pr%pw,pr%numpw)

      call free_memory_contraction_pola_state(cps)
   enddo
   do jw=1,pr%numpw 
      do iw=jw,pr%numpw   
         pr%pw(jw,iw)=pr%pw(iw,jw) 
      enddo
   enddo
!   pr%pw(:,:)=(0.d0,-1.d0)*pr%pw(:,:)
   pr%factor=(0.d0,-1.d0)
!now spin factor
   pr%pw(:,:)=2.d0*pr%pw(:,:)

   deallocate(expene)
   deallocate(exptmp)
   deallocate(outmp)
   deallocate(tmpc)
   return
 END SUBROUTINE create_polarization_contraction_state


 
 SUBROUTINE distribute_v_pot(vp,vpd)
!this subroutine distributes the coulomb matrix
! among processors

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, free_memory
    USE mp_world,         ONLY : nproc,mpime

    implicit none

    TYPE(v_pot), INTENT(in) :: vp !the descriptor of the coulomb matrix to be distributed
    TYPE(v_pot), INTENT(out) :: vpd!distributed coulomb matrix

    INTEGER :: l_blk,nbegin,nend,ii

    call free_memory(vpd)

    vpd%numpw = vp%numpw
 

     l_blk= vp%numpw/nproc
     if(l_blk*nproc < vp%numpw) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > vp%numpw) nend = vp%numpw

     allocate(vpd%vmat(vp%numpw,l_blk))

     do ii=nbegin,nend
        vpd%vmat(:,ii-nbegin+1)=vp%vmat(:,ii)
     enddo

    return

  END SUBROUTINE distribute_v_pot

    SUBROUTINE collect_v_pot(vp,vpd)
!this subroutine collects the coulomb matrix
! among processors

    USE io_global, ONLY : stdout
    USE basic_structures, ONLY : v_pot, free_memory
    USE mp_world,         ONLY : nproc,mpime,world_comm!group
    USE parallel_include

    implicit none

    TYPE(v_pot), INTENT(out) :: vp !the descriptor of the coulomb matrix to be distributed
    TYPE(v_pot), INTENT(in) :: vpd!distributed coulomb matrix

    INTEGER :: l_blk,nbegin,nend,ierr

    call free_memory(vp)

    vp%numpw = vpd%numpw

     l_blk= vp%numpw/nproc
     if(l_blk*nproc < vp%numpw) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > vp%numpw) nend = vp%numpw

     allocate(vp%vmat(vp%numpw,l_blk*nproc))

#if defined(__MPI)
      call MPI_ALLGATHER(vpd%vmat,l_blk*vp%numpw,MPI_DOUBLE_PRECISION, vp%vmat, &
           &    l_blk*vp%numpw, MPI_DOUBLE_PRECISION,world_comm, ierr)
#else
      vp%vmat(:,:)=vpd%vmat(:,:)
#endif



     return
   END SUBROUTINE collect_v_pot

   SUBROUTINE calculate_w_g(vp,pp,ww,xc_together,l_symm_epsilon,l_head_epsilon,agz,head,l_divergence,inv_epsi, &
                    &l_wing_epsilon, awing, awing_c)
!this subroutine calculates W=(1+vp)^{-1}v
!this is meaningful only on frequency domain
!lapack routines are used
!it use exteded basis for G=0 term


    USE io_global,            ONLY : stdout
    USE basic_structures,     ONLY : v_pot, head_epsilon

    implicit none
    TYPE(v_pot)  :: vp!coulomb potential
    TYPE(polaw)  :: pp!polarization on imaginary frequency, destroyed on exit
    TYPE(polaw)  :: ww!dressed interaction to be calculated
    LOGICAL :: xc_together!if true the entire W is taken, otherwise W-v
    LOGICAL :: l_symm_epsilon! if true uses the symmetrized form of the dielectric matrix
                             !for calculating W
    LOGICAL :: l_head_epsilon!if true add to the symmetrized form  of the dielectric matrix
                             !the head terms
    REAL(kind=DP), DIMENSION(:) :: agz!terms A_ij<\tilde{w^P_j}|G=0>
    REAL(kind=DP) :: head!term (G=0,G=0) of the symmetric dielectric matrix
    LOGICAL, INTENT(in) :: l_divergence!if true calculate the head of the inverse dielectric matrix
    REAL(kind=DP), INTENT(out) :: inv_epsi!head of the inverse dielectric matrix
    LOGICAL, INTENT(in) :: l_wing_epsilon!if true calculate the wings of the symmetrized dielectric matrix
    REAL(kind=DP), DIMENSION(:) :: awing!the terms A_ij wing_j
    REAL(kind=DP), DIMENSION(:) :: awing_c!the terms A_ij wing_j


    INTEGER iw,jw,kw
    REAL(kind=DP), ALLOCATABLE, DIMENSION(:,:) :: dtmp!temporary array
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
    INTEGER :: info
    REAL(kind=DP),ALLOCATABLE, DIMENSION(:) :: work
    INTEGER :: lwork
    REAL(kind=DP) sca
    REAL(kind=DP) :: workd

    REAL(kind=DP) alpha

!deallocate if the case
    call free_memory_polaw(ww)


!check and set
   if(pp%ontime) then
      write(stdout,*) 'Routine calculate_w: frequencies required'
      stop
   endif
   if(pp%numpw /= vp%numpw) then
      write(stdout,*) 'Routine calculate_w: basis set does not correspond',pp%numpw,vp%numpw
      stop
   endif

   ww%ontime=.false.
   ww%time=pp%time
   ww%label=pp%label
   ww%numpw=pp%numpw
   allocate(ww%pw(ww%numpw,ww%numpw))

   allocate(dtmp(ww%numpw+1,ww%numpw+1))
   allocate(ipiv(ww%numpw+1))

   dtmp(:,:)=0.d0

   if(.not.l_symm_epsilon) then
!not symmetric case calculates -vP
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,-1.d0*dble(pp%factor),&
           & vp%vmat,ww%numpw,pp%pw,ww%numpw,0.d0,dtmp,ww%numpw+1)
   else
!symmetric case calculates -v^1/2 P v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,-1.d0*dble(pp%factor),&
           & vp%vmat,ww%numpw,pp%pw,ww%numpw,0.d0,dtmp,ww%numpw+1)
      pp%pw(1:ww%numpw,1:ww%numpw)=dtmp(1:ww%numpw,1:ww%numpw)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & pp%pw,ww%numpw,vp%vmat,ww%numpw,0.d0,dtmp,ww%numpw+1)
   endif

!calculate normalization factor alpha

   sca=1.d0
   do iw=1,ww%numpw
      sca=sca-agz(iw)**2.d0
   enddo
   alpha=1.d0/sqrt(sca)
   
   write(stdout,*) 'ALPHA :', alpha
   FLUSH(stdout)

!calculate elements 0',0'   0',i   i,O'
   dtmp(ww%numpw+1,:)=0.d0
   dtmp(:,ww%numpw+1)=0.d0

   do iw=1,ww%numpw
      do jw=1,ww%numpw
         dtmp(ww%numpw+1,ww%numpw+1)=dtmp(ww%numpw+1,ww%numpw+1) &
                       & + alpha**2.d0 * agz(iw)* agz(jw)*dtmp(iw,jw)
      enddo
   enddo

   do iw=1,ww%numpw
      do jw=1,ww%numpw
         dtmp(ww%numpw+1,iw)=dtmp(ww%numpw+1,iw)-alpha*dtmp(jw,iw)*agz(jw)
      enddo
   enddo
   do iw=1,ww%numpw
      dtmp(iw,ww%numpw+1)=dtmp(ww%numpw+1,iw)
   enddo

!if required add the head
   if(l_symm_epsilon .and.l_head_epsilon) then

! terms i,j
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            dtmp(iw,jw)=dtmp(iw,jw)+agz(iw)*agz(jw)*head
         enddo
      enddo
   endif
!term O',0'  0',i  i,0'
   dtmp(ww%numpw+1,ww%numpw+1)= dtmp(ww%numpw+1,ww%numpw+1) + &
                        &head/alpha**2.d0
    do iw=1,ww%numpw
       dtmp(ww%numpw+1,iw)=dtmp(ww%numpw+1,iw)+head*agz(iw)/alpha
       dtmp(iw,ww%numpw+1)=dtmp(iw,ww%numpw+1)+head*agz(iw)/alpha
    enddo
   


!if required add the wings
!TODO ATTENZIONE
   if(l_symm_epsilon .and.l_wing_epsilon) then
!elements i,j
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            dtmp(iw,jw)=dtmp(iw,jw)+agz(iw)*awing_c(jw)+agz(jw)*awing_c(iw)
         enddo
      enddo
! 0',0'
      do iw=1,ww%numpw
         dtmp(ww%numpw+1,ww%numpw+1)=dtmp(ww%numpw+1,ww%numpw+1)  &
                 & - 2.d0*agz(iw)*awing_c(iw)
      enddo
!  0',i   ',0'
      do iw=1,ww%numpw
         dtmp(ww%numpw+1,iw)=dtmp(ww%numpw+1,iw)+awing(iw)/alpha
         do jw=1,ww%numpw
            dtmp(ww%numpw+1,iw)=dtmp(ww%numpw+1,iw)-alpha*agz(iw)*awing_c(jw)*agz(jw)
         enddo
         dtmp(iw,ww%numpw+1)=dtmp(ww%numpw+1,iw)
      enddo

   endif

   do iw=1,ww%numpw+1
      dtmp(iw,iw)=dtmp(iw,iw)+1.d0
   enddo


!inverse zmat

   write(stdout,*) 'Before inversion'
   FLUSH(stdout)
   call dgetrf(ww%numpw+1,ww%numpw+1,dtmp,ww%numpw+1,ipiv,info)
   if(info /= 0) then
     write(stdout,*) 'Routine calculate_w: problem with dgetrf :', info
     stop
   endif
   write(stdout,*) 'Before inversion2'
   FLUSH(stdout)

   call dgetri(ww%numpw+1,dtmp,ww%numpw+1,ipiv,workd,-1,info)
   write(stdout,*) 'Dimension', workd,ww%numpw,info!ATTENZIONE
   FLUSH(stdout)
   allocate(work(int(workd)))
   call dgetri(ww%numpw+1,dtmp,ww%numpw+1,ipiv,work,int(workd),info)

   write(stdout,*) 'Out of dgetri'
   FLUSH(stdout)


   if(info /= 0) then
     write(stdout,*) 'Routine calculate_w: problem with zgetri :', info
     stop
   endif


   if(.not.xc_together) then
      do iw=1,ww%numpw+1
         dtmp(iw,iw)=dtmp(iw,iw)-1.d0
      enddo
   endif

!if required calculates the head (G=0,G=0) of \epsilon^-1

   if(l_divergence) then
      inv_epsi=0.d0

!term i,j
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            inv_epsi = inv_epsi+dtmp(iw,jw)*agz(iw)*agz(jw)
         enddo
      enddo
!term 0',0'
      inv_epsi=inv_epsi+dtmp(ww%numpw+1,ww%numpw+1)/alpha**2.d0

!terms 0',i i,0'
      do iw=1,ww%numpw
         inv_epsi=inv_epsi+2.d0*agz(iw)*dtmp(ww%numpw+1,iw)/alpha
      enddo


      do iw=1,ww%numpw+1
         dtmp(iw,iw)=dtmp(iw,iw)-inv_epsi
      enddo
   endif

   write(stdout,*) 'INV EPSI G=0,G=0', inv_epsi, dtmp(ww%numpw+1,ww%numpw+1)
   FLUSH(stdout)
   if(l_symm_epsilon .and.l_head_epsilon) then!ATTENZIONE
!take away the G=0,G=0 term
! terms i,j
      write(stdout,*) 'Extract G=0,G=0 term'
      do jw=1,ww%numpw
         do iw=1,ww%numpw
            dtmp(iw,jw)=dtmp(iw,jw)-agz(iw)*agz(jw)*inv_epsi
         enddo
      enddo
   endif




   if(.not. l_symm_epsilon) then
!calculates (e-1 -1)v
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & dtmp,ww%numpw+1,vp%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)
   else
         !calculates v^1/2 (e-1-1)v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & vp%vmat,ww%numpw,dtmp,ww%numpw+1,0.d0,pp%pw,ww%numpw)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
              & pp%pw,ww%numpw,vp%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)
   endif


   ww%factor=(1.d0,0.d0)

!   if(.not.xc_together) then
!     do iw=1,ww%numpw
!       do jw=1,ww%numpw
!         ww%pw(iw,jw)=ww%pw(iw,jw)-vp%vmat(iw,jw)
!       enddo
!     enddo
!   endif


   deallocate(dtmp,ipiv,work)
   return
 END SUBROUTINE calculate_w_g



 SUBROUTINE create_polarization_file(uu, tf, prefix)
   
   USE basic_structures, ONLY : wannier_u, cprim_prod, free_memory, &
                              &initialize_memory_cprim_prod
   USE times_gw,  ONLY : times_freqs
   USE mp_world,   ONLY : world_comm, nproc,mpime
   USE io_global,  ONLY : stdout
   USE mp, ONLY : mp_barrier

   implicit none

   TYPE(wannier_u), INTENT(in) :: uu!for the energies
   TYPE(times_freqs) :: tf
   CHARACTER(LEN=256), INTENT(in) ::  prefix!to designate the PW files

   INTEGER :: l_blk, nbegin, nend
   INTEGER :: it,iv,ic, iw, jw
   TYPE(polaw) :: pp
   TYPE(cprim_prod) :: cpp
   LOGICAL :: ok_read
   REAL(kind=DP), ALLOCATABLE :: exp_table(:)
   LOGICAL :: l_first!trick for gaining access to numpw, quite BAD
   REAL(kind=DP), ALLOCATABLE :: cpmat_tmp(:,:)

   write(stdout,*) 'Routine create_polarization_file'

   allocate(exp_table(uu%nums-uu%nums_occ(1)))

!loop on time

    l_blk= (tf%n+1)/nproc
    if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1
    nbegin=mpime*l_blk
    nend=nbegin+l_blk-1
    
    do it=nbegin,nend


       if(it<= tf%n) then
          call initialize_polaw(pp)
          l_first=.true.
          do iv=1,uu%nums_occ(1)
 !loop on v
             write(stdout,*) 'STATE', iv
             FLUSH(stdout)

!set table of exponentials

             do ic=uu%nums_occ(1)+1,uu%nums
                exp_table(ic-uu%nums_occ(1))=exp((uu%ene(iv,1)-uu%ene(ic,1))*tf%times(it))
             enddo
             call mp_barrier( world_comm )
             call initialize_memory_cprim_prod(cpp)

   !!read in Svci 
             cpp%cprim=iv
              call read_data_pw_cprim_prod(cpp, prefix, .true., ok_read, .true.,.false.)
!if required allocate polaw

              if(l_first) then
                 allocate(pp%pw(cpp%numpw,cpp%numpw))
                 pp%pw(:,:)=0.d0
                 pp%numpw=cpp%numpw
                 l_first=.false.
              endif
!!!!the following for using blas routine
              allocate(cpmat_tmp(cpp%nums_cond,cpp%numpw))
              call mytranspose(cpp%cpmat, cpp%numpw, cpmat_tmp, cpp%nums_cond, cpp%numpw, cpp%nums_cond)
              do iw=1,cpp%numpw
                 cpmat_tmp(:,iw)=cpmat_tmp(:,iw)*exp_table(:)
              enddo

!!!
    !calculate term
!              do ic=1,cpp%nums_cond
!                 do jw=1,cpp%numpw
!                    do iw=1,cpp%numpw
!                       pp%pw(iw,jw)=pp%pw(iw,jw)+cpp%cpmat(iw,ic)*cpp%cpmat(jw,ic)*exp_table(ic)
!                    enddo
!                 enddo
!              enddo


              call dgemm('N','N',cpp%numpw,cpp%numpw,cpp%nums_cond,1.d0,cpp%cpmat,cpp%numpw,&
                   &cpmat_tmp,cpp%nums_cond,1.d0,pp%pw,cpp%numpw)

              deallocate(cpmat_tmp)
             call free_memory(cpp)
          enddo
 !write polaw on file
!
          pp%label=it
          pp%time=tf%times(it)
          pp%factor=(0.d0,-1.d0)
          pp%numpw=cpp%numpw
          pp%ontime=.true.
          pp%pw(:,:)=pp%pw(:,:)*2.d0!for spin multiplicity
          call write_polaw(pp,.false.)
          call free_memory_polaw(pp)
          else
!just parallelel reading of file
             do iv=1,uu%nums_occ(1)
                call mp_barrier( world_comm )
                call initialize_memory_cprim_prod(cpp)
                cpp%cprim=iv
                call read_data_pw_cprim_prod(cpp, prefix, .true., ok_read, .true.,.false.)
                call free_memory(cpp)
             enddo
          endif
       enddo

       deallocate(exp_table)
   return

 END SUBROUTINE create_polarization_file


 SUBROUTINE square_root_polaw(pw,numpw)
!this subroutine calculate the square root of the polaw matrix
!it is done by calculating the eigenvalues and eigenvectors
!it assumes that the matrix is symmetric and positive definite

   USE io_global, ONLY : stdout

   implicit none

   REAL(kind=DP) :: pw(numpw,numpw)!the matrix to be operated
   INTEGER :: numpw!dimension of the matrix

   REAL(kind=DP), ALLOCATABLE :: e_va(:), work(:)
   REAL(kind=DP), ALLOCATABLE :: tmp_pw(:,:)

   INTEGER :: lwork, info
   REAL(kind=DP) :: loptwork
   INTEGER :: iw,jw,kw

#if defined(__OPENMP)
   INTEGER :: ntids
   INTEGER :: omp_get_num_threads, omp_get_max_threads
   EXTERNAL omp_set_num_threads, omp_get_num_threads, omp_get_max_threads
#endif


   allocate(e_va(numpw))
   allocate(tmp_pw(numpw, numpw))

   tmp_pw(:,:)=pw(:,:)

#if defined(__OPENMP)
     ! go single-thread
     ntids = omp_get_max_threads()
     call omp_set_num_threads(1)
#endif


!check for optimal dimension
   call DSYEV( 'V', 'U', numpw, tmp_pw, numpw, e_va, loptwork, -1, info)

   lwork=int(loptwork)
   allocate(work(lwork))

!calculate the eigenvalues, eigenvectors

   call DSYEV( 'V', 'U', numpw, tmp_pw, numpw, e_va, work, lwork,info )
   if(info /= 0) then
      write(stdout,*) 'Problem with dsyev', info
      stop
   endif

#if defined(__OPENMP)
     ! go multi-thread
     call omp_set_num_threads(ntids)
#endif

!do square root of eigenvector
   do iw=1,numpw
      if(e_va(iw) < 0.d0) then
         write(stdout,*) 'Problem with eigenvalue', iw
         stop
      endif
      e_va(iw)=dsqrt(e_va(iw))
   enddo
   
!reform the matrix

   pw(:,:)=0.d0

! Carlo substitute with DGEMM
   do kw=1,numpw  
      do jw=1,numpw
         do iw=1,numpw
            pw(iw,jw)=pw(iw,jw)+tmp_pw(iw,kw)*tmp_pw(jw,kw)*e_va(kw)
         enddo
      enddo
   enddo
      
   deallocate(tmp_pw)
   deallocate(work)
   deallocate(e_va)
   return

 END SUBROUTINE square_root_polaw


  SUBROUTINE create_polarization_beta(time, pr, uu, qm)

!this subroutine create the polarization with the strategy beta:
!P_ij(\tau)=\int dr dr' <\omega^{P'}_i(r)U_{v,v'}exp(E_v\tau)\tilde{\omega_v'}(r)*U_{v,v''}\tilde{\omega_v''}(r')
!         *U_{c,c'}exp(E_c\tau)\tilde{\omega_c'}(r)U_{c,c''}\tilde{\omega_c''}(r)\omega^{P'}_j(r')
!it makes use of S_{i,vc}=\int dr \omega^{P'}_i(r)\tilde{\omega_v}(r)\tilde{\omega_c}(r)


   USE basic_structures, ONLY : wannier_u, free_memory,q_mat, wannier_P
   USE times_gw,  ONLY : times_freqs
   USE io_global,  ONLY : stdout

   implicit none

   REAL(kind=DP), INTENT(in) :: time! imaginary time tau
   TYPE(wannier_u), INTENT(in) :: uu!for the energies and trasformation matrix               
   TYPE(polaw), INTENT(out)   :: pr!polarization P(it) to be created
   TYPE(q_mat), INTENT(in) :: qm ! for S matrices
                      

   INTEGER :: i,j,k, ip, ii, jj 
   INTEGER :: nums_cond!number of conduction states
   INTEGER :: iw,jw
   REAL(kind=DP), ALLOCATABLE :: v_val(:,:), exp_table_v(:), v_cond(:,:), exp_table_c(:)
   REAL(kind=DP), ALLOCATABLE :: tmp_mat1(:,:), tmp_mat2(:,:)
   REAL(kind=DP) :: fermi_en
   REAL(kind=DP), ALLOCATABLE :: q(:,:),t(:,:), v(:)
   REAL(kind=DP), EXTERNAL :: ddot

   write(stdout,*) 'Routine create_polarization_beta'
!0)set polarization structure

   pr%ontime=.true.
   pr%time=time
   pr%numpw=qm%numpw
   pr%factor=(0.d0,-1.d0)
!allocate                                                                                                                  
   allocate(pr%pw( pr%numpw,pr%numpw))
   pr%pw(:,:) =(0.d0,0.d0)



!1)calculate V_v'v''= U_{vv'}U_{v,v''}exp(E_v \tau}
   allocate(v_val(uu%nums_occ(1),uu%nums_occ(1)))
   allocate(exp_table_v(uu%nums_occ(1)))

!fermi_en is used to reduce the numerical error
   fermi_en=(uu%ene(uu%nums_occ(1)+1,1)+uu%ene(uu%nums_occ(1),1))/2.d0
   exp_table_v(1:uu%nums_occ(1))=exp((uu%ene(1:uu%nums_occ(1),1)-fermi_en)*abs(time))

   v_val(:,:)=0.d0

   allocate(tmp_mat1(uu%nums_occ(1),uu%nums_occ(1)), tmp_mat2(uu%nums_occ(1),uu%nums_occ(1)))

   tmp_mat1(1:uu%nums_occ(1),1:uu%nums_occ(1))=dble(uu%umat(1:uu%nums_occ(1),1:uu%nums_occ(1),1))
   do i=1,uu%nums_occ(1)
      do j=1,uu%nums_occ(1)
         tmp_mat2(i,j)=dble(uu%umat(i,j,1))*exp_table_v(i)
      enddo
   enddo

   call dgemm('T','N',uu%nums_occ(1),uu%nums_occ(1),uu%nums_occ(1),1.d0,tmp_mat2,uu%nums_occ(1),tmp_mat1,uu%nums_occ(1),&
        &0.d0,v_val,uu%nums_occ(1))
   deallocate(tmp_mat1,tmp_mat2)
   deallocate(exp_table_v)

!2) calculate V_c'c''= U_{c,c'}U_{c,c''}exp(-E_c \tau}
   nums_cond=uu%nums-uu%nums_occ(1)

   allocate(v_cond(nums_cond,nums_cond))
   allocate(exp_table_c(nums_cond))

   exp_table_c(1:nums_cond)=exp((-uu%ene(uu%nums_occ(1)+1:uu%nums,1) +fermi_en)*abs(time))

   allocate(tmp_mat1(nums_cond,nums_cond), tmp_mat2(nums_cond,nums_cond))

   tmp_mat1(1:nums_cond,1:nums_cond)=dble(uu%umat(uu%nums_occ(1)+1:uu%nums,uu%nums_occ(1)+1:uu%nums,1))
   do i=1,nums_cond
      do j=1,nums_cond
         tmp_mat2(i,j)=dble(uu%umat(uu%nums_occ(1)+i,uu%nums_occ(1)+j,1))*exp_table_c(i)
      enddo
   enddo

   call dgemm('T','N',nums_cond,nums_cond,nums_cond,1.d0,tmp_mat2,nums_cond,tmp_mat1,nums_cond,&
        &0.d0,v_cond,nums_cond)
   deallocate(tmp_mat1,tmp_mat2)


   deallocate(exp_table_c)

   do iw=1,pr%numpw

!calculate T_v''c'=S_{i,v'c'}V_{v',v''}
      allocate(t(uu%nums_occ(1),nums_cond))
      t(:,:)=0.d0
      do ip=1,qm%wp(iw)%numij
         i=qm%wp(iw)%ij(1,ip)!valence only
         j=qm%wp(iw)%ij(2,ip)!valence and conduction
         if(i>uu%nums_occ(1)) then
            write(stdout,*) 'create_polarization_beta ERROR'
            FLUSH(stdout)
            stop
         endif
!only valence*conduction products are required
         if(j>uu%nums_occ(1)) then
            do ii=1,uu%nums_occ(1)
               t(ii,j-uu%nums_occ(1))=t(ii,j-uu%nums_occ(1))+qm%wp(iw)%o(ip)*v_val(i,ii)
            enddo
         endif
      enddo
   
!calculate Q v''c''=T_{v''c'}V_{c'c''}
      allocate( q(uu%nums_occ(1),nums_cond))
      call dgemm('N','N',uu%nums_occ(1),nums_cond,nums_cond,1.d0,t,uu%nums_occ(1),v_cond,nums_cond,0.d0,&
           &q,uu%nums_occ(1))
      deallocate(t)
   
!put q on a right order for multiplication with S_{j,v''c''}
!WARNING WARNING ATTENZIONE
!it suppose that the wp(:)%ij are all the same
      allocate(v(qm%wp(1)%numij))
      v(:)=0.d0
      do ip=1,qm%wp(iw)%numij
         i=qm%wp(iw)%ij(1,ip)!valence only                                                                        
         j=qm%wp(iw)%ij(2,ip)!valence and conduction
         if(j > uu%nums_occ(1)) then
            v(ip)=q(i,j-uu%nums_occ(1))
         endif
      enddo
    
      deallocate(q)
!product with jw
      do jw=iw,pr%numpw
         pr%pw(iw,jw)= ddot(qm%wp(iw)%numij,qm%wp(jw)%o,1,v,1)
         pr%pw(jw,iw)=pr%pw(iw,jw)
      enddo

      deallocate(v)
   enddo
!now spin factor 
                                                                                                             
    pr%pw(:,:)=2.d0*pr%pw(:,:)

    deallocate(v_val,v_cond)
   return

 END SUBROUTINE create_polarization_beta


 SUBROUTINE create_polarization_upper(uu, tf, prefix)
!this subroutine adds to the polarization the part from upper reduced states   

   USE basic_structures, ONLY : wannier_u, cprim_prod, free_memory, &
                              &initialize_memory, upper_states
   USE times_gw,   ONLY : times_freqs
   USE mp_world,   ONLY : nproc,mpime
   USE io_global,  ONLY : stdout

   implicit none

   TYPE(wannier_u), INTENT(in) :: uu!for the energies
   TYPE(times_freqs), INTENT(in) :: tf !for grid on imaginary time
   CHARACTER(LEN=256), INTENT(in) ::  prefix!to designate the PW files

   INTEGER :: l_blk, nbegin, nend
   INTEGER :: it,iv,ic, iw, jw
   TYPE(upper_states) :: us
   TYPE(polaw) :: pp
   TYPE(cprim_prod) :: cpp
   LOGICAL :: ok_read
   REAL(kind=DP), ALLOCATABLE :: exp_table(:)
   REAL(kind=DP), ALLOCATABLE :: cpmat_tmp(:,:)

   write(stdout,*) 'Routine create_polarization_upper'


!read-in upper states
   call initialize_memory(us)
   call read_data_pw_upper_states(us,prefix)

   allocate(exp_table(us%nums_reduced))





!loop on time

    l_blk= (tf%n+1)/nproc
    if(l_blk*nproc < (tf%n+1)) l_blk = l_blk+1
    nbegin=mpime*l_blk
    nend=nbegin+l_blk-1
    
    do it=nbegin,nend


       if(it<= tf%n) then
!read polarization
          call initialize_polaw(pp)
          call read_polaw(it,pp,.false.,.false.)
          do iv=1,uu%nums_occ(1)
 !loop on v
             write(stdout,*) 'STATE', iv
             FLUSH(stdout)

!set table of exponentials

             do ic=1,us%nums_reduced
                exp_table(ic)=exp((uu%ene(iv,1)-us%ene(ic))*tf%times(it))
             enddo

             call initialize_memory(cpp)

   !!read in Svci 
             cpp%cprim=iv
             call read_data_pw_cprim_prod(cpp, prefix, .true., ok_read, .true.,.true.)
!if required allocate polaw
                 
!!!!the following for using blas routine
              allocate(cpmat_tmp(us%nums_reduced,cpp%numpw))
              call mytranspose(cpp%cpmat, cpp%numpw, cpmat_tmp, us%nums_reduced, cpp%numpw, us%nums_reduced)
              do iw=1,cpp%numpw
                 cpmat_tmp(:,iw)=cpmat_tmp(:,iw)*exp_table(:)
              enddo

!!!
    !calculate term
!              do ic=1,cpp%nums_cond
!                 do jw=1,cpp%numpw
!                    do iw=1,cpp%numpw
!                       pp%pw(iw,jw)=pp%pw(iw,jw)+cpp%cpmat(iw,ic)*cpp%cpmat(jw,ic)*exp_table(ic)
!                    enddo
!                 enddo
!              enddo

!factor 2 for spin multiplicity
              call dgemm('N','N',cpp%numpw,cpp%numpw,us%nums_reduced,2.d0,cpp%cpmat,cpp%numpw,&
                   &cpmat_tmp,us%nums_reduced,1.d0,pp%pw,cpp%numpw)

              deallocate(cpmat_tmp)
             call free_memory(cpp)
          enddo
 !write polaw on file
!

          call write_polaw(pp,.false.)
          call free_memory_polaw(pp)
          else
!just parallelel reading of file
             do iv=1,uu%nums_occ(1)
                call initialize_memory(cpp)
                cpp%cprim=iv
                call read_data_pw_cprim_prod(cpp, prefix, .true., ok_read, .true.,.true.)
                call free_memory(cpp)
             enddo
          endif
       enddo

       deallocate(exp_table)
       call free_memory(us)
   return

 END SUBROUTINE create_polarization_upper




  SUBROUTINE calculate_w_g_l(vp,pp,ww,xc_together,l_head_epsilon,head,inv_epsi, &
                    &l_wing_epsilon, awing, awing_c, l_verbose)
!this subroutine calculates W=(1+vp)^{-1}v
!this is meaningful only on frequency domain
!lapack routines are used
!it use exteded basis for G=0 term
!version for lanczos chain scheme

    USE io_global,            ONLY : stdout
    USE basic_structures,     ONLY : v_pot, head_epsilon

    implicit none
    TYPE(v_pot)  :: vp!coulomb potential
    TYPE(polaw)  :: pp!polarization on imaginary frequency, destroyed on exit
    TYPE(polaw)  :: ww!dressed interaction to be calculated
    LOGICAL :: xc_together!if true the entire W is taken, otherwise W-v
    LOGICAL :: l_head_epsilon!if true add to the symmetrized form  of the dielectric matrix
                             !the head terms
    REAL(kind=DP) :: head(3)!term (G=0,G=0) of the symmetric dielectric matrix
      REAL(kind=DP), INTENT(out) :: inv_epsi!head of the inverse dielectric matrix    
    LOGICAL, INTENT(in) :: l_wing_epsilon!if true calculate the wings of the symmetrized dielectric matrix
    REAL(kind=DP), DIMENSION(pp%numpw,3) :: awing!the terms A_ij wing_j
    REAL(kind=DP), DIMENSION(pp%numpw,3) :: awing_c!the terms A_ij wing_j
    LOGICAL, INTENT(in) :: l_verbose

    INTEGER iw,jw,kw,ipol
    REAL(kind=DP), ALLOCATABLE, DIMENSION(:,:) :: dtmp!temporary array
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
    INTEGER :: info
    REAL(kind=DP),ALLOCATABLE, DIMENSION(:) :: work
    INTEGER :: lwork
    REAL(kind=DP) sca
    REAL(kind=DP) :: workd(1)

    REAL(kind=DP) alpha
    REAL(kind=DP) head_v
    REAL(kind=DP), ALLOCATABLE :: pw_save(:,:)

!deallocate if the case
    call free_memory_polaw(ww)


!check and set
   if(pp%ontime) then
      write(stdout,*) 'Routine calculate_w: frequencies required'
      stop
   endif
   if(pp%numpw /= vp%numpw) then
      write(stdout,*) 'Routine calculate_w: basis set does not correspond',pp%numpw,vp%numpw
      stop
   endif

   ww%ontime=.false.
   ww%time=pp%time
   ww%label=pp%label
   ww%numpw=pp%numpw
   allocate(ww%pw(ww%numpw,ww%numpw))

   allocate(dtmp(ww%numpw,ww%numpw))
   allocate(ipiv(ww%numpw))

   ww%pw(:,:)=0.d0

   allocate(pw_save(pp%numpw,pp%numpw))
   do ipol=1,3
      if(l_verbose) write(stdout,*) 'MAX P:', maxval(pp%pw(:,:)), 'MIN P:', minval(pp%pw(:,:))
      FLUSH(stdout)
      if(ipol==1) then
         pw_save(:,:)=pp%pw(:,:)
      else
         pp%pw(:,:)=pw_save(:,:)
      endif
      dtmp(:,:)=0.d0

      head_v=vp%vmat(vp%numpw,vp%numpw)

!DEBUG
      if(l_verbose) write(stdout,*) 'IRR POLA HEAD', pp%pw(pp%numpw,pp%numpw)
      if(l_verbose) write(stdout,*) 'IRR POLA FIRST', pp%pw(1,1), l_wing_epsilon

!      vp%vmat(vp%numpw,1:vp%numpw)=0.d0 !ATTENZIONE DEBUG
!      vp%vmat(1:vp%numpw,vp%numpw)=0.d0
      
      pp%pw(pp%numpw,1:pp%numpw)=0.d0
      pp%pw(1:pp%numpw,pp%numpw)=0.d0

      if(l_wing_epsilon) then!ATTENZIONE
!  0',i   ',0'                                                                   
         do iw=1,ww%numpw-1
            pp%pw(ww%numpw,iw)=awing(iw,ipol)
            pp%pw(iw,ww%numpw)=pp%pw(ww%numpw,iw)
         enddo

      endif


!symmetric case calculates -v^1/2 P v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,-1.d0*dble(pp%factor),&
           &vp%vmat,ww%numpw,pp%pw,ww%numpw,0.d0,dtmp,ww%numpw)
      pp%pw(1:ww%numpw,1:ww%numpw)=dtmp(1:ww%numpw,1:ww%numpw)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & pp%pw,ww%numpw,vp%vmat,ww%numpw,0.d0,dtmp,ww%numpw)
   




!if required add the head
      if(l_head_epsilon) then


!term O',0'  0',i  i,0'
         dtmp(ww%numpw,ww%numpw)= head(ipol)
         if(l_verbose) write(stdout,*) 'APPLYING HEAD', head!DEBUG
      endif

      do iw=1,ww%numpw
         dtmp(iw,iw)=dtmp(iw,iw)+1.d0
      enddo


!inverse zmat

      if(l_verbose) write(stdout,*) 'MAX D:', maxval(dtmp(:,:)), 'MIN D', minval(dtmp(:,:))

      if(l_verbose) write(stdout,*) 'Before inversion'
      FLUSH(stdout)
      call dgetrf(ww%numpw,ww%numpw,dtmp,ww%numpw,ipiv,info)
      if(info /= 0) then
         write(stdout,*) 'Routine calculate_w: problem with dgetrf :', info
         stop
      endif
      if(l_verbose) write(stdout,*) 'Before inversion2'
      FLUSH(stdout)

      call dgetri(ww%numpw,dtmp,ww%numpw,ipiv,workd,-1,info)
      if(l_verbose) write(stdout,*) 'Dimension', workd,ww%numpw,info!ATTENZIONE
      FLUSH(stdout)
      allocate(work(int(workd(1))))
      call dgetri(ww%numpw,dtmp,ww%numpw,ipiv,work,int(workd(1)),info)

      if(l_verbose) write(stdout,*) 'Out of dgetri'
      FLUSH(stdout)

      if(l_verbose)  write(stdout,*) 'MAX D1:', maxval(dtmp(:,:)), 'MIN D1:', minval(dtmp(:,:))

      if(info /= 0) then
         write(stdout,*) 'Routine calculate_w: problem with zgetri :', info
         stop
      endif


      if(.not.xc_together) then
         do iw=1,ww%numpw
            dtmp(iw,iw)=dtmp(iw,iw)-1.d0
         enddo
      endif


      inv_epsi=0.d0

!term i,j
      !term 0',0'
      inv_epsi=dtmp(ww%numpw,ww%numpw)
      
      write(stdout,*) 'INV EPSI G=0,G=0', inv_epsi, head_v
      FLUSH(stdout)
      
      vp%vmat(vp%numpw,vp%numpw)=head_v
      dtmp(ww%numpw,1:ww%numpw-1)=0.d0
      dtmp(1:ww%numpw-1,ww%numpw)=0.d0
    
!DEBUG   
!   dtmp(:,:)=0.d0
!   do iw=1,ww%numpw
!      dtmp(iw,iw)=inv_epsi
!   enddo
!   dtmp(ww%numpw,ww%numpw)=0.d0

!DEBUG
      if(l_verbose) write(stdout,*) 'MAX D2:', maxval(dtmp(:,:)), 'MIN D2:', minval(dtmp(:,:))
      FLUSH(stdout)

!calculates v^1/2 (e-1-1)v^1/2
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & vp%vmat,ww%numpw,dtmp,ww%numpw,0.d0,pp%pw,ww%numpw)
      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
           & pp%pw,ww%numpw,vp%vmat,ww%numpw,1.d0,ww%pw,ww%numpw)
  
      if(l_verbose) write(stdout,*) 'MAX W:', maxval(ww%pw(:,:)), 'MIN W:', minval(ww%pw(:,:))
      FLUSH(stdout)
   

      deallocate(work)
   enddo
   deallocate(pw_save)
   ww%pw(:,:)=ww%pw(:,:)/3.d0
!   ww%pw(1:ww%numpw-1,1:ww%numpw-1)=0.d0
!   ww%pw(:,ww%numpw)=0.d0
   ww%factor=(1.d0,0.d0)

   if(l_verbose) write(stdout,*) 'MAX:', maxval(ww%pw(:,:)), 'MIN:', minval(ww%pw(:,:))
   FLUSH(stdout)
   

   deallocate(dtmp,ipiv)
   return
 END SUBROUTINE calculate_w_g_l





  END MODULE polarization

