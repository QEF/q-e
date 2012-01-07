!P.Umari
!Program GWW

  MODULE green_function
!this module descibes the green function in imaginary time/frequency
!and contains subroutine to read/write from disk and to create

    USE kinds, ONLY : DP

      TYPE green
!this structure describe a generic green function
!usually in the space of wanniers
        INTEGER :: label!label to read/write to disk
        LOGICAL :: ontime!if .true. is on imaginary time, otherwise frequency
        REAL(kind=DP) :: time!imaginary time or frequency
        INTEGER :: nums!number of states
        LOGICAL :: zero_time_neg!if .true. the green function at t=0 is calculated as a negative time one
        COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: gf!green function
        LOGICAL :: l_part!if true the matrix is written as a real matrix times a sign
        REAL(kind=DP), DIMENSION(:,:), POINTER :: gf_p!green function
        COMPLEX(kind=DP) :: factor !complex factor for gf_p
     END TYPE green
  CONTAINS

    SUBROUTINE initialize_green(gr)
      implicit none
      TYPE(green) gr
       nullify(gr%gf)
       nullify(gr%gf_p)
      return
    END SUBROUTINE



    SUBROUTINE free_memory_green(gr)
!this subroutine deallocates the green descriptor
      implicit none
      TYPE(green) gr
      if(associated(gr%gf)) deallocate(gr%gf)
      nullify(gr%gf)
      if(associated(gr%gf_p)) deallocate(gr%gf_p)
      nullify(gr%gf_p)
      return
    END SUBROUTINE

    SUBROUTINE create_green(gr,wu,time,debug,zero_time_neg,l_hf_energies,ene_hf)
!this subroutine creates a green function on imagynary time
!on the basis of wanniers:
!the KS energies are fixed so that the fermi level is at 0
!  G_{i,j}=i*\sum_v U^{+}_{v,i}*U_{j,v}*exp(e_v*t)  t>=0
!         =-i*\sum_c U^{+}_{c,i}*U_{j,c}*exp(e_c*t) t<0
!if required uses HF energies


      USE kinds, ONLY : DP
      USE io_global, ONLY : stdout
      USE basic_structures, ONLY : wannier_u

      implicit none

      TYPE(green) :: gr!the green function on output
      TYPE(wannier_u) :: wu!data on U and e_i
      REAL(kind=DP) :: time!imaginary time
      LOGICAL :: debug!if true print debug information on stdout
      LOGICAL :: zero_time_neg!if true and time==0, the negative form is forced
      LOGICAL, INTENT(in) :: l_hf_energies!if true uses HF energies
      REAL(kind=DP), INTENT(in) :: ene_hf(:)


      INTEGER iw,jw,kw
      REAL(kind=DP) :: offset

!calculates energy offset

      if(.not.l_hf_energies) then
         if(wu%nums > wu%nums_occ) then
            offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
         else
            offset=-wu%ene(wu%nums_occ)
         endif
      else
         if(wu%nums > wu%nums_occ) then
            offset=-(ene_hf(wu%nums_occ+1)+ene_hf(wu%nums_occ))/2.d0
         else
            offset=-ene_hf(wu%nums_occ)
         endif
      endif
!sets data and allocate
!      call free_memory_green(gr)
      gr%nums=wu%nums
      allocate(gr%gf(gr%nums,gr%nums))
      gr%gf(:,:)=(0.d0,0.d0)
      gr%ontime=.TRUE.
      gr%time=time
      gr%zero_time_neg = zero_time_neg

      if((time < 0.d0).or. (time==0.d0 .and. zero_time_neg )) then
! only conduction states
        do iw=1,gr%nums
          do jw=iw,gr%nums
             do  kw=wu%nums_occ+1,wu%nums
                if(.not.l_hf_energies) then
                   gr%gf(iw,jw)=gr%gf(iw,jw)+wu%umat(kw,iw)*conjg(wu%umat(kw,jw))* &
                        & exp((wu%ene(kw)+offset)*time)
                else
                   gr%gf(iw,jw)=gr%gf(iw,jw)+wu%umat(kw,iw)*conjg(wu%umat(kw,jw))* &
                        & exp((ene_hf(kw)+offset)*time)
                endif
                 if(debug) then
                   write(stdout,*) 'Create green:' ,time,iw,jw,wu%ene(kw),wu%nums_occ+1,wu%umat(jw,kw)
                 endif
             enddo
             gr%gf(jw,iw)=conjg(gr%gf(iw,jw))
             gr%gf(iw,jw)=(0.d0,-1.d0)*gr%gf(iw,jw)
             if(iw /= jw) gr%gf(jw,iw)=(0.d0,-1.d0)*gr%gf(jw,iw)
          enddo
        enddo
      else
! only valence states
        do iw=1,gr%nums
          do jw=iw,gr%nums
             do  kw=1,wu%nums_occ
                if(.not. l_hf_energies) then
                   gr%gf(iw,jw)=gr%gf(iw,jw)+wu%umat(kw,iw)*conjg(wu%umat(kw,jw))* &
                        & exp((wu%ene(kw)+offset)*time)
                else
                   gr%gf(iw,jw)=gr%gf(iw,jw)+wu%umat(kw,iw)*conjg(wu%umat(kw,jw))* &
                        & exp((ene_hf(kw)+offset)*time)
                endif
                 if(debug) then
                   write(stdout,*) 'Create green:' ,time,iw,jw,wu%ene(kw),wu%umat(kw,iw),wu%umat(kw,jw)
                 endif
             enddo
             gr%gf(jw,iw)=conjg(gr%gf(iw,jw))
             gr%gf(iw,jw)=(0.d0,1.d0)*gr%gf(iw,jw)
             if(iw /= jw) gr%gf(jw,iw)=(0.d0,1.d0)*gr%gf(jw,iw)
             if(debug) write(stdout,*) 'Create green2:', iw,jw, gr%gf(iw,jw), offset
          enddo
        enddo

      endif


      return

   END SUBROUTINE

    SUBROUTINE create_green_part(gr,wu,time,debug,zero_time_neg,l_hf_energies,ene_hf)
!this subroutine creates a green function on imagynary time
!on the basis of wanniers:
!the KS energies are fixed so that the fermi level is at 0
!  G_{i,j}=i*\sum_v U^{+}_{v,i}*U_{j,v}*exp(e_v*t)  t>=0
!         =-i*\sum_c U^{+}_{c,i}*U_{j,c}*exp(e_c*t) t<0
!if required uses HF energies
!it uses consider a real part plus a factor


      USE kinds, ONLY : DP
      USE io_global, ONLY : stdout
      USE basic_structures, ONLY : wannier_u

      implicit none

      TYPE(green) :: gr!the green function on output
      TYPE(wannier_u) :: wu!data on U and e_i
      REAL(kind=DP) :: time!imaginary time
      LOGICAL :: debug!if true print debug information on stdout
      LOGICAL :: zero_time_neg!if true and time==0, the negative form is forced
      LOGICAL, INTENT(in) :: l_hf_energies!if true uses HF energies
      REAL(kind=DP), INTENT(in) :: ene_hf(:)


      INTEGER iw,jw,kw
      REAL(kind=DP) :: offset


      call free_memory_green(gr)

      gr%l_part=.true.

!calculates energy offset

      if(.not.l_hf_energies) then
         if(wu%nums > wu%nums_occ) then
            offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
         else
            offset=-wu%ene(wu%nums_occ)
         endif
      else
         if(wu%nums > wu%nums_occ) then
            offset=-(ene_hf(wu%nums_occ+1)+ene_hf(wu%nums_occ))/2.d0
         else
            offset=-ene_hf(wu%nums_occ)
         endif
      endif
!sets data and allocate
!      call free_memory_green(gr)
      gr%nums=wu%nums
      allocate(gr%gf_p(gr%nums,gr%nums))
      gr%gf_p(:,:)=0.d0
      gr%ontime=.TRUE.
      gr%time=time
      gr%zero_time_neg = zero_time_neg

      if((time < 0.d0).or. (time==0.d0 .and. zero_time_neg )) then
! only conduction states
        do iw=1,gr%nums
          do jw=iw,gr%nums
             do  kw=wu%nums_occ+1,wu%nums
                if(.not.l_hf_energies) then
                   gr%gf_p(iw,jw)=gr%gf_p(iw,jw)+dble(wu%umat(kw,iw))*dble(wu%umat(kw,jw))* &
                        & exp((wu%ene(kw)+offset)*time)
                  ! if(abs(aimag(wu%umat(kw,iw)) >=1.d-6)) write(stdout,*) 'PROBLEMA'
                else
                   gr%gf_p(iw,jw)=gr%gf_p(iw,jw)+dble(wu%umat(kw,iw))*dble(wu%umat(kw,jw))* &
                        & exp((ene_hf(kw)+offset)*time)
                endif
                 if(debug) then
                   write(stdout,*) 'Create green:' ,time,iw,jw,wu%ene(kw),wu%nums_occ+1,wu%umat(jw,kw)
                 endif
             enddo
             gr%gf_p(jw,iw)=gr%gf_p(iw,jw)
             gr%factor=(0.d0,-1.d0)
          enddo
        enddo
      else
! only valence states
        do iw=1,gr%nums
          do jw=iw,gr%nums
             do  kw=1,wu%nums_occ
                if(.not. l_hf_energies) then
                   gr%gf_p(iw,jw)=gr%gf_p(iw,jw)+dble(wu%umat(kw,iw))*dble(wu%umat(kw,jw))* &
                        & exp((wu%ene(kw)+offset)*time)
                else
                   gr%gf_p(iw,jw)=gr%gf_p(iw,jw)+dble(wu%umat(kw,iw))*dble(wu%umat(kw,jw))* &
                        & exp((ene_hf(kw)+offset)*time)
                endif
                 if(debug) then
                   write(stdout,*) 'Create green:' ,time,iw,jw,wu%ene(kw),wu%umat(kw,iw),wu%umat(kw,jw)
                 endif
             enddo
             gr%gf_p(jw,iw)=gr%gf_p(iw,jw)
             gr%factor=(0.d0,1.d0)
             if(debug) write(stdout,*) 'Create green2:', iw,jw, gr%gf_p(iw,jw), offset
          enddo
        enddo

      endif


      return

    END SUBROUTINE create_green_part






   SUBROUTINE write_green(gr, debug)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_files,             ONLY : find_free_unit, prefix
    implicit none
    TYPE(green) :: gr!the green function to be written
    LOGICAL :: debug!if true print formatted file

    INTEGER :: iw, jw, iung
    CHARACTER(5) :: nfile
    if(gr%label > 0 .or. (gr%label == 0 .and. .not.gr%zero_time_neg)) then
      write(nfile,'(5i1)') &
         & gr%label/10000,mod(gr%label,10000)/1000,mod(gr%label,1000)/100,mod(gr%label,100)/10,mod(gr%label,10)
      iung = find_free_unit()
      if(.not. debug) then
        open( unit=iung, file='green.'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iung, file='green.'// nfile, status='unknown',form='formatted')
      endif
    else
      write(nfile,'(5i1)') &
      & -gr%label/10000,mod(-gr%label,10000)/1000,mod(-gr%label,1000)/100,mod(-gr%label,100)/10,mod(-gr%label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file='green.-'// nfile, status='unknown',form='unformatted')
      else
        open( unit=iung, file='green.-'// nfile, status='unknown',form='formatted')
      endif
    endif
    if(.not.debug) then
      write(iung) gr%label
      write(iung) gr%ontime
      write(iung) gr%time
      write(iung) gr%nums
      write(iung) gr%zero_time_neg
      write(iung) gr%l_part
      write(iung) gr%factor
      if(.not.gr%l_part) then
         do iw=1,gr%nums
            write(iung)  gr%gf(1:gr%nums,iw)
         enddo
      else
         do iw=1,gr%nums
            write(iung)  gr%gf_p(1:gr%nums,iw)
         enddo
      endif
   else
      write(iung,*) gr%label
      write(iung,*) gr%ontime
      write(iung,*) gr%time
      write(iung,*) gr%nums
      write(iung,*) gr%zero_time_neg
      write(iung,*) gr%l_part
      write(iung,*) gr%factor
      if(.not.gr%l_part) then
         do iw=1,gr%nums
            do jw=1,gr%nums
               write(iung,*)  gr%gf(jw,iw)
            enddo
         enddo
      else
         do iw=1,gr%nums
            do jw=1,gr%nums
               write(iung,*)  gr%gf_p(jw,iw)
            enddo
         enddo
      endif
   endif
   close(iung)
   return
 END SUBROUTINE write_green

   SUBROUTINE  read_green(label, gr, debug,zero_time_neg)
!this subroutine reads the green function from disk
!the file name is taken from the label

    USE io_files,             ONLY : find_free_unit, prefix
    implicit none
    TYPE(green) :: gr!the green function to be read
    INTEGER :: label! the label identifing the required green function
    LOGICAL :: debug!if true print formatted file
    LOGICAL :: zero_time_neg !if true and time == 0, a negative kind of green function is considered

    INTEGER :: iw, jw, iung
    CHARACTER(5) :: nfile

!first deallocate
    call free_memory_green(gr)
    if(label > 0 .or. (label == 0 .and. .not.zero_time_neg))  then
      write(nfile,'(5i1)') label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file='green.'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file='green.'// nfile, status='old',form='formatted')
      endif
    else
      write(nfile,'(5i1)') -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
      iung = find_free_unit()
      if(.not.debug) then
        open( unit=iung, file='green.-'// nfile, status='old',form='unformatted')
      else
        open( unit=iung, file='green.-'// nfile, status='old',form='formatted')
      endif
    endif
    if(.not.debug) then
      read(iung) gr%label
      read(iung) gr%ontime
      read(iung) gr%time
      read(iung) gr%nums
      read(iung) gr%zero_time_neg
      read(iung) gr%l_part
      read(iung) gr%factor

!now allocate
      if(.not. gr%l_part) then
         allocate(gr%gf(gr%nums,gr%nums))
         nullify(gr%gf_p)
      else
         allocate(gr%gf_p(gr%nums,gr%nums))
         nullify(gr%gf)
      endif
      if(.not. gr%l_part) then
         do iw=1,gr%nums
            read(iung)  gr%gf(1:gr%nums,iw)
         enddo
      else
         do iw=1,gr%nums
            read(iung)  gr%gf_p(1:gr%nums,iw)
         enddo
      endif
   else
      read(iung,*) gr%label
      read(iung,*) gr%ontime
      read(iung,*) gr%time
      read(iung,*) gr%nums
      read(iung,*) gr%zero_time_neg
      read(iung,*) gr%l_part
      read(iung,*) gr%factor

!now allocate

      if(.not. gr%l_part) then
         allocate(gr%gf(gr%nums,gr%nums))
         nullify(gr%gf_p)
      else
         allocate(gr%gf_p(gr%nums,gr%nums))
         nullify(gr%gf)
      endif

      if(.not. gr%l_part) then
         do iw=1,gr%nums
            do jw=1,gr%nums
               read(iung,*)  gr%gf(jw,iw)
            enddo
         enddo
      else
         do iw=1,gr%nums
            do jw=1,gr%nums
               read(iung,*)  gr%gf_p(jw,iw)
            enddo
         enddo
      endif
   endif
   close(iung)
   return
 END SUBROUTINE read_green


  END MODULE green_function


