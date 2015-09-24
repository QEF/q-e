!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


  MODULE compact_product
!this module describes the contracted products O^{P}_{n,kl}U_{ki}=Q^{P}_{n,l,i}
!and of the contracted products  O^{P}_{i,km}O^{P}_{j,ln}U^{+}_{vk}U_{lv}U^{+}_{cm}U_{nc}

    USE kinds, ONLY : DP

     TYPE contraction
!this structure described the localized and normalized products of wanniers  w_P
      INTEGER :: numpw!number of wannier-products
      INTEGER :: nums!number of KS or wannier states
      INTEGER :: max_i!maximum number of KS states to be addresses
      INTEGER, DIMENSION(:), POINTER :: numl!array for number of functionsl (numl,:)
      INTEGER, DIMENSION(:,:), POINTER :: l !array for l indices (index,:,:)
      COMPLEX(kind=DP),DIMENSION(:,:,:), POINTER :: q!contraction terms
    END TYPE contraction

!for treating large systems contraction has been split in two parts

    TYPE contraction_index
!this structure described the localized and normalized products of wanniers  w_P
!index part
      INTEGER :: numpw!number of wannier-products
      INTEGER :: nums!number of KS or wannier states
      INTEGER :: max_i!maximum number of KS states to be addresses
      INTEGER, DIMENSION(:), POINTER :: numl!array for number of functionsl (numl,:)
      INTEGER, DIMENSION(:,:), POINTER :: l !array for l indices (index,:,:)
   END TYPE contraction_index

    TYPE contraction_state
!this structure described the localized and normalized products of wanniers  w_P
      INTEGER :: numpw!number of wannier-products
      INTEGER :: nums!number of KS or wannier states
      INTEGER :: max_i!maximum number of KS states to be addresses
      INTEGER :: state!state for which the contraction corresponds
      REAL(kind=DP),DIMENSION(:,:), POINTER :: q!contraction terms
   END TYPE contraction_state




      
     TYPE contraction_pola
!this structure described the localized and normalized products of wanniers with U matrices
      INTEGER :: numpw!number of wannier-products
      INTEGER :: nums!number of KS or wannier states
      INTEGER :: nums_occ!number of occupied states
      COMPLEX(kind=DP),DIMENSION(:,:,:), POINTER :: ou!contraction terms
    END TYPE contraction_pola     

    TYPE contraction_pola_state
!this structure described the localized and normalized products of wanniers with U matrices
!just for one occupied state
      INTEGER :: numpw!number of wannier-products
      INTEGER :: nums!number of KS or wannier states
      INTEGER :: nums_occ!number of occupied states
      INTEGER :: state!occupied state relative to this data
      REAL(kind=DP),DIMENSION(:,:), POINTER :: ou!contraction terms
   END TYPE contraction_pola_state



  CONTAINS

    SUBROUTINE  free_memory_contraction_pola(cp)
      implicit none
      
      TYPE(contraction_pola) :: cp

      if(associated(cp%ou)) then
        deallocate(cp%ou)
        nullify(cp%ou)
      endif
      
      return
    END SUBROUTINE free_memory_contraction_pola



   SUBROUTINE  free_memory_contraction_pola_state(cp)
      implicit none

      TYPE(contraction_pola_state) :: cp

      if(associated(cp%ou)) then
        deallocate(cp%ou)
        nullify(cp%ou)
      endif

      return
    END SUBROUTINE free_memory_contraction_pola_state



    SUBROUTINE free_memory_contraction(cr)
      implicit none

      TYPE(contraction) :: cr

      if(associated(cr%numl)) then
        deallocate(cr%numl)
        nullify(cr%numl)
      endif
      if(associated(cr%l)) then
        deallocate(cr%l)
        nullify(cr%l)
      endif
      if(associated(cr%q)) then
        deallocate(cr%q)
        nullify(cr%q)
      endif

      return
    END SUBROUTINE

    SUBROUTINE free_memory_contraction_index(cr)
      implicit none

      TYPE(contraction_index) :: cr

      if(associated(cr%numl)) then
        deallocate(cr%numl)
        nullify(cr%numl)
      endif
      if(associated(cr%l)) then
        deallocate(cr%l)
        nullify(cr%l)
      endif

      return
    END SUBROUTINE

    SUBROUTINE free_memory_contraction_state(cr)
      implicit none

      TYPE(contraction_state) :: cr

      if(associated(cr%q)) then
        deallocate(cr%q)
        nullify(cr%q)
      endif

      return
    END SUBROUTINE


   SUBROUTINE write_contraction(cr, options)
!this subroutine writes the contracted products on disk
!in parallel case only ionode writes

   
    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction) :: cr!the contraction descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun

    if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction', status='unknown',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction', status='unknown',form='formatted')
       endif

       if(.not.options%debug) then
          write(iun) cr%numpw
          write(iun) cr%nums
          write(iun) cr%max_i
          write(iun) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             write(iun) cr%l(1:cr%numl(iw),iw)
          enddo
          do iw=1,cr%numpw
             write(iun) cr%q(iw,1:cr%numl(iw),1:cr%max_i)
          enddo
       else
          write(iun,*) cr%nums
          write(iun,*) cr%max_i
          write(iun,*) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                write(iun,*) cr%l(jw,iw)
             enddo
          enddo
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                do kw=1,cr%max_i
                   write(iun,*) cr%q(iw,jw,kw)
                enddo
             enddo
          enddo
       endif
       close(iun)
    endif
    return         

  END SUBROUTINE  write_contraction

  SUBROUTINE read_contraction(cr, options)
!this subroutine reads the contracted products from disk
!in parallel case only ionode reads


    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : stdout, ionode, ionode_id
    USE mp,                 ONLY : mp_bcast
    USE mp_world,           ONLY : world_comm
    USE io_files,             ONLY : prefix,tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction) :: cr!the contraction descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun, ii
    INTEGER maxl

    if(ionode) then 
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction', status='old',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction', status='old',form='formatted')
       endif
    endif

    !call free_memory_contraction(cr)



    if(ionode) then
       if(.not.options%debug) then
          read(iun) cr%numpw
          read(iun) cr%nums
          read(iun) cr%max_i
       else
          read(iun,*) cr%numpw
          read(iun,*) cr%nums
          read(iun,*) cr%max_i
       endif
    endif

    call mp_bcast(cr%numpw, ionode_id,world_comm)
    call mp_bcast(cr%nums, ionode_id,world_comm)
    call mp_bcast(cr%max_i, ionode_id,world_comm)

    maxl=cr%numpw!TEMPORARY SOLUTION ATTENZIONE
    maxl=cr%nums
    allocate(cr%numl(cr%numpw))
    allocate(cr%l(maxl,cr%numpw))
    allocate(cr%q(cr%numpw,maxl,cr%max_i))

    if(ionode) then
       write(stdout,*) 'CR-READ',cr%numpw,maxl,cr%max_i

       if(.not.options%debug) then    
          read(iun) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             read(iun) cr%l(1:cr%numl(iw),iw)
          enddo
           write(stdout,*) 'CR-READ L'
          do iw=1,cr%numpw
             read(iun) cr%q(iw,1:cr%numl(iw),1:cr%max_i)
             write(stdout,*) 'CR-READ Q', iw
          enddo
       else
          read(iun,*) cr%nums
          read(iun,*) cr%max_i
          read(iun,*) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                read(iun,*) cr%l(jw,iw)
             enddo
          enddo
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                do kw=1,cr%max_i
                   read(iun,*) cr%q(iw,jw,kw)
                enddo
             enddo
          enddo
       endif
    endif
    call mp_bcast(cr%numl(:), ionode_id,world_comm)
    call mp_bcast(cr%l(:,:), ionode_id,world_comm)
    write(stdout,*) 'CR-SEND L'
    do ii=1,cr%max_i
       call mp_bcast(cr%q(:,:,ii), ionode_id,world_comm)
       write(stdout,*) 'CR-SEND Q',ii
    enddo
    if(ionode) close(iun)
  END SUBROUTINE read_contraction




  SUBROUTINE write_contraction_index(cr, options)
!this subroutine writes the contracted products on disk
!in parallel case only ionode writes

    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_index) :: cr!the contraction index descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun

    if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_index', status='unknown',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_index', status='unknown',form='formatted')
       endif

       if(.not.options%debug) then
          write(iun) cr%numpw
          write(iun) cr%nums
          write(iun) cr%max_i
          write(iun) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             write(iun) cr%l(1:cr%numl(iw),iw)
          enddo
       else
          write(iun,*) cr%nums
          write(iun,*) cr%max_i
          write(iun,*) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                write(iun,*) cr%l(jw,iw)
             enddo
          enddo
       endif
       close(iun)
    endif
    return

  END SUBROUTINE  write_contraction_index

  SUBROUTINE read_contraction_index(cr, options)
!this subroutine reads the contracted products from disk
!in parallel case only ionode reads


    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : stdout, ionode, ionode_id
    USE mp,                 ONLY : mp_bcast
    USE mp_world,           ONLY : world_comm
    USE io_files,             ONLY : prefix,tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_index) :: cr!the contraction descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun, ii
    INTEGER maxl

    if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_index', status='old',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_index', status='old',form='formatted')
       endif
    endif

    !call free_memory_contraction(cr)



    if(ionode) then
       if(.not.options%debug) then
          read(iun) cr%numpw
          read(iun) cr%nums
          read(iun) cr%max_i
       else
          read(iun,*) cr%numpw
          read(iun,*) cr%nums
          read(iun,*) cr%max_i
       endif
    endif

    call mp_bcast(cr%numpw, ionode_id,world_comm)
    call mp_bcast(cr%nums, ionode_id,world_comm)
    call mp_bcast(cr%max_i, ionode_id,world_comm)

    maxl=cr%nums
    allocate(cr%numl(cr%numpw))
    allocate(cr%l(maxl,cr%numpw))
    
    if(ionode) then
       write(stdout,*) 'CR-READ',cr%numpw,maxl,cr%max_i

       if(.not.options%debug) then
          read(iun) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             read(iun) cr%l(1:cr%numl(iw),iw)
          enddo
           write(stdout,*) 'CR-READ L'
        else
          read(iun,*) cr%nums
          read(iun,*) cr%max_i
          read(iun,*) cr%numl(1:cr%numpw)
          do iw=1,cr%numpw
             do jw=1,cr%numl(iw)
                read(iun,*) cr%l(jw,iw)
             enddo
          enddo
       endif
    endif
    call mp_bcast(cr%numl(:), ionode_id,world_comm)
    call mp_bcast(cr%l(:,:), ionode_id,world_comm)
    write(stdout,*) 'CR-SEND L'
    if(ionode) close(iun)
  END SUBROUTINE read_contraction_index

   SUBROUTINE write_contraction_state(cri,crs, options)
!this subroutine writes the contracted products on disk
!in parallel case only ionode writes

     USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_index), INTENT(in) :: cri!the contraction index descriptor
    TYPE(contraction_state)             :: crs!the contraction state to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun
    CHARACTER(5) :: nfile

    write(nfile,'(5i1)') &
         & crs%state/10000,mod(crs%state,10000)/1000,mod(crs%state,1000)/100,mod(crs%state,100)/10,mod(crs%state,10)
    iun = find_free_unit()
    if(.not. options%debug) then
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='unknown',form='unformatted')
    else
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='unknown',form='formatted')
    endif

    if(.not.options%debug) then
       write(iun) crs%numpw
       write(iun) crs%nums
       write(iun) crs%max_i
       write(iun) crs%state
       do iw=1,crs%nums
          write(iun) crs%q(1:cri%numpw,iw)
       enddo
    else
       write(iun,*) crs%numpw
       write(iun,*) crs%nums
       write(iun,*) crs%max_i
       write(iun,*) crs%state

       do iw=1,crs%numpw
             do jw=1,cri%nums
                write(iun,*) crs%q(iw,jw)
             enddo
          enddo
       endif
       close(iun)
       return

     END SUBROUTINE  write_contraction_state


   SUBROUTINE read_contraction_state(cri,crs, options)
!this subroutine writes the contracted products on disk
!in parallel case only ionode writes

     USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_index), INTENT(in) :: cri!the contraction index descriptor
    TYPE(contraction_state)             :: crs!the contraction state to be read from file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun
    CHARACTER(5) :: nfile
    INTEGER :: maxl


    write(nfile,'(5i1)') &
         & crs%state/10000,mod(crs%state,10000)/1000,mod(crs%state,1000)/100,mod(crs%state,100)/10,mod(crs%state,10)
    iun = find_free_unit()
    if(.not. options%debug) then
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='old',form='unformatted')
    else
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='old',form='formatted')
    endif

    if(.not.options%debug) then
       read(iun) crs%numpw
       read(iun) crs%nums
       read(iun) crs%max_i
       read(iun) crs%state
    else
       read(iun,*) crs%numpw
       read(iun,*) crs%nums
       read(iun,*) crs%max_i
       read(iun,*) crs%state
    endif

    maxl=crs%nums
    allocate(crs%q(crs%numpw,maxl))

    if(.not.options%debug) then
       do iw=1,crs%nums
          read(iun) crs%q(1:cri%numpw, iw)
       enddo
    else
       do iw=1,crs%numpw
          do jw=1,cri%nums
             read(iun,*) crs%q(iw,jw)
          enddo
       enddo
    endif


    close(iun)
    return

  END SUBROUTINE read_contraction_state




    SUBROUTINE do_contraction(qm,uu,cr, max_i)
!this subroutine creates the product O*U

      USE io_global,            ONLY : stdout
      USE basic_structures,     ONLY : wannier_u, q_mat

      implicit none

      TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
      TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
      TYPE(contraction) :: cr! the contraction product descriptor to be calculated
      INTEGER :: max_i !maximum number of states to be clauclates

      INTEGER :: ii,jj,kk,maxl, num_l
      INTEGER, ALLOCATABLE :: posi(:)


!free and allocates arrays
      !call free_memory_contraction(cr)

      cr%numpw=qm%numpw
      cr%nums=uu%nums
      cr%max_i=max_i



      allocate(posi(cr%nums))

      maxl=cr%nums

      write(stdout,*) 'routine do_contraction allocate dimension', cr%nums,maxl,max_i

      allocate(cr%numl(cr%numpw))
      allocate(cr%l(maxl,cr%numpw))
      allocate(cr%q(cr%numpw,maxl,max_i))

  

!do contractions
      do ii=1,cr%numpw
        posi(:)=0
        kk=0      
        cr%q(ii,:,:)=(0.d0,0.d0)
        do jj=1,qm%wp(ii)%numij
!first index
          if(posi(qm%wp(ii)%ij(1,jj))==0) then
            kk=kk+1
            posi(qm%wp(ii)%ij(1,jj))=kk
            cr%l(kk,ii)=qm%wp(ii)%ij(1,jj)
          endif
          cr%q(ii,posi(qm%wp(ii)%ij(1,jj)),1:max_i) =  cr%q(ii,posi(qm%wp(ii)%ij(1,jj)),1:max_i)+&
                             &qm%wp(ii)%o(jj)*conjg(uu%umat( 1:max_i,qm%wp(ii)%ij(2,jj),1))
!second index
          if(qm%wp(ii)%ij(1,jj)/=qm%wp(ii)%ij(2,jj)) then
            if(posi(qm%wp(ii)%ij(2,jj))==0) then
              kk=kk+1
              posi(qm%wp(ii)%ij(2,jj))=kk
              cr%l(kk,ii)=qm%wp(ii)%ij(2,jj)
            endif
            cr%q(ii,posi(qm%wp(ii)%ij(2,jj)),1:max_i) =  cr%q(ii,posi(qm%wp(ii)%ij(2,jj)),1:max_i)+&
                               &qm%wp(ii)%o(jj)*conjg(uu%umat(1:max_i, qm%wp(ii)%ij(1,jj),1))
            endif 
        enddo
        cr%numl(ii)=kk
      enddo



      deallocate(posi)
      
    END SUBROUTINE



    SUBROUTINE do_contraction_index_state(qm,uu, max_i, options)
!this subroutine creates the product O*U
!writes separately index part and states on disk
!is parallel on states


      USE io_global,            ONLY : stdout
      USE basic_structures,     ONLY : wannier_u, q_mat
      USE para_gww,             ONLY : is_my_state
      USE input_gw,             ONLY : input_options

      implicit none

      TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
      TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
      INTEGER :: max_i !maximum number of states to be clauclates
      TYPE(input_options) :: options!for calling I/O routines 

      INTEGER :: ii,jj,kk,maxl, num_l, is
      INTEGER, ALLOCATABLE :: posi(:)
      
      TYPE(contraction_index) :: cri! the contraction index descriptor to be calculated
      TYPE(contraction_state) :: crs!the contraction state to be calculated


!free and allocates arrays
      
      cri%numpw=qm%numpw
      cri%nums=uu%nums
      cri%max_i=max_i

      crs%numpw=qm%numpw
      crs%nums=uu%nums
      crs%max_i=max_i
      


      allocate(posi(cri%nums))


      maxl=cri%nums

      write(stdout,*) 'routine do_contraction_state_index allocate dimension', cri%nums,maxl,max_i
      FLUSH(stdout)
      allocate(cri%numl(cri%numpw))
      allocate(cri%l(maxl,cri%numpw))

      allocate(crs%q(cri%numpw,maxl))

      write(stdout,*) 'DO CONT INDEX 1'
      FLUSH(stdout)
!set index descriptor

!do contractions
      do ii=1,cri%numpw
         posi(:)=0
         kk=0
         do jj=1,qm%wp(ii)%numij
            !first index
            if(posi(qm%wp(ii)%ij(1,jj))==0) then
               kk=kk+1
               posi(qm%wp(ii)%ij(1,jj))=kk
               cri%l(kk,ii)=qm%wp(ii)%ij(1,jj)
            endif
!second index
            if(qm%wp(ii)%ij(1,jj)/=qm%wp(ii)%ij(2,jj)) then
               if(posi(qm%wp(ii)%ij(2,jj))==0) then
                  kk=kk+1
                  posi(qm%wp(ii)%ij(2,jj))=kk
                  cri%l(kk,ii)=qm%wp(ii)%ij(2,jj)
               endif
            endif
         enddo
         cri%numl(ii)=kk
      enddo
!write index descriptor on file
      call write_contraction_index(cri, options)


!do contractions for states
      do is=1,max_i
         if(is_my_state(is) )then
            crs%state=is
            do ii=1,crs%numpw
               posi(:)=0
               kk=0
               crs%q(ii,:)=0.d0
               do jj=1,qm%wp(ii)%numij
                  !first index
                  if(posi(qm%wp(ii)%ij(1,jj))==0) then
                     kk=kk+1
                     posi(qm%wp(ii)%ij(1,jj))=kk
                  endif
                  crs%q(ii,posi(qm%wp(ii)%ij(1,jj))) =  crs%q(ii,posi(qm%wp(ii)%ij(1,jj)))+&
                             &qm%wp(ii)%o(jj)*dble(uu%umat( is,qm%wp(ii)%ij(2,jj),1))
!second index
                  if(qm%wp(ii)%ij(1,jj)/=qm%wp(ii)%ij(2,jj)) then
                     if(posi(qm%wp(ii)%ij(2,jj))==0) then
                        kk=kk+1
                        posi(qm%wp(ii)%ij(2,jj))=kk
                     endif
                     crs%q(ii,posi(qm%wp(ii)%ij(2,jj))) =  crs%q(ii,posi(qm%wp(ii)%ij(2,jj)))+&
                               &qm%wp(ii)%o(jj)*dble(uu%umat(is, qm%wp(ii)%ij(1,jj),1))
                  endif
               enddo
            enddo
         
!writes of file
            call write_contraction_state(cri, crs, options)
         endif
      enddo

      call free_memory_contraction_index(cri)
      call free_memory_contraction_state(crs)
      deallocate(posi)

    END SUBROUTINE do_contraction_index_state








    SUBROUTINE do_contraction_pola(qm,uu,cp)
!this subroutine creates the product O*U

      USE io_global,            ONLY : stdout
      USE basic_structures,     ONLY : wannier_u, q_mat

     implicit none

      TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
      TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
      TYPE(contraction_pola) :: cp!the contraction product descriptor to be calculated


     INTEGER :: iw,jw,vv,cc,k,m,l,n,ii,jj
     INTEGER :: nums_con
     REAL(kind=DP) :: o_ii,o_jj

!free memory and set up parameters

!      call free_memory_contraction_pola(cp)
      cp%numpw=qm%numpw
      cp%nums=uu%nums
      cp%nums_occ=uu%nums_occ(1)
      nums_con=cp%nums-cp%nums_occ
      allocate(cp%ou(cp%numpw,cp%nums_occ,nums_con))

      cp%ou(:,:,:)=(0.d0,0.d0)

      do iw=1,cp%numpw
         do vv=1,cp%nums_occ
            do cc=cp%nums_occ+1,cp%nums
               do ii=1,qm%wp(iw)%numij
                  k=qm%wp(iw)%ij(1,ii)
                  m=qm%wp(iw)%ij(2,ii)
                  cp%ou(iw,vv,cc-cp%nums_occ)=cp%ou(iw,vv,cc-cp%nums_occ)+qm%wp(iw)%o(ii)*&
                      &uu%umat(vv,k,1)*uu%umat(cc,m,1)
                  if(k /= m) then
                    cp%ou(iw,vv,cc-cp%nums_occ)=cp%ou(iw,vv,cc-cp%nums_occ)+qm%wp(iw)%o(ii)*&
                      &uu%umat(vv,m,1)*uu%umat(cc,k,1)
                  endif
                enddo
             enddo
          enddo
       enddo                   

!      do iw=1,cp%numpw
!        do jw=iw,cp%numpw
!           do vv=1,cp%nums_occ
!              do cc=cp%nums_occ+1,cp%nums
!                 do ii=1,qm%wp(iw)%numij 
!                    do jj=1,qm%wp(jw)%numij
               
!                       k=qm%wp(iw)%ij(1,ii)
!                       m=qm%wp(iw)%ij(2,ii)
!                       l=qm%wp(jw)%ij(1,jj)
!                       n=qm%wp(jw)%ij(2,jj)
                       
!                       o_ii=qm%wp(iw)%o(ii)
!                       o_jj=qm%wp(jw)%o(jj)

!                       cp%q(iw,jw,vv,cc-cp%nums_occ)=cp%q(iw,jw,vv,cc-cp%nums_occ)+o_ii*o_jj*&
!  &conjg(uu%umat(k,vv))*uu%umat(l,vv)*conjg(uu%umat(m,cc))*uu%umat(n,cc)

!                       if(k /= m) then
!                         cp%q(iw,jw,vv,cc-cp%nums_occ)=cp%q(iw,jw,vv,cc-cp%nums_occ)+o_ii*o_jj*&
!  &conjg(uu%umat(m,vv))*uu%umat(l,vv)*conjg(uu%umat(k,cc))*uu%umat(n,cc)
!                       endif

!                       if(l /= n) then
!                         cp%q(iw,jw,vv,cc-cp%nums_occ)=cp%q(iw,jw,vv,cc-cp%nums_occ)+o_ii*o_jj*&
!  &conjg(uu%umat(k,vv))*uu%umat(n,vv)*conjg(uu%umat(m,cc))*uu%umat(l,cc)
!                       endif

!                       if( k /= m  .and. l /= n) then
!                         cp%q(iw,jw,vv,cc-cp%nums_occ)=cp%q(iw,jw,vv,cc-cp%nums_occ)+o_ii*o_jj*&
!  &conjg(uu%umat(m,vv))*uu%umat(n,vv)*conjg(uu%umat(k,cc))*uu%umat(l,cc)
!                       endif

!                    enddo
!                 enddo
!              enddo
!           enddo
!           cp%q(jw,iw,:,:)=conjg(cp%q(iw,jw,:,:))
!        enddo
!      enddo
    END SUBROUTINE 

    SUBROUTINE do_contraction_pola_state(qm,uu, options)
!this routine calculates contraction for all states and writes on disk
      USE io_global,            ONLY : stdout, ionode
      USE basic_structures,     ONLY : wannier_u, q_mat
      USE input_gw,             ONLY : input_options
      USE mp_world,             ONLY : mpime, nproc, world_comm
      USE mp,                   ONLY : mp_barrier
      
      implicit none

      TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
      TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
      TYPE(input_options) :: options!for i/o purpose

      INTEGER :: vv
      TYPE(contraction_pola_state) :: cps

      do vv=1,uu%nums_occ(1)
         if(mod(vv,nproc)==mpime) then         
            write(stdout,*) 'Contracting occupied state :', vv
            call do_contraction_pola_state_single(vv,qm,uu,cps)
         !if(ionode) call write_contraction_pola_state(cps, options)
            call write_contraction_pola_state(cps, options)
            call free_memory_contraction_pola_state(cps)
         endif
      enddo
      call mp_barrier( world_comm )

      return
    END SUBROUTINE do_contraction_pola_state




    SUBROUTINE do_contraction_pola_state_single(state,qm,uu,cps)
!this subroutine creates the product O*U
!for state state
!parallel

      USE io_global,            ONLY : stdout
      USE basic_structures,     ONLY : wannier_u, q_mat
      USE mp_world,             ONLY : mpime, nproc
      USE mp,                   ONLY : mp_sum

     implicit none

      INTEGER :: state!state for which the contraction will be calculated
      TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
      TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
      TYPE(contraction_pola_state) :: cps!the contraction product descriptor to be calculated


     INTEGER :: iw,jw,cc,k,m,l,n,ii,jj
     INTEGER :: nums_con
     REAL(kind=DP) :: o_ii,o_jj

!free memory and set up parameters

!      call free_memory_contraction_pola(cp)
      cps%numpw=qm%numpw
      cps%nums=uu%nums
      cps%nums_occ=uu%nums_occ(1)
      cps%state=state
      nums_con=cps%nums-cps%nums_occ
      allocate(cps%ou(nums_con,cps%numpw))

      cps%ou(:,:)=0.d0
 

      do iw=1,cps%numpw
!         do cc=cps%nums_occ+1,cps%nums
!            if(mod(cc,nproc)==mpime) then
         do ii=1,qm%wp(iw)%numij
            k=qm%wp(iw)%ij(1,ii)
            m=qm%wp(iw)%ij(2,ii)
            do cc=cps%nums_occ+1,cps%nums
               cps%ou(cc-cps%nums_occ,iw)=cps%ou(cc-cps%nums_occ,iw)+qm%wp(iw)%o(ii)*&
                    &dble(uu%umat(cps%state,k,1))*dble(uu%umat(cc,m,1))
               if(k /= m) then
                  cps%ou(cc-cps%nums_occ,iw)=cps%ou(cc-cps%nums_occ,iw)+qm%wp(iw)%o(ii)*&
                       &dble(uu%umat(cps%state,m,1))*dble(uu%umat(cc,k,1))
               endif
            enddo
!            endif
         enddo
!         call mp_sum(cps%ou(:,iw))
      enddo
       

    END SUBROUTINE do_contraction_pola_state_single

    SUBROUTINE read_contraction_pola_state(cps, options)
!this subroutine writes the contracted pola products on disk
!in parallel case only ionode writes

    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_pola_state)             :: cps!the contraction pola state to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun
    CHARACTER(5) :: nfile



    write(nfile,'(5i1)') &
         & cps%state/10000,mod(cps%state,10000)/1000,mod(cps%state,1000)/100,mod(cps%state,100)/10,mod(cps%state,10)
    iun = find_free_unit()
    if(.not. options%debug) then
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_pola'// nfile, status='old',form='unformatted')
    else
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_pola'// nfile, status='old',form='formatted')
    endif

    if(.not.options%debug) then
       read(iun) cps%numpw
       read(iun) cps%nums
       read(iun) cps%nums_occ
       read(iun) cps%state
    else
       read(iun,*) cps%numpw
       read(iun,*) cps%nums
       read(iun,*) cps%nums_occ
       read(iun,*) cps%state
    endif
    allocate(cps%ou(cps%nums-cps%nums_occ,cps%numpw))

    if(.not.options%debug) then
       do iw=1,cps%numpw
          read(iun) cps%ou(1:(cps%nums-cps%nums_occ),iw)
       enddo
    else
       do iw=1,cps%numpw
          do jw=1,cps%nums-cps%nums_occ
             read(iun,*) cps%ou(jw,iw)
          enddo
       enddo
    endif
    close(iun)
    return
  END SUBROUTINE read_contraction_pola_state

  SUBROUTINE write_contraction_pola_state(cps, options)
!this subroutine writes the contracted pola products on disk
!in parallel case only ionode writes

    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_pola_state)             :: cps!the contraction pola state to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun
    CHARACTER(5) :: nfile



    write(nfile,'(5i1)') &
         & cps%state/10000,mod(cps%state,10000)/1000,mod(cps%state,1000)/100,mod(cps%state,100)/10,mod(cps%state,10)
    iun = find_free_unit()
    if(.not. options%debug) then
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_pola'// nfile, status='unknown',form='unformatted')
    else
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction_pola'// nfile, status='unknown',form='formatted')
    endif

    if(.not.options%debug) then
       write(iun) cps%numpw
       write(iun) cps%nums
       write(iun) cps%nums_occ
       write(iun) cps%state
       do iw=1,cps%numpw
          write(iun) cps%ou(1:(cps%nums-cps%nums_occ),iw)
       enddo
    else
       write(iun,*) cps%numpw
       write(iun,*) cps%nums
       write(iun,*) cps%nums_occ
       write(iun,*) cps%state
       do iw=1,cps%numpw
          do jw=1,cps%nums-cps%nums_occ
             write(iun,*) cps%ou(jw,iw)
          enddo
       enddo
    endif
    close(iun)
    return
  END SUBROUTINE write_contraction_pola_state

   SUBROUTINE read_contraction_state_central(cri,crs, options)
!this subroutine writes the contracted products on disk
!in parallel case only ionode writes

    USE input_gw,           ONLY : input_options
    USE io_global,          ONLY : ionode, ionode_id
    USE mp,                 ONLY : mp_bcast
    USE mp_world,           ONLY : world_comm
    USE io_files,           ONLY : prefix, tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(contraction_index), INTENT(in) :: cri!the contraction index descriptor
    TYPE(contraction_state)             :: crs!the contraction state to be read from file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun
    CHARACTER(5) :: nfile
    INTEGER :: maxl


    if(ionode) then
       write(nfile,'(5i1)') &
         & crs%state/10000,mod(crs%state,10000)/1000,mod(crs%state,1000)/100,mod(crs%state,100)/10,mod(crs%state,10)
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='old',form='unformatted')
       else
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'contraction'// nfile, status='old',form='formatted')
       endif

       if(.not.options%debug) then
          read(iun) crs%numpw
          read(iun) crs%nums
          read(iun) crs%max_i
          read(iun) crs%state
       else
          read(iun,*) crs%numpw
          read(iun,*) crs%nums
          read(iun,*) crs%max_i
          read(iun,*) crs%state
       endif
    endif
    call mp_bcast(crs%numpw, ionode_id,world_comm)
    call mp_bcast(crs%nums, ionode_id,world_comm)
    call mp_bcast( crs%max_i, ionode_id,world_comm)
    call mp_bcast(crs%state, ionode_id,world_comm)

    maxl=crs%nums
    allocate(crs%q(crs%numpw,maxl))

    if(ionode) then
       if(.not.options%debug) then
          do iw=1,crs%numpw
             read(iun) crs%q(iw,1:cri%numl(iw))
          enddo
       else
          do iw=1,crs%numpw
             do jw=1,cri%numl(iw)
                read(iun,*) crs%q(iw,jw)
             enddo
          enddo
       endif
       close(iun)
    endif
    do iw=1,crs%numpw
       call mp_bcast( crs%q(iw,1:cri%numl(iw)), ionode_id,world_comm)
    enddo
    return

  END SUBROUTINE read_contraction_state_central





   END MODULE compact_product


