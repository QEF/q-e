!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

  MODULE lanczos
!this module describes the structures for the calculation
!of the polarization and of the self-energy through an 
!lanczos chain style                                                                   

    USE kinds, ONLY : DP

    TYPE compact_q_lanczos
!this structure describes the "compact" term:
! Q^v_in=\sum U_{vv'}V^v'_{i,l}T^v'_{l,n}

       INTEGER :: ii!corresponding KS state 
       INTEGER :: numpw!dimension of polarization basis
       INTEGER :: numt!dimension of the basis {t_n} 
       REAL(kind=DP), POINTER, DIMENSION(:,:) :: qlm!matrix Q(numpw,numt)
    END TYPE compact_q_lanczos

    TYPE lanczos_matrix
!this structure describes the (H-i\alpha)^-1 matrix
       INTEGER :: iw!corresponding imaginary frequency
       INTEGER :: numt!dimension of the basis {t_n}
       COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: e_mat
    END TYPE lanczos_matrix

  CONTAINS

    SUBROUTINE initialize_compact_q_lanczos(cql)
!this subroutine initializes compact_q_lanczos   

      implicit none
      TYPE(compact_q_lanczos) :: cql
      nullify(cql%qlm)
      return
    END SUBROUTINE initialize_compact_q_lanczos

    SUBROUTINE free_memory_compact_q_lanczos(cql)
!this subroutine initializes compact_q_lanczos   

      implicit none
      TYPE(compact_q_lanczos) :: cql
      if(associated(cql%qlm)) deallocate(cql%qlm)
      nullify(cql%qlm)
      return
    END SUBROUTINE free_memory_compact_q_lanczos

    SUBROUTINE initialize_lanczos_matrix(lm)
!this subroutine initializes compact_q_lanczos   

      implicit none
      TYPE(lanczos_matrix) :: lm
      nullify(lm%e_mat)
      return
    END SUBROUTINE initialize_lanczos_matrix

    SUBROUTINE free_memory_lanczos_matrix(lm)
!this subroutine initializes compact_q_lanczos

      implicit none
      TYPE(lanczos_matrix) :: lm
      if(associated(lm%e_mat)) deallocate(lm%e_mat)
      nullify(lm%e_mat)
      return
    END SUBROUTINE free_memory_lanczos_matrix

     SUBROUTINE write_compact_q_lanczos(cql)
!this subroutine writes the compact_q_lanczos function on disk  
!the file name is taken from the label                                                                     

    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(compact_q_lanczos) :: cql!the compact_q_lanczos function to be written   

    INTEGER :: iunq, ii
    CHARACTER(5) :: nfile

    write(nfile,'(5i1)') &
        & cql%ii/10000,mod(cql%ii,10000)/1000,mod(cql%ii,1000)/100,mod(cql%ii,100)/10,mod(cql%ii,10)
    iunq = find_free_unit()    
    open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'-'//'q_lanczos.'// nfile, status='unknown',form='unformatted')
    
    write(iunq) cql%ii
    write(iunq) cql%numpw
    write(iunq) cql%numt
    do ii=1,cql%numt
       write(iunq) cql%qlm(1:cql%numpw,ii)
    enddo

    close(iunq)
    return
  END SUBROUTINE write_compact_q_lanczos

     SUBROUTINE read_compact_q_lanczos(cql, iv)
!this subroutine reads the compact_q_lanczos function from disk     

    USE io_files,             ONLY : prefix,tmp_dir
    USE mp,                   ONLY : mp_barrier,mp_bcast, mp_sum
    USE mp_world,             ONLY : world_comm
    USE io_global,            ONLY : ionode,ionode_id
    implicit none

    INTEGER, EXTERNAL :: find_free_unit

    TYPE(compact_q_lanczos) :: cql!the compact_q_lanczos function to be read 
    INTEGER, INTENT(in) :: iv!the index of the file to be read

    INTEGER :: iunq, ii
    CHARACTER(5) :: nfile

    call free_memory_compact_q_lanczos(cql)

    cql%ii=iv
    write(nfile,'(5i1)') &
        & cql%ii/10000,mod(cql%ii,10000)/1000,mod(cql%ii,1000)/100,mod(cql%ii,100)/10,mod(cql%ii,10)
    if(ionode) then
       iunq = find_free_unit()
       open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'-'//'q_lanczos.'// nfile, status='old',form='unformatted')       
       read(iunq) cql%ii
       read(iunq) cql%numpw
       read(iunq) cql%numt
    endif
    call mp_bcast(cql%ii,ionode_id,world_comm)
    call mp_bcast(cql%numpw,ionode_id,world_comm)
    call mp_bcast(cql%numt,ionode_id,world_comm)

    allocate(cql%qlm(cql%numpw,cql%numt))

    do ii=1,cql%numt
       if(ionode) then
          read(iunq) cql%qlm(1:cql%numpw,ii)
       else
          cql%qlm(1:cql%numpw,ii)=0.d0
       endif
       !call mp_barrier
       !call mp_bcast(cql%qlm(1:cql%numpw,ii),ionode_id,world_comm)
       !call mp_sum(cql%qlm(1:cql%numpw,ii))
   enddo
   call mp_bcast(cql%qlm(:,:), ionode_id,world_comm)
   
    if(ionode) close(iunq)

    return

  END SUBROUTINE read_compact_q_lanczos


  SUBROUTINE do_compact_q_lanczos(vtl,ttl,cql,alpha)
!this subroutines performs the calculation:
! Q^v'_in= Q^v'_in +alpha*V^v'_{i,l}T^v'_{l,n}

    USE kinds, ONLY : DP
    USE basic_structures, ONLY : tt_mat_lanczos, vt_mat_lanczos
    USE io_global, ONLY : stdout, ionode, ionode_id

    implicit none

    TYPE(vt_mat_lanczos), INTENT(in) :: vtl!V matrix
    TYPE(tt_mat_lanczos), INTENT(in) :: ttl!T matrix 
    TYPE(compact_q_lanczos), INTENT(out) :: cql!Q matrix to be calculated
    REAL(kind=DP), INTENT(in) :: alpha!constant alpha

    INTEGER il,it,ip




    if(ttl%ii /= vtl%ii) then
       write(stdout,*) 'Routine do_compact_q_lanczos: state v not equal'
       FLUSH(stdout)
       stop
    else
       cql%ii=ttl%ii
    endif
    cql%numpw=vtl%numpw
    cql%numt=ttl%numt

   

    call dgemm('N','T',cql%numpw,cql%numt,vtl%numl,alpha,vtl%vt_mat,vtl%numpw,ttl%tt_mat,ttl%numt,1.d0,cql%qlm,cql%numpw)

!   cql%qlm(:,:)=0.d0
!    do ip=1,cql%numpw
!       do it=1,cql%numt
!          do il=1,vtl%numl
!             cql%qlm(ip,it)=cql%qlm(ip,it)+vtl%vt_mat(ip,il)*ttl%tt_mat(it,il)
!          enddo
!       enddo
!    enddo

    return
  END SUBROUTINE do_compact_q_lanczos


     SUBROUTINE write_lanczos_matrix(lm)
!this subroutine writes the lanczos matrix  on disk  
!the file name is taken from the label                                                                     

    USE io_files,             ONLY : prefix,tmp_dir
    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(lanczos_matrix) :: lm!the matrix to be written   

    INTEGER :: iunm, ii
    CHARACTER(5) :: nfile

    if(lm%iw >= 0) then
       write(nfile,'(5i1)') &
         & lm%iw/10000,mod(lm%iw,10000)/1000,mod(lm%iw,1000)/100,mod(lm%iw,100)/10,mod(lm%iw,10)
       iunm = find_free_unit()    
       open( unit=iunm, file=trim(tmp_dir)//trim(prefix)//'-'//'emat_lanczos.'// nfile, status='unknown',form='unformatted')
    else
       write(nfile,'(5i1)') &
         & -lm%iw/10000,mod(-lm%iw,10000)/1000,mod(-lm%iw,1000)/100,mod(-lm%iw,100)/10,mod(-lm%iw,10)
       iunm = find_free_unit()
       open( unit=iunm, file=trim(tmp_dir)//trim(prefix)//'-'//'emat_lanczos.-'// nfile, status='unknown',form='unformatted')
    endif
    write(iunm) lm%iw
    write(iunm) lm%numt
    do ii=1,lm%numt
       write(iunm) lm%e_mat(1:lm%numt,ii)
    enddo

    close(iunm)
    return
  end SUBROUTINE write_lanczos_matrix

    SUBROUTINE read_lanczos_matrix(lm,iw)
!this subroutine reads the lanczos matrix  from disk
!the file name is taken from the label
!it does not allocate the matrix  
    USE io_files,             ONLY : prefix,tmp_dir
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(lanczos_matrix) :: lm!the matrix to be read
    INTEGER :: iw!index of matrix to be read

    INTEGER :: iunm, ii
    CHARACTER(5) :: nfile

    lm%iw=iw
    
    if(lm%iw >= 0) then
       write(nfile,'(5i1)') &
            & lm%iw/10000,mod(lm%iw,10000)/1000,mod(lm%iw,1000)/100,mod(lm%iw,100)/10,mod(lm%iw,10)
       iunm = find_free_unit()
       open( unit=iunm, file=trim(tmp_dir)//trim(prefix)//'-'//'emat_lanczos.'// nfile, status='old',form='unformatted')
    else
       write(nfile,'(5i1)') &
            & -lm%iw/10000,mod(-lm%iw,10000)/1000,mod(-lm%iw,1000)/100,mod(-lm%iw,100)/10,mod(-lm%iw,10)
       iunm = find_free_unit()
       open( unit=iunm, file=trim(tmp_dir)//trim(prefix)//'-'//'emat_lanczos.-'// nfile, status='unknown',form='unformatted')
    endif
    read(iunm) lm%iw
    read(iunm) lm%numt
    do ii=1,lm%numt
       read(iunm) lm%e_mat(1:lm%numt,ii)
    enddo

    close(iunm)
    return


  END SUBROUTINE read_lanczos_matrix




    
  END MODULE lanczos
