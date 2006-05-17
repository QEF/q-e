!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!     
!---------------------------------------------------------------------
      subroutine read_newpseudo (ios)
!---------------------------------------------------------------------
!
!     This routine reads from input the quantities which defines 
!     a multiprojector pseudopotential. It can be in the
!     Vanderbilt form or in the norm-conserving form
!
use ld1inc
use funct, only: set_dft_from_indices
   implicit none

      integer :: &
             nb,mb, &  ! counters on beta functions
             n,     &  ! counter on mesh points
             ir,    &  ! counters on mesh points
             ios,   &  ! I/O control: ios /= 0 means error
             iunps    ! the unit with the pseudopotential
      integer :: iexch, icorr, igcx, igcc

      logical :: reldum

      if (file_pseudo.eq.' ') return

      iunps=29
      open(unit=iunps,file=file_pseudo,status='unknown', &
         &  form='formatted', err=50, iostat=ios)
50    call errore('read_newpseudo','opening file_pseudo',abs(ios))

      read( iunps, '(a75)', err=100, iostat=ios ) title

      read( iunps, '(i5)',err=100, iostat=ios ) pseudotype
      if (pseudotype /= 2 .and. pseudotype /= 3) &
         call errore('read_newpseudo','pseudotype is wrong',1)

      read( iunps, '(2l5)',err=100, iostat=ios ) reldum, nlcc
      if (reldum.and.rel.eq.0) call errore('read_newpseudo', &
   &    'relativistic pseudopotential and non relativistic calculation',-1)
      if (.not.reldum.and.rel.gt.0) call errore('read_newpseudo', &
   &    'non relativistic pseudopotential and relativistic calculation',-1)
         
      read( iunps, '(4i5)',err=100, iostat=ios ) iexch, icorr, igcx, igcc
      call set_dft_from_indices(iexch, icorr, igcx, igcc)
      
      read( iunps, '(2e17.11,i5)') zval, etots, lmax

      read( iunps, '(4e17.11,i5)',err=100, iostat=ios ) &
                                   xmin,rmax,zmesh,dx,mesh

      call do_mesh(rmax,zmesh,xmin,dx,0,ndm,mesh,r,r2,rab,sqr)

      read( iunps, '(2i5)', err=100, iostat=ios ) nwfs, nbeta
      read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                                    ( rcut(nb), nb=1,nwfs )
      read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                                    ( rcutus(nb), nb=1,nwfs )

      do nb=1,nwfs
         read(iunps,'(a2,2i3,f6.2)',err=100,iostat=ios) &
                          els(nb), nns(nb), lls(nb), ocs(nb)
      enddo
      do nb=1,nbeta
         read ( iunps, '(i6)',err=100, iostat=ios ) ikk(nb)
         read ( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
                            ( betas(ir,nb), ir=1,ikk(nb))
         do ir=ikk(nb)+1,mesh
            betas(ir,nb)=0.0_dp
         enddo
         do mb=1,nb
            read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                  bmat(nb,mb)
            bmat(mb,nb)=bmat(nb,mb)
            if (pseudotype.eq.3) then
              read(iunps,'(1p4e19.11)',err=100,iostat=ios) & 
                  qq(nb,mb)
              qq(mb,nb)=qq(nb,mb)
              read(iunps,'(1p4e19.11)',err=100,iostat=ios)   &  
                 (qvan(n,nb,mb),n=1,mesh)
              do n=1,mesh
                 qvan(n,mb,nb)=qvan(n,nb,mb)
              enddo
            else
              qq(nb,mb)=0.0_dp
              qq(mb,nb)=0.0_dp
              do n=1,mesh
                 qvan(n,mb,nb)=0.0_dp
                 qvan(n,nb,mb)=0.0_dp
              enddo
            endif
         enddo
      enddo
!
!   reads the local potential 
!
      read( iunps, '(1p4e19.11)',err=100, iostat=ios ) rcloc, &
                             ( vpsloc(ir), ir=1,mesh )
!
!     reads the atomic charge
!
      read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                               ( rhos(ir,1), ir=1,mesh )
!
!  if present reads the core charge
!
      if ( nlcc ) then 
         read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                               ( rhoc(ir), ir=1,mesh )
      else
         rhoc = 0.0_dp
      endif
!
!    read the pseudo wavefunctions of the atom
!      
      read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
                   ((phis(ir,nb),ir=1,mesh),nb=1,nwfs)
100   continue
      !
      ! do not stop with error message here: return error code instead
      !
      !call errore('read_newpseudo','Reading pseudo file', &
      !                                         abs(ios))
      close(iunps)

      return
      end subroutine read_newpseudo
