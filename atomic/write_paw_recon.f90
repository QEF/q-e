!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
      subroutine write_paw_recon
!--------------------------------------------------------------
use ld1inc
implicit none

      integer :: i, j, n, m, l, ios, iae, isign
      real(DP) ::  wmax

      open(unit=51,file=file_recon,status='unknown', err=1111, &
           iostat=ios,form='formatted')
1111     call errore('write_result','opening 51',abs(ios))

      iae=nwf

      write (51,*) '<PP_PAW>'
      write (51,*) nwfts
      write (51,*) '</PP_PAW>'
      do i=nwfts,1,-1
         write (51,*) '<PP_REC>'
         write (51,*) '<PP_kkbeta>'
         write (51,*) mesh
         write (51,*) '</PP_kkbeta>'
         write (51,*) '<PP_L>'
         write (51,*) llts(i)
         write (51,*) '</PP_L>'
         write (51,*) '<PP_REC_AE>'

! Check sign of ae-wfct (max should be >0)

         wmax=0.0_dp
         do n=1,mesh
            if(abs(psi(n,1,iae)).gt.wmax.and.r(n).lt.4.0_dp)then
               wmax=abs(psi(n,1,iae))
               if(psi(n,1,iae).lt.0.0_dp)then
                  isign=-1
               else
                  isign=+1
               endif
            endif
         enddo

         write (51,'(1p4e19.11)') (isign*psi(n,1,iae) , n=1,mesh )
         write (51,*) '</PP_REC_AE>'
         write (51,*) '<PP_REC_PS>'

! check sign of pseudo wfct (max should be >0)

         wmax=0.0_dp
         do n=1,mesh
            if(abs(phits(n,i)).gt.wmax.and.r(n).lt.4.0_dp)then
               wmax=abs(phits(n,i))
               if(phits(n,i).lt.0.0_dp)then
                  isign=-1
               else
                  isign=+1
               endif
            endif
         enddo


         write (51,'(1p4e19.11)') (isign*phits(n,i) , n=1,mesh )
         write (51,*) '</PP_REC_PS>'
         write (51,*) '</PP_REC>'
         iae=iae-1
      enddo
      close(51)


      return
    end subroutine write_paw_recon
