! Copyright (C) 2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)


subroutine wannier_check()
  use io_global,  only : stdout
  use kinds, only : DP 
  use klist, only : nks, nkstot
  use ions_base, only : nat, ityp, atm,tau
  use wvfct, only: nbnd
  USE basis, only: natomwfc
  use wannier_new, only: nwan, wan_in, use_energy_int
  use lsda_mod, only: nspin
  USE control_flags,    ONLY : gamma_only
  USE uspp_param,       ONLY : upf 

  implicit none
  integer :: nwfc, lmax_wfc, na, nt, n, l, m, i, iwan, ispin

  ! number of k points in this pool .ne. total number of k points
  if (nks.ne.nkstot) call errore ('wannier_check', 'not implemented 1', 1)
  
  if ( gamma_only ) call errore ('wannier_check', 'gamma_only calculation not implemented', 1) 

  ! here we will write to stdout source of wannier functions (atomic functions from which wannier are generated) 
  do ispin=1, nspin
     !
     write(stdout,'(5x,a4,i2)') 'Spin',ispin
     do iwan=1,nwan
        write(stdout,'(7x,"Wannier #",i3," centered on atom ",a3," (position ",3f8.5," )")') &
             iwan, atm(ityp(wan_in(iwan,ispin)%iatom)), (tau(i,wan_in(iwan,ispin)%iatom),i=1,3)
        
        if( use_energy_int) then
           write(stdout,'(9x,"Bands for generation: from",f6.3," to",f6.3)') &
                wan_in(iwan,ispin)%bands_from,wan_in(iwan,ispin)%bands_to
        else 
           write(stdout,'(9x,"Bands for generation: from",i4," to",i4)') &
                INT(wan_in(iwan,ispin)%bands_from),INT(wan_in(iwan,ispin)%bands_to)
        end if
        
        write(stdout,'(9x,a31)') 'Trial wavefunction ingredients:'
        do i=1,wan_in(iwan,ispin)%ning
           nwfc=0 
           lmax_wfc = 0 
           write(stdout,'(10x,f12.10," of l=",i1,", m=",i1)') &
                wan_in(iwan,ispin)%ing(i)%c, wan_in(iwan,ispin)%ing(i)%l, wan_in(iwan,ispin)%ing(i)%m
           
           ! now we shoud associate every ingridient of trial wavefunction with atomic orbital
           ! it will be done only once - for future using in wannier_proj
           DO na = 1, nat
              nt = ityp (na)
              DO n = 1, upf(nt)%nwfc
                 IF (upf(nt)%oc (n) >= 0.d0) THEN
                    l = upf(nt)%lchi (n)
                    lmax_wfc = max (lmax_wfc, l ) 
                    DO m=1, 2*l+1          
                       nwfc=nwfc+1 
                       ! the most important part
                       if ( &
                            (na == wan_in(iwan,ispin)%iatom) .AND. &
                            (l == wan_in(iwan,ispin)%ing(i)%l) .AND. &
                            (m == wan_in(iwan,ispin)%ing(i)%m) ) &
                            wan_in(iwan,ispin)%ing(i)%iatomwfc = nwfc
                    enddo
                 endif
              enddo
           enddo
           
        end do ! ingredients
     end do ! iwannier
  end do !ispin
  
  ! do iwan=1,nwan
  !   write(stdout,'(7x,"Wannier #",i3," atomic wavefunction", i3)') iwan, wan_in(iwan,1)%ing(1)%iatomwfc
 ! end do ! iwannier
     
  if (lmax_wfc > 3) call errore ('wannier_check', 'l > 3 not yet implemented', 1) 
  if (nwfc /= natomwfc) call errore ('wannier_check', 'wrong # of atomic wfcs?', 1)
  if (nwan > nbnd ) call errore( 'wannier_check','too few bands', nwan-nbnd)

  return

end subroutine wannier_check
